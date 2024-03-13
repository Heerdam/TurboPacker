#pragma once

#include <MQT2.hpp>

#include <atomic>
#include <condition_variable>
#include <thread>
#include <exception>
#include <concepts>
#include <type_traits>
#include <functional>
#include <random>
#include <syncstream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/ext/scalar_constants.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>

namespace TP {

    namespace Detail {
        template<class, typename, int32_t, typename> 
        class Promise;
    }

    template<class, template<typename> class, class, int32_t, class>
    struct Config;

    //----------------------

    template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
    [[nodiscard]] Detail::Promise<T, H_T, BS, HA> solve(const Config<T, CF, H_T, BS, HA>& _config);

    namespace Detail {

        [[nodiscard]] inline int32_t get_power_of_2(const int32_t _target, const int32_t _base) noexcept {
            int32_t t = _base;
            while (t < _target)
                t  = t << 1;
            return t;
        }//get_power_2

        template<class T>
        struct FBox {
            glm::vec<3, T> min_;
            glm::vec<3, T> max_;
            [[nodiscard]] glm::vec<3, T> GetExtent() const { return (max_ - min_) * T(0.5); }
            [[nodiscard]] glm::vec<3, T> GetCenter() const { return min_ + GetExtent(); }
        };//FBox

        //-------------------------

        class TaskGraph {

            std::vector<std::thread> t_;
            std::queue<std::function<void()>> tasks_;
            std::condition_variable cv_, cv_finished_;
            std::mutex m_;
            bool run_ = false;
            int active_tasks_ = 0;

            void run_impl() {
                while(true){
                    std::unique_lock<std::mutex> lock(m_);
                    cv_.wait(lock, [&]{ return !run_ || !tasks_.empty(); });
                    if(!run_) break;
                    if(tasks_.empty()) continue;
                    const auto tt = std::move(tasks_.front());
                    tasks_.pop();
                    ++active_tasks_;
                    lock.unlock();
                    tt();
                    lock.lock();     
                    --active_tasks_;
                    if (tasks_.empty() && active_tasks_ == 0){
                        lock.unlock();
                        cv_finished_.notify_all();
                    }
                }
            }

        public:

            ~TaskGraph() {
                {
                    std::lock_guard<std::mutex> lock(m_);
                    run_ = false;
                }
                cv_.notify_all();
                for(size_t i = 0; i < t_.size(); ++i)
                    t_[i].join();
            }

            TaskGraph() = default;
            TaskGraph(const int32_t _num_threads) : run_(true){
                t_.resize(_num_threads);
                for(size_t i = 0; i < t_.size(); ++i){
                    t_[i] = std::thread(&TaskGraph::run_impl, this);
                }
            }

            template<class FUNC>
            void dispatch(FUNC _task){
                {
                    std::lock_guard<std::mutex> lock(m_);
                    tasks_.push(std::function<void()>(std::move(_task)));
                }
                cv_.notify_one();
            }

            void wait() {
                std::unique_lock<std::mutex> lock(m_);
                cv_finished_.wait(lock, [&] { return tasks_.empty() && active_tasks_ == 0; });
            }

        };//TaskGraph

        //-------------------------

        enum class EAxisPerm : uint8_t {
            Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
            Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
            X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,
        };//EAxisPerm

        template<class T>
        [[nodiscard]] glm::mat<4, 4, T> make_transform(
            const EAxisPerm _perm,
            const FBox<T>& _target,
            const glm::vec<3, T>& _pivot_offset,
            const glm::vec<3, T>& _ext
        );//make_transform

        //-------------------------

        template<class T, typename H_T, int32_t BS, typename HA>
        struct SolverContext {

            std::vector<H_T> map_;
            std::unique_ptr<MQT2::MedianQuadTree<H_T, BS, HA>> tree_;

            std::thread mt_;
            TaskGraph tg_;

            std::vector<glm::mat<4, 4, T>> data_;
            std::mutex m_data;

            int32_t N_;
            T tot_vol_;

            std::chrono::high_resolution_clock::time_point time_start_;
            std::atomic<bool> isDone_ = true;
            std::atomic<double> vol_ = 0.;
            std::atomic<int32_t> bcc_ = 0;
            std::atomic<int32_t> mcc_ = 0;

            SolverContext(const int32_t _num_threads) : tg_(TaskGraph(_num_threads)) {}
            SolverContext(const SolverContext&) = delete;
            SolverContext(SolverContext&&) = delete;
            ~SolverContext() {
                mt_.join();
            }

        };//SolverContext

        //-------------------------

        template<class T, typename H_T, int32_t BS, typename HA>
        class Promise {

            std::unique_ptr<SolverContext<T, H_T, BS, HA>> context_;      
            //-------------
            Promise() = default;
            Promise(Promise&&) = default;
            Promise(const Promise&) = delete;
            Promise& operator=(Promise&&) = default;
            Promise& operator=(const Promise&) = delete;
            //-------------
            template<typename T_, template<typename> class CF_, typename H_T_ , int32_t BS_, typename HA_>
            friend Promise<T_, H_T_, BS_, HA_> TP::solve(const ::TP::Config<T_, CF_, H_T_, BS_, HA_>& _conf);
        
        public:
            //locks the thread until isDone() == true
            void wait();

            //stops the packing and shutdowns all threads safely
            void stop();

            //NOT thread safe! only call when isDone() returns true!
            [[nodiscard]] const std::vector<glm::mat<4, 4, T>>& data() const;

            //NOT thread safe! only call when isDone() returns true!
            [[nodiscard]] std::vector<glm::mat<4, 4, T>>& data();
            
            //thread safe!
            [[nodiscard]] bool isDone() const;

            //thread safe! returns the current pack density [0, 1]
            [[nodiscard]] T getPackDensity() const;

            //thread safe!
            [[nodiscard]] int32_t getBoxCount() const;

            //thread safe!
            [[nodiscard]] int32_t getMissedCount() const;

            //thread safe! seconds since start
            [[nodiscard]] double getTime() const;

        };//Promise

        //-------------------------

        enum class BoxGenerationType : uint8_t {
            RANDOM, LIST
        };//PackType

        //-------------------------

        struct BoxEntry {
            glm::vec3 Size;
            int32_t Count;
        };//BoxEntry

        struct BoxList {
            std::vector<BoxEntry> List;
            bool ShuffleBoxes = false;
        };//BoxList

        //-------------------------
        template<class T>
        struct Result {
            bool isRandomBox;
            T weight;
            int32_t n0, n1, h;
            int32_t l, b_l, b_m, b_h;
            int32_t set_index;
            glm::vec<3, T> ext;
            glm::vec<3, T> ext_org;
            EAxisPerm perm;
        };//Result

    }//Detail

    namespace CostFunction {

        template<class T>
        struct CF_Basic {
            static T eval(const Detail::Result<T>& _res) { return 0.; };
        };//CF_Basic

        //------------------------------

        template<typename T, template<typename> class CF>
        concept isValidCF = requires(Detail::Result<T> a) { { CF<T>::eval(a) } -> std::same_as<T>; };

    };//CostFunction

    template<class T, template<typename> class COSTFUNCTION, class HEIGHTMAP_T = int16_t, int32_t BUCKET_SIZE = 15, class HEIGHTMAP_ALLOCATOR = std::allocator<HEIGHTMAP_T>>
    struct Config {
        using T_ = T;
        using HEIGHTMAP_T_ = HEIGHTMAP_T;
        constexpr static int32_t BUCKET_SIZE_ = BUCKET_SIZE;
        using HEIGHTMAP_ALLOCATOR_ = HEIGHTMAP_ALLOCATOR;
        using COSTFUNCTION_ = COSTFUNCTION<T>;
        //----------------------

        //if the solver should use worker threads [default: true]
        bool MultiThreading = true;

        //number of worker threads [default: 4]
        int32_t NumThreads = 4;

        //if the solver should use a random seed [default: true]
        bool UseRandomSeed_ = true;

        //the seed [default: 1234567890]
        uint64_t Seed = 1234567890;

        //the bounds of the domain [default: none]
        glm::vec<2, T> Bounds;

        //the max height [default: 150]
        int32_t Height = 150;

        //if random boxes or a pre-defined list of boxes should be used [default: Random]
        Detail::BoxGenerationType BoxType = Detail::BoxGenerationType::RANDOM;

        //if boxes are allowed to overlap [default: false]
        bool AllowOverlap = false;

        //how many misses in a row before terminating [default: 25]
        int32_t EmptryTries = 25;

        //how many misses in total before terminating [default: 50]
        int32_t MaxEmptryTries = 50;

        //how large the test set is. The same as how far the solver can look ahead [default: 1]
        int32_t SetSize = 1;

        //--------------------------

        //if the random boxes should be cubes [default: true]
        bool CubeRandomBoxes = true;

        //the minimum volume of the random boxes [default: 250]
        double MinBoxVolume = 250.;

        //the maximum volume of the random boxes [default: 1500]
        double MaxBoxVolume = 1500.;

        //--------------------------

        //the list fo predefined boxes
        Detail::BoxList BoxList;

        //if it should enforce emptry tries. if false the solver runs until the list is empty [default: false]
        bool EnforceMisses = false;

    };//Config

    namespace Detail {

        template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
        void run_impl(
            Config<T, CF, H_T, BS, HA> _conf,
            SolverContext<T, H_T, BS, HA>* _cont
        );//run_impl

        template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
        [[nodiscard]] std::vector<Result<T>> overlap_impl(
            const Config<T, CF, H_T, BS, HA>& _conf,
            SolverContext<T, H_T, BS, HA>* _cont,
            const int32_t _ext0, 
            const int32_t _ext1, 
            const int32_t _h
        );//overlap_impl

        template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
        void dispatch_impl(
            const Config<T, CF, H_T, BS, HA>& _conf,
            SolverContext<T, H_T, BS, HA>* _cont,
            std::vector<Result<T>>& _res,
            std::mutex& _m,
            std::atomic<double>& _minc,
            const int32_t _set_index,
            const int32_t _n0, const int32_t _n1,
            const glm::vec<3, T>& _ext_org,
            const int32_t _h, EAxisPerm _perm
        );//dispatch_impl       

    }//Detail

}//TP

//------------------------------

template<class T, typename H_T, int32_t BS, typename HA>
void TP::Detail::Promise<T, H_T, BS, HA>::wait() {
    using namespace std::chrono_literals;
    while(!isDone())
        std::this_thread::sleep_for(0.25s);
}//TP::Detail::Promise::wait

template<class T, typename H_T, int32_t BS, typename HA>
void TP::Detail::Promise<T, H_T, BS, HA>::stop() {
    assert(context_);
    context_->isPacking_ = false;
}//TP::Detail::Promise::wait

template<class T, typename H_T, int32_t BS, typename HA>
const std::vector<glm::mat<4, 4, T>>& TP::Detail::Promise<T, H_T, BS, HA>::data() const {
    assert(context_);
    std::lock_guard<std::mutex> lock(context_->m_data);
    return context_->data_;
}//TP::Detail::Promise::data

template<class T, typename H_T, int32_t BS, typename HA>
std::vector<glm::mat<4, 4, T>>& TP::Detail::Promise<T, H_T, BS, HA>::data() {
    assert(context_);
    std::lock_guard<std::mutex> lock(context_->m_data);
    return context_->data_;
}//TP::Detail::Promise::data

template<class T, typename H_T, int32_t BS, typename HA>
bool TP::Detail::Promise<T, H_T, BS, HA>::isDone() const {
    assert(context_);
    return context_->isDone_.load();
}//TP::Detail::Promise::isDone

template<class T, typename H_T, int32_t BS, typename HA>
T TP::Detail::Promise<T, H_T, BS, HA>::getPackDensity() const {
    assert(context_);
    return context_->vol_.load() / context_->tot_vol_;
}//TP::Detail::Promise::getPackDensity

template<class T, typename H_T, int32_t BS, typename HA>
int32_t TP::Detail::Promise<T, H_T, BS, HA>::getBoxCount() const {
    assert(context_);
    return context_->bcc_.load();
}//TP::Detail::Promise::getBoxCount

template<class T, typename H_T, int32_t BS, typename HA>
int32_t TP::Detail::Promise<T, H_T, BS, HA>::getMissedCount() const {
    assert(context_);
    return context_->mcc_.load();
}//TP::Detail::Promise::getMissedCount

template<class T, typename H_T, int32_t BS, typename HA>
double TP::Detail::Promise<T, H_T, BS, HA>::getTime() const {
    assert(context_);
    const auto now = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> ee = now - context_->time_start_;
    return ee.count();
}//TP::Detail::Promise::getTime

//------------------------------

template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
TP::Detail::Promise<T, H_T, BS, HA> TP::solve(const TP::Config<T, CF, H_T, BS, HA>& _conf) {

    using namespace TP;
    static_assert(std::is_floating_point_v<T>, "T needs to be floating point type!");
    static_assert(CostFunction::isValidCF<T, CF>, "Costfunction has not a valid signature!");
    
    Detail::Promise<T, H_T, BS, HA> out;
    out.context_ = std::make_unique<Detail::SolverContext<T, H_T, BS, HA>>(_conf.NumThreads);
    auto c = out.context_.get();

    c->isDone_ = false;
    c->tot_vol_ = _conf.Bounds.x * _conf.Bounds.y * _conf.Height;
    c->mt_ = std::thread(Detail::run_impl<T, CF, H_T, BS, HA>, std::ref(_conf), c);

    return out;    

}//TP::solve

//------------------------------

template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
void TP::Detail::run_impl(
    Config<T, CF, H_T, BS, HA> _conf,
    SolverContext<T, H_T, BS, HA>* _cont
) {

    using namespace MQT2;
    using Tree = MQT2::MedianQuadTree<H_T, BS, HA>;

    using Dist = ::std::uniform_int_distribution<>;
    using DistD = ::std::uniform_real_distribution<T>;
    using Rand = std::mt19937_64;

    const int32_t ee0 = _conf.Bounds.y + 2;
    const int32_t ee1 = _conf.Bounds.x + 2;

    _cont->N_ = get_power_of_2(std::max(ee0, ee1), Tree::BUCKET_SIZE);

    const bool useRandomBox = _conf.BoxType == BoxGenerationType::RANDOM;

    _cont->map_.resize(_cont->N_ * _cont->N_);
    std::fill(_cont->map_.begin(), _cont->map_.end(), 0);

    for (int32_t n0 = 0; n0 < ee0; ++n0) {
        const int32_t i1 = n0 * _cont->N_;
        const int32_t i2 = (ee1 - 1) + n0 * _cont->N_;
        _cont->map_[i1] = _conf.Height;
        _cont->map_[i2] = _conf.Height;
    }

    for (int32_t n1 = 0; n1 < ee1; ++n1) {
        const int32_t i1 = n1;
        const int32_t i2 = n1 + (ee0 - 1) * _cont->N_;
        _cont->map_[i1] = _conf.Height;
        _cont->map_[i2] = _conf.Height;
    }

    _cont->tree_ = std::make_unique<Tree>(_cont->map_, _cont->N_);

    //--------------------------

    const int32_t bc = (_cont->N_ / Tree::BUCKET_SIZE);
    std::vector<bool> mm;
    mm.resize(bc * bc);
    std::fill(mm.begin(), mm.end(), true);

    int32_t et = 0;

    std::random_device rd;
    int64_t Seed = _conf.UseRandomSeed_ ? rd() : _conf.Seed;
    Rand g(Seed);
    DistD dist = DistD(_conf.MinBoxVolume, _conf.MaxBoxVolume);

    const auto rbox = [&](const T _vol, const bool _c) -> glm::vec<3, T> {
        const T s = std::pow(_vol, 1. / 3.);
        if (_c) return glm::vec<3, T>(s) * T(0.5);
        const T ss = 3 * s;
        const T d1 = DistD(0.2 * ss, ss)(g);
        const T d2 = DistD(0.2 * (ss - d1), ss - d1)(g);
        const T d3 = ss - d1 - d2;
        return { d1, d2, d3 };
    };

    std::vector<glm::vec<3, T>> next_set(_conf.SetSize);
    std::deque<glm::vec<3, T>> next_q;

    if (_conf.BoxType == BoxGenerationType::LIST) {
        const BoxList& list = _conf.BoxList;
        std::vector<glm::vec<3, T>> tmp;
        for (const BoxEntry& p : list.List) {
            for (int32_t i = 0; i < p.Count; ++i)
                tmp.push_back(p.Size);
        }
        if (list.ShuffleBoxes)
            std::shuffle(tmp.begin(), tmp.end(), g);
        for (const auto& f : tmp)
            next_q.push_back(f);
    }

    while (true) {

        if (_conf.BoxType == BoxGenerationType::LIST && next_q.empty()) break;
        if (_cont->isDone_) break;

        std::mutex mut;
        std::vector<Detail::Result<T>> res;
        std::atomic<double> minc = std::numeric_limits<double>::infinity();

        next_set.clear();
        
        switch (_conf.BoxType) {
            case BoxGenerationType::RANDOM:
            {
                for (int32_t i = 0; i < _conf.SetSize; ++i)
                    next_set.push_back(rbox(dist(g), _conf.CubeRandomBoxes));
            }
            break;
            case BoxGenerationType::LIST:
            {
                for (int32_t i = 0; i < _conf.SetSize; ++i) {
                    if (next_q.empty()) break;
                    next_set.push_back(next_q.front() * T(0.5));
                    next_q.pop_front();
                }
            }
            break;
        }
        
        for (size_t i = 0; i < next_set.size(); ++i) {

            const glm::vec<3, T>& nextSize = next_set[i];

            const FBox aabb = FBox{-nextSize, nextSize};
            
            //Z_XY
            dispatch_impl(_conf, _cont, res, mut, minc, i,
                int32_t(std::round(aabb.GetExtent().x)),
                int32_t(std::round(aabb.GetExtent().y)),
                nextSize,
                int32_t(std::round(aabb.GetExtent().z)),
                EAxisPerm::Z_XY_0);

            //Z_YX
            dispatch_impl(_conf, _cont, res, mut, minc, i,
                int32_t(std::round(aabb.GetExtent().y)),
                int32_t(std::round(aabb.GetExtent().x)),
                nextSize,
                int32_t(std::round(aabb.GetExtent().z)),
                EAxisPerm::Z_XY_1);

            //Y_XZ
            dispatch_impl(_conf, _cont, res, mut, minc, i,
                int32_t(std::round(aabb.GetExtent().x)),
                int32_t(std::round(aabb.GetExtent().z)),
                nextSize,
                int32_t(std::round(aabb.GetExtent().y)),
                EAxisPerm::Y_XZ_0);

            //Y_ZX
            dispatch_impl(_conf, _cont, res, mut, minc, i,
                int32_t(std::round(aabb.GetExtent().z)),
                int32_t(std::round(aabb.GetExtent().x)),
                nextSize,
                int32_t(std::round(aabb.GetExtent().y)),
                EAxisPerm::Y_XZ_1);

            //X_YZ
            dispatch_impl(_conf, _cont, res, mut, minc, i,
                int32_t(std::round(aabb.GetExtent().y)),
                int32_t(std::round(aabb.GetExtent().z)),
                nextSize,
                int32_t(std::round(aabb.GetExtent().x)),
                EAxisPerm::X_YZ_0);

            //X_ZY
            dispatch_impl(_conf, _cont, res, mut, minc, i,
                int32_t(std::round(aabb.GetExtent().z)),
                int32_t(std::round(aabb.GetExtent().y)),
                nextSize,
                int32_t(std::round(aabb.GetExtent().x)),
                EAxisPerm::X_YZ_1);

            _cont->tg_.wait();

            if (_cont->isDone_) break;

        }

        if (_cont->isDone_) break;

        if (res.empty()) {
            et++;
            _cont->mcc_++;
            if (useRandomBox && _cont->mcc_ > _conf.MaxEmptryTries) break;
            if (useRandomBox && et > _conf.EmptryTries) break;
            if (!useRandomBox && _conf.EnforceMisses && _cont->mcc_ > _conf.MaxEmptryTries) break;
            if (!useRandomBox && _conf.EnforceMisses && et > _conf.EmptryTries) break;
            continue;
        }

        et = 0;

        std::sort(res.begin(), res.end(), [](const Result<T>& _e1, const Result<T>& _e2) {
            return _e1.weight < _e2.weight;
        });

        auto& r = res[0];

        for (size_t i = 0; i < next_set.size(); ++i) {
            if (i == r.set_index) continue;
            next_q.push_front(next_set[i] * T(2.));
        }

        const FBox tar = FBox{
            glm::vec<3, T>(r.n0 - r.ext.x, r.n1 - r.ext.y, r.h),
            glm::vec<3, T>(r.n0 + r.ext.x, r.n1 + r.ext.y, r.h + 2 * r.ext.z)};

        const auto tr = make_transform<T>(
            r.perm,
            tar,
            glm::vec<3, T>(0.),
            r.ext * T(2.)
        ); 

        //std::osyncstream(std::cout) << "---------" << std::endl;
        //std::osyncstream(std::cout) << glm::to_string(tr) << std::endl;
        //std::osyncstream(std::cout) << "[" << tar.min_.x << ", " << tar.min_.y << ", " << r.h << "][" << tar.max_.x << ", " << tar.max_.y << "]" << std::endl;

        {
            std::lock_guard<std::mutex> lock(_cont->m_data);
            _cont->data_.push_back(tr);
            _cont->bcc_ += 1;
            _cont->vol_ += 4. * r.ext.x * r.ext.y * r.ext.z;
            //std::cout << r.ext.x << ", " << r.ext.y << ", " << r.ext.z << std::endl;
        }

        std::fill(mm.begin(), mm.end(), false);
        for (int32_t n0 = int32_t(tar.min_.x) / Tree::BUCKET_SIZE; n0 < int32_t(tar.max_.x) / Tree::BUCKET_SIZE; ++n0) {
            for (int32_t n1 = int32_t(tar.min_.y) / Tree::BUCKET_SIZE; n1 < int32_t(tar.max_.y) / Tree::BUCKET_SIZE; ++n1) {
                const int32_t iid = n0 + n1 * bc;
                mm[iid] = true;
            }
        }

        for (int32_t n0 = int32_t(tar.min_.x); n0 <= int32_t(tar.max_.x); ++n0) {
            for (int32_t n1 = int32_t(tar.min_.y); n1 <= int32_t(tar.max_.y); ++n1) {
                const int32_t i = n1 + n0 * _cont->N_;
                _cont->map_[i] = tar.max_.z;
            }
        }

        _cont->tree_->recompute(mm);

    }

    _cont->isDone_ = true;
    
}//TP::Detail::run_impl

//------------------------------

template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
void TP::Detail::dispatch_impl(
    const Config<T, CF, H_T, BS, HA>& _conf,
    SolverContext<T, H_T, BS, HA>* _cont,
    std::vector<TP::Detail::Result<T>>& _res,
    std::mutex& _m,
    std::atomic<double>& _minc,
    const int32_t _set_index,
    const int32_t _ext0, const int32_t _ext1,
    const glm::vec<3, T>& _ext_org,
    const int32_t _h, TP::Detail::EAxisPerm _perm
) {
        auto tr = [=, &_conf, &_res, &_m, &_minc] () -> void {
        auto ro = overlap_impl<T, CF, H_T, BS, HA>(_conf, _cont, _ext0, _ext1, _h);
        if (!ro.empty()) {
            for (auto& r : ro) {
                r.set_index = _set_index;
                r.ext = glm::vec<3, T>(_ext0, _ext1, _h);
                r.ext_org = _ext_org;
                r.perm = _perm;
                r.isRandomBox = _conf.BoxType == BoxGenerationType::RANDOM;
                r.weight = CF<T>::eval(r);
                if (r.weight < _minc) {
                    _minc = std::min<T>(r.weight, _minc.load());
                    std::lock_guard<std::mutex> l(_m);
                    _res.push_back(r);
                }
            }
        }
    };

    if (_conf.MultiThreading) _cont->tg_.dispatch(std::move(tr));
    else tr();
}//TP::Detail::dispatch_impl

//------------------------------

template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
std::vector<TP::Detail::Result<T>> TP::Detail::overlap_impl(
    const Config<T, CF, H_T, BS, HA>& _conf,
    SolverContext<T, H_T, BS, HA>* _cont,
    const int32_t _ext0, 
    const int32_t _ext1, 
    const int32_t _h
) {

    using namespace MQT2;

    const int32_t ee0 = _conf.Bounds.y + 2;
    const int32_t ee1 = _conf.Bounds.x + 2;

    std::vector<Result<T>> res;
    for (int32_t n0 = _ext0 + 1; n0 < ee0 - _ext0 - 1; ++n0) {
        for (int32_t n1 = _ext1 + 1; n1 < ee1 - _ext1 - 1; ++n1) {

            const int32_t i = n1 + n0 * _cont->N_;
            if (int32_t(_cont->map_[i]) + 2 * _h >= _conf.Height) continue;

            const auto [l1, m1, h1] = _cont->tree_->check_overlap(
                Vec2{ int32_t(n0 - _ext0), int32_t(n1 - _ext1) },
                Vec2{ int32_t(n0 + _ext0), int32_t(n1 + _ext1) },
                _cont->map_[i]);

            const auto [l2, m2, h2] = _cont->tree_->check_border_overlap(
                Vec2{ int32_t(n0 - _ext0) - 1, int32_t(n1 - _ext1) + 1 },
                Vec2{ int32_t(n0 + _ext0) - 1, int32_t(n1 + _ext1) + 1 },
                _cont->map_[i]);

            if (_conf.AllowOverlap) {
                if (h1 != 0) continue;
            } else {
                if (h1 != 0 || l1 != 0) continue;
            }

            Result<T> out;
            out.n0 = n0;
            out.n1 = n1;
            out.l = l1;
            out.b_l = l2;
            out.b_m = m2;
            out.b_h = h2;
            out.h = _cont->map_[i];

            res.push_back(out);
        }
    }
    return res;
}//TP::Detail::overlap_impl

//------------------------------

template<class T>
glm::mat<4, 4, T> TP::Detail::make_transform (
    const EAxisPerm _perm,
    const FBox<T>& _target,
    const glm::vec<3, T>& _pivot_offset,
    const glm::vec<3, T>& _ext
) {

    const glm::vec<3, T>& ctr = _target.GetCenter();

    switch (_perm) {
        case EAxisPerm::Z_XY_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr - _pivot_offset);
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::Z_XY_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ -_pivot_offset.y, _pivot_offset.x, -_pivot_offset.z });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::Z_XY_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ _pivot_offset.x, _pivot_offset.y, -_pivot_offset.z });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::Z_XY_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ _pivot_offset.y, _pivot_offset.x, -_pivot_offset.z });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        //---------------------
        case EAxisPerm::Y_XZ_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ -_pivot_offset.x, -_pivot_offset.z, _pivot_offset.y });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::Y_XZ_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ _pivot_offset.z, -_pivot_offset.x, _pivot_offset.y });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::Y_XZ_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ _pivot_offset.x, _pivot_offset.z, _pivot_offset.y });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::Y_XZ_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ -_pivot_offset.z, _pivot_offset.x, _pivot_offset.y });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        //---------------------
        case EAxisPerm::X_YZ_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ -_pivot_offset.y, -_pivot_offset.z, -_pivot_offset.x });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::X_YZ_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ -_pivot_offset.z, _pivot_offset.y, -_pivot_offset.x });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::X_YZ_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ _pivot_offset.y, _pivot_offset.z, -_pivot_offset.x });
            tr = glm::scale(tr, _ext);
            return tr;
        }
        case EAxisPerm::X_YZ_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::translate(tr, ctr + glm::vec<3, T>{ _pivot_offset.z, -_pivot_offset.y, -_pivot_offset.x });
            tr = glm::scale(tr, _ext);
            return tr;
        }
    }
    assert(false);
    return {};
}//TP::Detail::make_transform
