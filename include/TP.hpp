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
//#include <syncstream>
#include <queue>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include <optional>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/ext/scalar_constants.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>

#include <lodepng.h>

namespace TP {

    using Rand = std::mt19937_64;

    template<class T>
    using Dist = std::uniform_real_distribution<T>;

    template<class T>
    [[nodiscard]] constexpr bool aeq(const T _v1, const T _v2) noexcept {
        if constexpr(std::is_integral_v<T>) return _v1 == _v2;
        else return std::abs( _v1 - _v2 ) <= std::numeric_limits<T>::epsilon() * std::max<T>(T(1.), std::max<T>(_v1, _v2));
    }//aeq

    namespace Detail {
        template<class, typename, typename, uint32_t, typename> 
        class Promise;
    }

    template<class, template<typename> class, class, class, uint32_t, class>
    struct Config;

    //-------------------------

    constexpr uint32_t PF_Z_XY = 0x1;
    constexpr uint32_t PF_Z_YX = 0x2;

    constexpr uint32_t PF_Y_XZ = 0x4;
    constexpr uint32_t PF_Y_ZX = 0x8;

    constexpr uint32_t PF_X_YZ = 0x10;
    constexpr uint32_t PF_X_ZY = 0x20;

    constexpr uint32_t PF_ALL = PF_Z_XY | PF_Z_YX | PF_Y_XZ | PF_Y_ZX | PF_X_YZ | PF_X_ZY;

    //----------------------

    template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
    [[nodiscard]] Detail::Promise<T, H_T, R_T, BS, HA> solve(Config<T, CF, H_T, R_T, BS, HA>& _config);

    namespace Detail {

        template<class R, bool LogScaling>
        void image_real(
            const int32_t _n0,
            const int32_t _n1,
            R* _buffer,
            const std::string& _filename,
            const std::string& _folder
        ) {
            {
                const std::filesystem::path path = std::filesystem::path(_folder);
                if (!std::filesystem::exists(path))
                    std::filesystem::create_directory(path);
            }

            double max = -std::numeric_limits<double>::infinity();
            for (int32_t i = 0; i < _n0 * _n1; ++i)
                max = std::max<double>(max, std::abs(double(_buffer[i])));

            const double frac = 1. / max;
            std::vector<unsigned char> img;
            for (int32_t n1 = 0; n1 < _n1; ++n1) {
                for (int32_t n0 = 0; n0 < _n0; ++n0) {	
                    const int32_t i1 = n0 + n1 * _n0;
                    unsigned char c = 0;

                    if constexpr (LogScaling) {
                        c = std::abs(1. - std::abs(std::log(std::abs(1. + double(_buffer[i1])))) / std::log(max));
                    } else {
                        c = 255 * (std::abs(double(_buffer[i1])) * frac);
                    }

                    img.push_back(c);
                    img.push_back(c);
                    img.push_back(c);
                    img.push_back(255);
                }
            }
            const std::filesystem::path path = std::filesystem::path(_folder) / _filename;
            lodepng::encode(path.string(), img.data(), _n0, _n1);
        }//Layouter::Detail::image_real

        [[nodiscard]] inline uint32_t get_power_of_2(const uint32_t _target, const uint32_t _base) noexcept {
            uint32_t t = _base;
            while (t < _target)
                t  = t << 1;
            return t;
        }//get_power_2

        template<int32_t DIM, class T>
        struct FBox {
            glm::vec<DIM, T> min_;
            glm::vec<DIM, T> max_;
            //---------------------------
            [[nodiscard]] glm::vec<DIM, T> GetExtent() const noexcept { return (max_ - min_) * T(0.5); }
            [[nodiscard]] glm::vec<DIM, T> GetCenter() const noexcept { return min_ + GetExtent(); }
            [[nodiscard]] glm::vec<DIM, T> getSize() const noexcept { return (max_ - min_); }
            [[nodiscard]] T minWidth() const noexcept { return std::min(max_.x - min_.x, max_.y - min_.y); }
            [[nodiscard]] T getArea() const noexcept { return (max_.x - min_.x) * (max_.y - min_.y); }
            [[nodiscard]] T getAspectRatio() const noexcept { return (max_.x - min_.x) / (max_.y - min_.y); } 
        };//FBox

        template<int32_t DIM, class T>
        [[nodiscard]] inline FBox<DIM, T> operator*(const FBox<DIM, T>& _q, const T _frac) noexcept {
            const T f = std::sqrt(_frac);
            const auto c = _q.GetCenter() * f;
            const auto e = _q.GetExtent() * f;
            return FBox<DIM, T>{ c - e, c + e };
        }//operator*

        template<int32_t DIM, class T>
        [[nodiscard]] inline bool operator<(const FBox<DIM, T>& _q1, const FBox<DIM, T>& _q2) noexcept {
            return _q1.getArea() < _q2.getArea();
        }//operator>

        //-------------------------

        template<class T>
        struct Interval {
            T min_, max_;
        };//Interval

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

        template<class T, bool NEGATIVE_WEIGHTS = true>
        class Annealer {

            static_assert(std::is_floating_point_v<T>);

            std::mt19937_64 rand_;
            std::vector<T> weights_;
            std::vector<T> l_weights_;

        public:
            Annealer() = default;
            Annealer(const uint64_t _seed, const size_t _weights_count) : rand_(std::mt19937_64(_seed)) {
                weights_.resize(_weights_count);
                l_weights_.resize(_weights_count);

                std::uniform_real_distribution<T> dist(NEGATIVE_WEIGHTS ? T(-1.) : T(0.), T(1.));
                for(size_t i = 0; i < _weights_count; ++i){
                    weights_[i] = dist(rand_);
                }
            }

            [[nodiscard]] size_t get_size() const noexcept { return weights_.size(); }

            void step(const size_t sub_set_size, const T _max_delta) {

                assert(sub_set_size <= get_size());
                std::memcpy(l_weights_.data(), weights_.data(), sizeof(T) * get_size());

                std::vector<size_t> is;
                is.reserve(get_size());
                for(size_t i = 0; i < get_size(); ++i)
                    is[i] = i;
                std::shuffle(is.begin(), is.end(), rand_);

                std::uniform_real_distribution<T> dist(NEGATIVE_WEIGHTS ? -_max_delta : T(0.), _max_delta);

                for(size_t i = 0; i < sub_set_size; ++i){  
                    weights_[is[i]] = weights_[is[i]] + dist(rand_);
                }
            }

            void punch_through(const T _min_delta, const T _max_delta) {
                std::uniform_real_distribution<T> dist(_min_delta, _max_delta);
                std::uniform_int_distribution<> sign (0, 1);
                for(size_t i = 0; i < get_size(); ++i){  
                    weights_[i] = weights_[i] + (sign(rand_) ? 1 : -1) * dist(rand_);
                }
            }

            void reverse() {
                std::memcpy(weights_.data(), l_weights_.data(), sizeof(T) * get_size());
            }

            [[nodiscard]] const T& operator[](const size_t _idx) const {
                assert(_idx < get_size());
                return weights_[_idx];
            }

            [[nodiscard]] std::vector<T> cpy() const { return weights_; }

            void set(const std::vector<T>& _w) {
                weights_ = l_weights_ = _w;
            }

        };//Annealer

        //-------------------------

        template<class T>
        class Topology {
            std::vector<FBox<3, T>> bxs_;
        public:
            void push(const FBox<3, T>& _bx);
            [[nodiscard]] std::optional<FBox<3, T>> sample(const glm::vec<3, T>& _pos) const;
            [[nodiscard]] std::optional<FBox<3, T>> under(const glm::vec<3, T>& _pos) const;
            [[nodiscard]] std::vector<std::pair<bool, FBox<3, T>>> under(const FBox<3, T>& _box) const;
        };//Topology

        //-------------------------

        enum class EAxisPerm : uint8_t {
            Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
            Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
            X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,
        };//EAxisPerm

        template<class T>
        [[nodiscard]] glm::mat<4, 4, T> make_transform(
            const EAxisPerm _perm,
            const FBox<3, T>& _target,
            const glm::vec<3, T>& _pivot_offset,
            const glm::vec<3, T>& _ext
        );//make_transform

        //-------------------------

        template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
        struct Bin {

            int32_t id_ = -1;
            //-------------------------
            std::vector<H_T> map_;
            std::unique_ptr<MQT2::MedianQuadTree<H_T, R_T, BS, HA>> tree_;
            Topology<T> topo_;
            //-------------------------
            int32_t N_ = 0;
            std::atomic<int32_t> bcc_ = 0;
            std::atomic<double> vol_ = 0.;
            T tot_vol_ = T(0.);
            //-------------------------
            uint64_t seed_ = 0;

            Bin() = default;
            Bin(Bin<T, H_T, R_T, BS, HA>&&) = default;
            Bin<T, H_T, R_T, BS, HA>& operator=(Bin<T, H_T, R_T, BS, HA>&&) = default;

            // kiff alert: dont use these and dont make copies. 
            // Atomics can't be moved or copied. This works as a default constructor
            // and are needed for std::vector resize and reserve to work
            explicit Bin(const Bin<T, H_T, R_T, BS, HA>&) {}   
            Bin<T, H_T, R_T, BS, HA>& operator=(const Bin<T, H_T, R_T, BS, HA>&) { return *this; }
            
        };//Bin

        //-------------------------

        template<class T>
        struct Entry {
            int32_t id_ = -1;
            int32_t bin_id_ = 0;
            glm::mat<4, 4, T> tf_;
        };//Entry

        //-------------------------

        template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
        struct SolverContext {

            std::thread mt_;
            TaskGraph tg_;

            std::vector<Entry<T>> data_;
            std::mutex m_data;

            std::vector<Bin<T, H_T, R_T, BS, HA>> bins_;
            T tot_vol_ = T(0.);

            uint64_t seed_;

            std::condition_variable cv_;

            std::chrono::high_resolution_clock::time_point time_start_;
            std::chrono::high_resolution_clock::time_point time_end_;
            std::atomic<bool> isDone_ = true;
            std::atomic<bool> terminate_ = false;

            std::atomic<int32_t> mcc_ = 0;

            SolverContext(const int32_t _num_threads) : tg_(TaskGraph(_num_threads)) {}
            SolverContext(const SolverContext&) = delete;
            SolverContext(SolverContext&&) = delete;
            ~SolverContext() { mt_.join(); }

        };//SolverContext

        //-------------------------

        template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
        class Promise {

            std::unique_ptr<SolverContext<T, H_T, R_T, BS, HA>> context_;      
            //-------------
            Promise() = default;
            //-------------
            template<typename T_, template<typename> class CF_, typename H_T_, typename R_T_, uint32_t BS_, typename HA_>
            friend Promise<T_, H_T_, R_T_, BS_, HA_> TP::solve(::TP::Config<T_, CF_, H_T_, R_T_, BS_, HA_>& _conf);
        
        public:

            Promise(Promise&&) = default;
            Promise(const Promise&) = delete;
            Promise& operator=(Promise&&) = default;
            Promise& operator=(const Promise&) = delete;

            //-------------

            //locks the thread until isDone() == true
            void wait();

            //stops the packing and shutdowns all threads safely
            void stop();

            //NOT thread safe! only call when isDone() returns true!
            [[nodiscard]] std::vector<Entry<T>>& data() const;

            //NOT thread safe! only call when isDone() returns true!
            [[nodiscard]] std::vector<Entry<T>>& data();

            //thread safe! returns a copy of the data
            [[nodiscard]] std::vector<Entry<T>> data_cpy();

            //thread safe!
            [[nodiscard]] bool isDone() const;

            //thread safe! returns the current pack density [0, 1]
            [[nodiscard]] T getTotalPackDensity() const;

            //thread safe! returns the current pack density [0, 1]
            [[nodiscard]] T getPackDensity(const size_t _bin) const;

            //thread safe!
            [[nodiscard]] uint32_t getTotalBoxCount() const;

            //thread safe!
            [[nodiscard]] uint32_t getBoxCount(const size_t _bin) const;

            //thread safe!
            [[nodiscard]] uint32_t getMissedCount() const;

            //thread safe! seconds since start
            [[nodiscard]] double getTime() const;

            //thread safe!
            [[nodiscard]] uint64_t getSeed() const;

        };//Promise

        //-------------------------

        enum class BoxGenerationType : uint8_t {
            RANDOM, LIST, VALIDATE
        };//PackType

        //-------------------------

        template<class T>
        struct BoxEntry {
            glm::vec<3, T> size_;
            uint32_t count_ = 0;
            uint32_t perms_ = 0x3F; //permutations
        };//BoxEntry

        template<class T>
        struct BoxList {
            std::vector<BoxEntry<T>> list_;
            bool shuffle_boxes_ = false;
        };//BoxList

        //-------------------------

        template<class T>
        struct BinInfo {

            //the bounds of the domain [default: none]
            glm::vec<2, T> Bounds;

            //the max height [default: 150]
            uint32_t Height = 150;

        };//BinInfo

        //------------------------- TODO: document those members!!
        template<class T>
        struct NormalizedResult {
            T bin;
            T n0, n1, h;
            T l, m, b_l, b_m, b_h;
        };//NormalizedResult

        template<class T>
        struct Result {
            uint32_t ext0_, ext1_, height_;
            uint32_t bin_count_;
            int32_t id;
            bool isRandomBox;
            //-----------
            Topology<T>* topo;
            //-----------
            int32_t bin;          
            uint32_t n0, n1, h;
            uint32_t l, b_l, b_m, b_h;
            glm::vec<3, T> ext;
            glm::vec<3, T> ext_org;
            EAxisPerm perm;
            //-----------
            T weight;
            size_t set_index_; //internal use only

            [[nodiscard]] NormalizedResult<T> normalize() const {
                NormalizedResult<T> out;
                out.bin = T(1.) / T(bin_count_) * T(bin);
                out.n0 = T(1.) / T(ext0_) * T(n0);
                out.n1 = T(1.) / T(ext1_) * T(n1);
                out.h = T(1.) / T(height_) * T(h);

                out.h = T(1.) / T(ext.x * ext.y) * T(l);
                out.m = T(1.) - out.h;

                out.b_l = T(1.) / T(2 * (ext.x + ext.y)) * T(b_l);
                out.b_m = T(1.) / T(2 * (ext.x + ext.y)) * T(b_m);
                out.b_h = T(1.) / T(2 * (ext.x + ext.y)) * T(b_h);
                return out;
            }
        };//Result

        //-------------------------

        template<class T>
        void impl_squarify(
            std::vector<std::pair<int32_t, FBox<2, T>>>& _out, 
            const FBox<2, T>& _r,
            std::deque<std::pair<int32_t, T>>& _children, 
            std::vector<std::pair<int32_t, T>>& _row, 
            const int32_t _dir,
            const T _width
        );//impl_squarify

    }//Detail

    template<class T>
    [[nodiscard]] std::vector<std::pair<int32_t, Detail::FBox<2, int32_t>>> 
    squarify(
        const Detail::FBox<2, T>& _bounds, 
        const std::vector<std::pair<int32_t, Detail::Interval<T>>>& _boxes, 
        const uint64_t _seed
    );//squarify

    namespace CostFunction {

        template<class T>
        struct CF_Constant {
            static T eval(const Detail::Result<T>& _r) { 
                return 0.; 
            };
        };//CF_Constant

        template<class T>
        struct CF_BottomLeft {
            static T eval(const Detail::Result<T>& _r) { 
                return std::pow(_r.n0 + _r.n1, 2); 
            };
        };//CF_BottomLeft

        template<class T>
        struct CF_BottomLeftHeight {
            static T eval(const Detail::Result<T>& _r) { 
                return std::pow(double(_r.h), 3) + std::pow(_r.n0 + _r.n1, 2); 
            };
        };//CF_BottomLeftHeight

        template<class T>
        struct CF_Krass {
            static T eval(const Detail::Result<T>& _r) {
                const double c1 = std::pow(_r.n0 + _r.n1, 1);
                const double c2 = std::pow(double(_r.h), 3);
                const double c3 = std::pow(_r.b_l + _r.b_h, 3);
                const double c4 = std::pow(_r.l, 3);
                const double c5 = std::pow(_r.ext.x * _r.ext.y, 2);
                const double c6 = std::pow(double(_r.bin * 10), 3);
                return c1 + c2 - c3 + c4 - c5 + c6;
             };
        };//CF_Krass

        //------------------------------

        template<typename T, template<typename> class CF>
        concept isValidCF = requires(Detail::Result<T> a) { { CF<T>::eval(a) } -> std::same_as<T>; };

    };//CostFunction

    template<class T, template<typename> class COSTFUNCTION, class HEIGHTMAP_T = uint16_t, class R_T = uint32_t, uint32_t BUCKET_SIZE = 15, class HEIGHTMAP_ALLOCATOR = std::allocator<HEIGHTMAP_T>>
    struct Config {
        using T_ = T;
        using HEIGHTMAP_T_ = HEIGHTMAP_T;
        using R_T_ = R_T;
        constexpr static uint32_t BUCKET_SIZE_ = BUCKET_SIZE;
        using HEIGHTMAP_ALLOCATOR_ = HEIGHTMAP_ALLOCATOR;
        using COSTFUNCTION_ = COSTFUNCTION<T>;
        //----------------------

        //if the solver should use worker threads [default: true]
        bool MultiThreading = true;

        //number of worker threads [default: 4]
        uint32_t NumThreads = 4;

        //if the solver should use a random seed [default: true]
        bool UseRandomSeed = true;

        //the seed [default: 1234567890]
        uint64_t Seed = 1234567890;

        std::vector<Detail::BinInfo<T>> Bins;

        //if random boxes or a pre-defined list of boxes should be used [default: Random]
        Detail::BoxGenerationType BoxType = Detail::BoxGenerationType::RANDOM;

        //if boxes are allowed to overlap [default: false]
        bool AllowOverlap = false;

        //how many misses in a row before terminating [default: 25]
        uint32_t EmptryTries = 25;

        //how many misses in total before terminating [default: 50]
        uint32_t MaxEmptryTries = 50;

        //How far the solver can look ahead [default: 1]
        uint32_t LookAheadSize = 1;

        //--------------------------

        //
        uint32_t AllowedPermutations = PF_ALL;

        //if the random boxes should be cubes [default: true]
        bool CubeRandomBoxes = true;

        //the minimum volume of the random boxes [default: 250]
        T MinBoxVolume = 250.;

        //the maximum volume of the random boxes [default: 1500]
        T MaxBoxVolume = 1500.;

        //--------------------------

        //the list fo predefined boxes
        Detail::BoxList<T> BoxList;

        //if it should enforce emptry tries. if false the solver runs until the list is empty [default: false]
        bool EnforceMisses = false;

        //--------------------------

        //How many boxes the ground truth should contain [default: 12]
        uint32_t EvalBoxCount = 12;

    };//Config

    namespace Detail {

        template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
        void run_impl(
            const Config<T, CF, H_T, R_T, BS, HA> _conf,
            SolverContext<T, H_T, R_T, BS, HA>* _cont
        );//run_impl

        template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
        [[nodiscard]] std::vector<Result<T>> overlap_impl(
            const Config<T, CF, H_T, R_T, BS, HA>& _conf,
            SolverContext<T, H_T, R_T, BS, HA>* _cont,
            const size_t _bin_index,
            const uint32_t _ext0, 
            const uint32_t _ext1, 
            const uint32_t _h
        );//overlap_impl

        template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
        void dispatch_impl(
            const Config<T, CF, H_T, R_T, BS, HA>& _conf,
            SolverContext<T, H_T, R_T, BS, HA>* _cont,
            std::vector<Result<T>>& _res,
            std::mutex& _m,
            std::atomic<double>& _minc,
            const size_t _set_index,
            const size_t _bin_index,
            const uint32_t _n0, const uint32_t _n1,
            const glm::vec<3, T>& _ext_org, const int32_t _id,
            const uint32_t _h, EAxisPerm _perm
        );//dispatch_impl       

    }//Detail

}//TP

//------------------------------

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
void TP::Detail::Promise<T, H_T, R_T, BS, HA>::wait() {
    std::unique_lock<std::mutex> lock (context_->m_data);
    context_->cv_.wait(lock, [&]{ return context_->isDone_.load(); });
}//TP::Detail::Promise::wait

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
void TP::Detail::Promise<T, H_T, R_T, BS, HA>::stop() {
    assert(context_);
    context_->terminate_ = true;
}//TP::Detail::Promise::wait

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
std::vector<TP::Detail::Entry<T>>& TP::Detail::Promise<T, H_T, R_T, BS, HA>::data() const {
    assert(context_);
    return context_->data_;
}//TP::Detail::Promise::data

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
std::vector<TP::Detail::Entry<T>>& TP::Detail::Promise<T, H_T, R_T, BS, HA>::data() {
    assert(context_);
    return context_->data_;
}//TP::Detail::Promise::data

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
std::vector<TP::Detail::Entry<T>> TP::Detail::Promise<T, H_T, R_T, BS, HA>::data_cpy() {
    assert(context_);
    std::lock_guard<std::mutex> lock(context_->m_data);
    std::vector<TP::Detail::Entry<T>> out = context_->data_;
    return out;
}//TP::Detail::Promise::data

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
bool TP::Detail::Promise<T, H_T, R_T, BS, HA>::isDone() const {
    assert(context_);
    return context_->isDone_.load();
}//TP::Detail::Promise::isDone

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
T TP::Detail::Promise<T, H_T, R_T, BS, HA>::getTotalPackDensity() const {
    assert(context_);
    float tv = 0.f;
    for(const auto& b : context_->bins_)
        tv += b.vol_.load();
    return tv * context_->tot_vol_;
}//TP::Detail::Promise::getPackDensity

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
T TP::Detail::Promise<T, H_T, R_T, BS, HA>::getPackDensity(const size_t _bin) const {
    assert(context_);
    return context_->bins_[_bin].tot_vol_ * context_->bins_[_bin].vol_.load();
}//TP::Detail::Promise::getPackDensity

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
uint32_t TP::Detail::Promise<T, H_T, R_T, BS, HA>::getTotalBoxCount() const {
    assert(context_);
    uint32_t c = 0;
    for(const auto& b : context_->bins_)
        c += b.bcc_.load();
    return c;
}//TP::Detail::Promise::getBoxCount

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
uint32_t TP::Detail::Promise<T, H_T, R_T, BS, HA>::getBoxCount(const size_t _bin) const {
    assert(context_);
    return context_->bins_[_bin].bcc_.load();
}//TP::Detail::Promise::getBoxCount

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
uint32_t TP::Detail::Promise<T, H_T, R_T, BS, HA>::getMissedCount() const {
    assert(context_);
    return context_->mcc_.load();
}//TP::Detail::Promise::getMissedCount

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
double TP::Detail::Promise<T, H_T, R_T, BS, HA>::getTime() const {
    assert(context_);
    const auto now = context_->isDone_.load() ? context_->time_end_ : std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> ee = now - context_->time_start_;
    return ee.count();
}//TP::Detail::Promise::getTime

template<class T, typename H_T, typename R_T, uint32_t BS, typename HA>
uint64_t TP::Detail::Promise<T, H_T, R_T, BS, HA>::getSeed() const {
    assert(context_);
    return context_->seed_;
}//TP::Detail::Promise::getTime

//------------------------------

template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
TP::Detail::Promise<T, H_T, R_T, BS, HA> TP::solve(TP::Config<T, CF, H_T, R_T, BS, HA>& _conf) {

    using namespace TP;
    static_assert(std::is_floating_point_v<T>, "T needs to be floating point type!");
    static_assert(CostFunction::isValidCF<T, CF>, "Costfunction has not a valid signature!");
    
    Detail::Promise<T, H_T, R_T, BS, HA> out;
    out.context_ = std::make_unique<Detail::SolverContext<T, H_T, R_T, BS, HA>>(_conf.NumThreads);
    auto c = out.context_.get();

    std::random_device rd;
    _conf.Seed = c->seed_ = _conf.UseRandomSeed ? rd() : _conf.Seed;
    c->isDone_ = false;
    c->terminate_ = false;

    Rand rng (_conf.Seed);

    assert(!_conf.Bins.empty());
    c->bins_.resize(_conf.Bins.size());
    for(size_t i = 0; i < _conf.Bins.size(); ++i){
        Detail::Bin<T, H_T, R_T, BS, HA> b;
        b.id_ = int32_t(i);
        b.seed_ = rng();
        c->bins_[i] = std::move(b);
    }

    c->mt_ = std::thread(Detail::run_impl<T, CF, H_T, R_T, BS, HA>, _conf, c);

    return out;    

}//TP::solve

//------------------------------

template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
void TP::Detail::run_impl(
    const Config<T, CF, H_T, R_T, BS, HA> _conf,
    SolverContext<T, H_T, R_T, BS, HA>* _cont
) {

    using namespace MQT2;
    using Tree = MQT2::MedianQuadTree<H_T, R_T, BS, HA>;

    using Dist = ::std::uniform_int_distribution<>;
    using DistD = ::std::uniform_real_distribution<T>;
    using Rand = std::mt19937_64;

    _cont->time_start_ = std::chrono::high_resolution_clock::now();

    const bool useRandomBox = _conf.BoxType == BoxGenerationType::RANDOM;

    std::vector<glm::vec<2, uint32_t>> ee;
    ee.resize(_conf.Bins.size());

    for(size_t i = 0; i < _conf.Bins.size(); ++i){

        const auto& bc = _conf.Bins[i];
        auto& bi = _cont->bins_[i];

        const int32_t ee0 = bc.Bounds.y + 2;
        const int32_t ee1 = bc.Bounds.x + 2;

        ee[i] = { ee0, ee1 };

        _cont->tot_vol_ += T(bc.Height * bc.Bounds.x * bc.Bounds.y);
        bi.tot_vol_ = T(1.) / T(bc.Height * bc.Bounds.x * bc.Bounds.y);
        bi.N_ = Detail::get_power_of_2(std::max(ee0, ee1), Tree::BUCKET_SIZE);

        bi.map_.resize(bi.N_ * bi.N_);
        std::fill(bi.map_.begin(), bi.map_.end(), 0);

        for (size_t n0 = 0; n0 < ee0; ++n0) {
            const size_t i1 = n0 * bi.N_;
            const size_t i2 = (ee1 - 1) + n0 * bi.N_;
            bi.map_[i1] = bc.Height;
            bi.map_[i2] = bc.Height;
        }

        for (size_t n1 = 0; n1 < ee1; ++n1) {
            const size_t i1 = n1;
            const size_t i2 = n1 + (ee0 - 1) * bi.N_;
            bi.map_[i1] = bc.Height;
            bi.map_[i2] = bc.Height;
        }

        bi.tree_ = std::make_unique<Tree>(bi.map_, bi.N_);

    }

    _cont->tot_vol_ = T(1.) / _cont->tot_vol_;

    //--------------------------

    std::vector<std::vector<bool>> mm;
    mm.resize(_conf.Bins.size());

    for(size_t i = 0; i < _conf.Bins.size(); ++i){
        const auto& bc = _conf.Bins[i];
        const auto& bi = _cont->bins_[i];

        const uint32_t bbc = (bi.N_ / Tree::BUCKET_SIZE);
        
        mm[i].resize(bbc * bbc);
        std::fill(mm[i].begin(), mm[i].end(), true);
    }

    uint32_t et = 0;

    std::random_device rd;
    const int64_t Seed = _cont->seed_;
    Rand g(Seed);
    DistD dist = DistD(_conf.MinBoxVolume, _conf.MaxBoxVolume);

    const auto rbox = [&](const T _vol, const bool _c) -> glm::vec<3, T> {
        const T s = std::pow(_vol, 1. / 3.);
        if (_c) return glm::vec<3, T>(s) * T(0.5);
        const T ss = 3 * s;
        const T d1 = DistD(0.2 * ss, ss)(g);
        const T d2 = DistD(0.2 * (ss - d1), ss - d1)(g);
        const T d3 = ss - d1 - d2;
        return { std::max(T(2.), d1), std::max(T(2.), d2), std::max(T(2.), d3) };
    };

    const uint32_t dp = _conf.AllowedPermutations;

    std::vector<std::tuple<glm::vec<3, T>, uint32_t, int32_t>> next_set(_conf.LookAheadSize);
    std::deque<std::tuple<glm::vec<3, T>, uint32_t, int32_t>> next_q; //<bounds, perms, id>

    if (_conf.BoxType == BoxGenerationType::LIST) {
        const BoxList<T>& list = _conf.BoxList;
        std::vector<std::tuple<glm::vec<3, T>, uint32_t, int32_t>> tmp;
        for(size_t l = 0; l < list.list_.size(); ++l){
            const BoxEntry<T>& p = list.list_[l];
            for (uint32_t i = 0; i < p.count_; ++i)
                tmp.push_back({ p.size_, p.perms_, int32_t(l) });
        }
        if (list.shuffle_boxes_)
            std::shuffle(tmp.begin(), tmp.end(), g);
        for (const auto& f : tmp)
            next_q.push_back(f);
    }

    while (true) {

        if (_conf.BoxType == BoxGenerationType::LIST && next_q.empty()) break;
        if (_cont->terminate_) break;

        std::mutex mut;
        std::vector<Detail::Result<T>> res;
        std::atomic<double> minc = std::numeric_limits<double>::infinity();

        next_set.clear();
        
        switch (_conf.BoxType) {
            case BoxGenerationType::RANDOM:
            {
                for (uint32_t i = 0; i < _conf.LookAheadSize; ++i)
                    next_set.push_back({ rbox(dist(g), _conf.CubeRandomBoxes), dp, -1 });
            }
            break;
            case BoxGenerationType::LIST:
            {
                for (uint32_t i = 0; i < _conf.LookAheadSize; ++i) {
                    if (next_q.empty()) break;
                    const auto& [bb, pp, ii] = next_q.front();
                    next_set.push_back({ bb * T(0.5), pp, ii });
                    next_q.pop_front();
                }
            }
            break;
        }
        
        for(size_t j = 0; j < _conf.Bins.size(); ++j){

            for (size_t i = 0; i < next_set.size(); ++i) {

                const auto& [nextSize, nextPerm, nextId] = next_set[i];

                const FBox<3, T> aabb = FBox<3, T>{-nextSize, nextSize};

                const uint32_t n0 = uint32_t(std::round(aabb.GetExtent().x));
                const uint32_t n1 = uint32_t(std::round(aabb.GetExtent().y));
                const uint32_t n2 = uint32_t(std::round(aabb.GetExtent().z));
                
                //Z_XY
                const int32_t Z_XY = (nextPerm & PF_Z_XY); //todo: check on gcc if strange behaviour can be reproduced
                if(Z_XY){
                    dispatch_impl(_conf, _cont, res, mut, minc, i, j,
                        n0, n1, nextSize, nextId, n2, EAxisPerm::Z_XY_0);
                }

                //Z_YX
                const int32_t Z_YX = (nextPerm & PF_Z_YX);
                if(Z_YX){
                    dispatch_impl(_conf, _cont, res, mut, minc, i, j,
                        n1, n0, nextSize, nextId, n2, EAxisPerm::Z_XY_1);
                }

                //Y_XZ
                const int32_t Y_XZ = (nextPerm & PF_Y_XZ);
                if(Y_XZ) {
                    dispatch_impl(_conf, _cont, res, mut, minc, i, j,
                        n0, n2, nextSize, nextId, n1, EAxisPerm::Y_XZ_0);
                }

                //Y_ZX
                const int32_t Y_ZX = (nextPerm & PF_Y_ZX);
                if(Y_ZX){
                    dispatch_impl(_conf, _cont, res, mut, minc, i, j,
                        n2, n0, nextSize, nextId, n1,  EAxisPerm::Y_XZ_1);
                }

                //X_YZ
                const int32_t X_YZ = (nextPerm & PF_X_YZ);
                if(X_YZ){
                    dispatch_impl(_conf, _cont, res, mut, minc, i, j,
                        n1, n2, nextSize, nextId, n0, EAxisPerm::X_YZ_0);
                }

                //X_ZY
                const int32_t X_ZY = (nextPerm & PF_X_ZY);
                if(X_ZY){
                    dispatch_impl(_conf, _cont, res, mut, minc, i, j,
                        n2, n1, nextSize, nextId, n0, EAxisPerm::X_YZ_1);
                }

                if (_cont->terminate_) break;

            }
            if (_cont->terminate_) break;
        }
 
        _cont->tg_.wait();

        if (_cont->terminate_) break;

        //TODO misses when list mode
        if(res.empty()){
            if (((!useRandomBox && _conf.EnforceMisses) || useRandomBox)) {
                et++;
                _cont->mcc_ += next_set.size();
                if ( _cont->mcc_ > _conf.MaxEmptryTries) break;
                if (et > _conf.EmptryTries) break;   
            }
            continue;
        }

        et = 0;

        std::sort(res.begin(), res.end(), [](const Result<T>& _e1, const Result<T>& _e2) {
            return _e1.weight < _e2.weight;
        });

        const auto& r = res[0];

        for (size_t i = 0; i < next_set.size(); ++i) {
            if (i == r.set_index_) continue;
            const auto& [ns, np, ni] = next_set[i];
            next_q.push_front({ ns * T(2.), np, ni} );
        }

        const FBox<3, T> tar = FBox<3, T>{
            glm::vec<3, T>(r.n0 - r.ext.x, r.n1 - r.ext.y, r.h),
            glm::vec<3, T>(r.n0 + r.ext.x, r.n1 + r.ext.y, r.h + 2 * r.ext.z)};

        Entry<T> e;
        e.id_ = r.id;
        e.bin_id_ = r.bin;
        e.tf_ = make_transform<T>(r.perm, tar, glm::vec<3, T>(0.), r.ext * T(2.)); 

        const auto& bi = _conf.Bins[e.bin_id_];
        auto& bn = _cont->bins_[e.bin_id_];

        {
            auto& bi = _cont->bins_[e.bin_id_];
            bi.topo_.push(tar);
            std::lock_guard<std::mutex> lock(_cont->m_data);
            _cont->data_.push_back(e);
            bi.bcc_ += 1;
            const auto si = r.ext * T(2.);
            bi.vol_ = bi.vol_ +  si.x * si.y * si.z;
        }

        for (uint32_t n0 = uint32_t(tar.min_.x); n0 < uint32_t(tar.max_.x); ++n0) {
            for (uint32_t n1 = uint32_t(tar.min_.y); n1 < uint32_t(tar.max_.y); ++n1) {
                const size_t i = n1 + n0 * bn.N_;
                bn.map_[i] = tar.max_.z;
            }
        }

        {  
            auto& mmm = mm[e.bin_id_];
            const uint32_t bbc = (bn.N_ / Tree::BUCKET_SIZE);
            std::fill(mmm.begin(), mmm.end(), false);
            for (uint32_t n0 = uint32_t(tar.min_.x) / Tree::BUCKET_SIZE; n0 <= uint32_t(tar.max_.x) / Tree::BUCKET_SIZE; ++n0) {
                for (uint32_t n1 = uint32_t(tar.min_.y) / Tree::BUCKET_SIZE; n1 <= uint32_t(tar.max_.y) / Tree::BUCKET_SIZE; ++n1) {
                    const size_t iid = n0 + n1 * bbc;
                    mmm[iid] = true;
                }
            }
            bn.tree_->recompute(mmm);
        }

    }

    _cont->time_end_ = std::chrono::high_resolution_clock::now();
    _cont->terminate_ = true;
    _cont->isDone_ = true;
    _cont->cv_.notify_all();
    
}//TP::Detail::run_impl

//------------------------------

template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
void TP::Detail::dispatch_impl(
    const Config<T, CF, H_T, R_T, BS, HA>& _conf,
    SolverContext<T, H_T, R_T, BS, HA>* _cont,
    std::vector<TP::Detail::Result<T>>& _res,
    std::mutex& _m,
    std::atomic<double>& _minc,
    const size_t _set_index,
    const size_t _bin_index,
    const uint32_t _ext0, const uint32_t _ext1,
    const glm::vec<3, T>& _ext_org, const int32_t _id,
    const uint32_t _h, TP::Detail::EAxisPerm _perm
) {
        auto tr = [=, &_conf, &_res, &_m, &_minc] () -> void {
        auto ro = overlap_impl<T, CF, H_T, R_T, BS, HA>(_conf, _cont, _bin_index,  _ext0, _ext1, _h);
        if (!ro.empty()) {
            for (auto& r : ro) {
                r.id = _id;
                r.set_index_ = _set_index;
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

template<typename T, template<typename> class CF, typename H_T, typename R_T, uint32_t BS, typename HA>
std::vector<TP::Detail::Result<T>> TP::Detail::overlap_impl(
    const Config<T, CF, H_T, R_T, BS, HA>& _conf,
    SolverContext<T, H_T, R_T, BS, HA>* _cont,
    const size_t _bin_index,
    const uint32_t _ext0, 
    const uint32_t _ext1, 
    const uint32_t _h
) {

    using namespace MQT2;
    using Vec2i = Vec2<R_T>;

    const auto& bc = _conf.Bins[_bin_index];
    const auto& bi = _cont->bins_[_bin_index];

    const uint32_t ee0 = bc.Bounds.y + 2;
    const uint32_t ee1 = bc.Bounds.x + 2;

    if(_ext0 >= ee0 || _ext1 >= ee1) return {};

    std::vector<Result<T>> res;
    for (int32_t n0 = _ext0 + 1; n0 < ee0 - _ext0 - 1; ++n0) {
        for (int32_t n1 = _ext1 + 1; n1 < ee1 - _ext1 - 1; ++n1) {

            const size_t i = n1 + n0 * bi.N_;
            assert(i < bi.map_.size());
            if (int32_t(bi.map_[i]) + 2 * _h >= bc.Height) continue;

            const auto [l1, m1, h1] = bi.tree_->check_overlap(
                Vec2i{ R_T(n0 - _ext0), R_T(n1 - _ext1) },
                Vec2i{ R_T(n0 + _ext0), R_T(n1 + _ext1) },
                bi.map_[i]);

            const auto [l2, m2, h2] = bi.tree_->check_border_overlap(
                Vec2i{ R_T(n0 - _ext0) - 1, R_T(n1 - _ext1) + 1 },
                Vec2i{ R_T(n0 + _ext0) - 1, R_T(n1 + _ext1) + 1 },
                bi.map_[i]);

            if (_conf.AllowOverlap) {
                if (h1 != 0) continue;
            } else {
                if (h1 != 0 || l1 != 0) continue;
            }

            Result<T> out;
            out.bin_count_ = _conf.Bins.size();
            out.ext0_ = bc.Bounds.y;
            out.ext1_ = bc.Bounds.x;
            out.height_ = bc.Height;
            out.topo = const_cast<Topology<T>*>(&bi.topo_);
            out.bin = _bin_index;
            out.n0 = n0;
            out.n1 = n1;
            out.l = l1;
            out.b_l = l2;
            out.b_m = m2;
            out.b_h = h2;
            out.h = bi.map_[i];

            res.push_back(out);
        }
    }
    return res;
}//TP::Detail::overlap_impl

//------------------------------

template<class T>
glm::mat<4, 4, T> TP::Detail::make_transform (
    const EAxisPerm _perm,
    const FBox<3, T>& _target,
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

//-----------------------------------------

template<class T>
std::vector<std::pair<int32_t, TP::Detail::FBox<2, int32_t>>> TP::squarify(
    const Detail::FBox<2, T>& _bounds, 
    const std::vector<std::pair<int32_t, Detail::Interval<T>>>& _boxes,
    const uint64_t _seed
) {

    using namespace Detail;

    Rand rng (_seed);
    std::vector<std::pair<int32_t, T>> bxs;
    T vl = 0.;
    for(const auto& [ii, t] : _boxes){
        const T v = Dist<T>(t.min_, t.max_)(rng);
        vl += v;
        bxs.push_back({ ii, v });
    }

    const T ifrac = vl / _bounds.getArea();
    const T frac = 1. / ifrac;

    std::sort(bxs.begin(), bxs.end(), [](const auto& _v1, const auto& _v2){ return _v1.second > _v2.second; });
    std::deque<std::pair<int32_t, T>> c;
    for(const auto t : bxs)
        c.push_back(t);

    const auto sbounds = _bounds * ifrac;
    const auto ext = sbounds.GetExtent();
    const int32_t nd = int32_t(ext.x > ext.y);

    std::vector<std::pair<int32_t, T>> row;
    std::vector<std::pair<int32_t, FBox<2, T>>> res;
    
    Detail::impl_squarify(res, sbounds, c, row, nd, sbounds.minWidth());

    const auto iquad = [](const FBox<2, T>& _b){
        FBox<2, int32_t> out;
        out.min_ = { (int32_t)std::round(_b.min_.x), (int32_t)std::round(_b.min_.y) };
        out.max_ = { (int32_t)std::round(_b.max_.x), (int32_t)std::round(_b.max_.y) };
        return out;
    };

    std::vector<std::pair<int32_t, FBox<2, int32_t>>> out;
    out.reserve(res.size());

    for(auto& r : res)
        out.push_back( { r.first, iquad(r.second * frac) });

    return out;

}//Layouter::Detail::SquarifiedTreemap::squarify

template<class T>
void TP::Detail::impl_squarify(
    std::vector<std::pair<int32_t, FBox<2, T>>>& _out, 
    const FBox<2, T>& _r,
    std::deque<std::pair<int32_t, T>>& _c, 
    std::vector<std::pair<int32_t, T>>& _row, 
    const int32_t _dir,
    const T _w
) {

    const auto worst = [](const std::vector<std::pair<int32_t, T>>& _r, const T _w) {
        T min = std::numeric_limits<T>::infinity();
        T max = -std::numeric_limits<T>::infinity();
        T s = T(0);
        for(const auto&[i, r] : _r){
            min = std::min(min, r);
            max = std::max(max, r);
            s += r;
        }
        if (s == 0.) return std::numeric_limits<T>::infinity();
        const T ww = std::pow(_w, 2);
        const T ss = std::pow(s, 2);
        return std::max( (ww * max) / ss, ss / (ww * min));
    };

    const auto worst_ext = [](const std::vector<std::pair<int32_t, T>>& _r, const T _c, const T _w) {
        T min = std::min(std::numeric_limits<T>::infinity(), _c);
        T max = std::max(-std::numeric_limits<T>::infinity(), _c);
        T s = _c;
        for(const auto&[i, r] : _r){
            min = std::min(min, r);
            max = std::max(max, r);
            s += r;
        }
        if (s == 0.) return std::numeric_limits<T>::infinity();
        const T ww = std::pow(_w, 2);
        const T ss = std::pow(s, 2);
        return std::max( (ww * max) / ss, ss / (ww * min));
    };

    //------------------------------------

    const bool last = _c.empty();

    const auto[idx, c] = last ? std::make_pair(int32_t(0), T(0.)) : _c.front();
    if(!last) _c.pop_front();

    const T r1 = worst(_row, _w);
    const T r2 = worst_ext(_row, c, _w);

    if(!last && r1 > r2){
        _row.push_back({ idx, c });
        impl_squarify(_out, _r, _c, _row, _dir, _w);
    } else {

        if(!last) _c.push_front({ idx, c });

        const T vl = std::accumulate(_row.begin(), _row.end(), T(0.), [](const T _v1, const auto& _v2){
            return _v1 + _v2.second;
        });
        const T ff = vl / _w;
        const T frac = 1. / ff;

        glm::vec<2, T> mi = _r.min_;

        //layout
        switch(_dir){
            case 0://x
            {
                for(const auto& [ii, t] : _row){
                    const T s = t * frac;
                    const glm::vec<2, T> ma = mi + glm::vec<2, T>{s, ff}; 
                    _out.push_back( { ii, FBox<2, T>{ mi, ma } });
                    mi = glm::vec<2, T>{ ma.x, mi.y };
                }
            }
            break;
            case 1://y
            {
                for(const auto& [ii, t] : _row){
                    const T s = t * frac;
                    const glm::vec<2, T> ma = mi + glm::vec<2, T>{ff, s}; 
                    _out.push_back({ ii, FBox<2, T>{ mi, ma } });
                    mi = glm::vec<2, T>{ mi.x, ma.y };
                }
            }
            break;
        }
        
        //next row
        _row.clear();

        FBox<2, T> nr;
        nr.max_ = _r.max_;
        
        switch(_dir){
            case 0:
                nr.min_ = { _r.min_.x, _r.min_.y + ff };
            break;
            case 1:
                nr.min_ = { _r.min_.x + ff, _r.min_.y };
            break;
        }

        const auto ext = nr.GetExtent();
        const int32_t nd = int32_t(ext.x > ext.y);

        if(!last) impl_squarify(_out, nr, _c, _row, nd, nr.minWidth());
    }

}//Layouter::Detail::SquarifiedTreemap::impl_squarify

template<class T>
void TP::Detail::Topology<T>::push(const TP::Detail::FBox<3, T>& _bx) {
    bxs_.push_back(_bx);
}//TP::Detail::Topology::push

template<class T>
std::optional<TP::Detail::FBox<3, T>> TP::Detail::Topology<T>::sample(const glm::vec<3, T>& _pos) const {
    return {};
}//TP::Detail::Topology::sample

template<class T>
std::optional<TP::Detail::FBox<3, T>> TP::Detail::Topology<T>::under(const glm::vec<3, T>& _pos) const {
    return {};
}//TP::Detail::Topology::under

template<class T>
std::vector<std::pair<bool, TP::Detail::FBox<3, T>>> TP::Detail::Topology<T>::under(const TP::Detail::FBox<3, T>& _box) const {
    return {};
}//TP::Detail::Topology::under
