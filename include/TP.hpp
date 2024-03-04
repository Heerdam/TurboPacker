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

#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>
#include <glm/ext/scalar_constants.hpp>

namespace TP {

    template<class, typename, int32_t, typename> 
    class Promise;

    template<class, template<typename> class, class, int32_t, class>
    struct Config;

    //----------------------

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
            [[nodiscard]] glm::vec<3, T> GetExtend() const { return (max_ - min_) * T(0.5); }
            [[nodiscard]] glm::vec<3, T> GetCenter() const { return min_ + GetExtend(); }
        };//FBox

        //-------------------------

        enum class EAxisPerm : uint8_t {
            Z_XY_0, Z_XY_1, Z_XY_2, Z_XY_3,
            Y_XZ_0, Y_XZ_1, Y_XZ_2, Y_XZ_3,
            X_YZ_0, X_YZ_1, X_YZ_2, X_YZ_3,

            Z_n_XY_0, Z_n_XY_1, Z_n_XY_2, Z_n_XY_3,
            Y_n_XZ_0, Y_n_XZ_1, Y_n_XZ_2, Y_n_XZ_3,
            X_n_YZ_0, X_n_YZ_1, X_n_YZ_2, X_n_YZ_3
        };//EAxisPerm

        template<class T>
        [[nodiscard]] glm::mat<4, 4, T> make_transform(
            const EAxisPerm _perm,
            const FBox<T>& _target,
            const glm::vec<3, T>& _pivot_offset
        );

        //-------------------------

        template<class T, typename H_T, int32_t BS, typename HA>
        struct SolverContext {

            SolverContext() = default;
            SolverContext(const SolverContext&) = delete;
            SolverContext(SolverContext&&) = delete;

            std::vector<H_T> map_;
            MQT2::MedianQuadTree<H_T, BS, HA> tree_;

            std::thread mt_;
            std::vector<std::thread> workers_;

            std::vector<glm::mat<4, 4, T>> data_;
            std::mutex m_;

            bool isPacking_ = false;
            double last_time_ = 0.;
            double vol_ = 0.;
            int32_t bcc_ = 0;
            int32_t mcc_ = 0;

        };//SolverContext

        //-------------------------

        template<class T, typename H_T, int32_t BS, typename HA>
        class Promise {
            bool done_;
            std::unique_ptr<std::mutex> m_;
            std::unique_ptr<std::condition_variable> cv_;
            std::unique_ptr<SolverContext<T, H_T, BS, HA>> context_;
            
            //-------------
            Promise();
            Promise(Promise&&) = default;
            Promise(const Promise&) = delete;
            Promise& operator=(Promise&&) = default;
            Promise& operator=(const Promise&) = delete;
            //-------------
            template<typename, template<typename> class CF, typename, int32_t, typename>
            friend Promise<T, H_T, BS, HA> ::TG::solve(const Config<T, CF, H_T, BS, HA>& _conf)
        public:

            //locks the thread until isDone() == true
            void wait();

            //NOT thread safe! only call when isDone() returns true!
            [[nodiscard]] const std::vector<glm::mat<4, 4, T>>& data() const { return parent_->result_; }
            
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

        enum class PackMethod : uint8_t {
            ONLINE
        };//PackType

        enum class BoxGenerationType : uint8_t {
            RANDOM, LIST
        };//PackType

        enum class CostFunction : uint8_t {
            SIMPLE, SIMPLE_HEIGHT, CUSTOM
        };//CostFunction

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
            glm::mat<4, 4, T> trans;
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

    template<class T, template<typename> class COSTFUNCTION, class HEIGHTMAP_T = int16_t, int32_t BUCKET_SIZE = 15, class HEIGHTMAP_ALLOCATOR = std::allocator<T>>
    struct Config {
        using T_ = T;
        using HEIGHTMAP_T_ = HEIGHTMAP_T;
        constexpr static int32_t BUCKET_SIZE_ = BUCKET_SIZE;
        using HEIGHTMAP_ALLOCATOR_ = HEIGHTMAP_ALLOCATOR;
        using COSTFUNCTION_<T> = COSTFUNCTION<T>;
        //----------------------
        bool MultiThreading = true;
        int32_t NumThreads = 4;
        uint64_t Seed = 1234567890;
        glm::vec2 Bounds;
        int32_t Height = 150;
        Detail::PackMethod Method = Detail::PackMethod::ONLINE;
        Detail::BoxGenerationType BoxType = Detail::BoxGenerationType::RANDOM;
        Detail::CostFunction CostFunction = Detail::CostFunction::SIMPLE;
        bool AllowOverlap = true;
        //--------------------------
        int32_t EmptryTries = 25;
        int32_t MaxEmptryTries = 50;
        int32_t SetSize = 1;
        //--------------------------
        bool CubeRandomBoxes = true;
        double MinBoxVolume = 250.;
        double MaxBoxVolume = 1500.;
        //--------------------------
        Detail::BoxList BoxList;
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
            const int32_t _ext0, 
            const int32_t _ext1, 
            const int32_t _h
        );//overlap_impl

        template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
        void dispatch_impl(
            std::vector<Result<T>>& _res,
            std::mutex& _m,
            std::atomic<double>& _minc,
            std::atomic<int32_t>& _c,
            const int32_t _set_index,
            const int32_t _n0, const int32_t _n1,
            const glm::vec<3, T>& _ext_org,
            const int32_t _h, EAxisPerm _perm
        );//dispatch_impl       

    }//Detail

    template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
    [[nodiscard]] Detail::Promise<T, H_T, BS, HA> solve(const Config<T, CF, H_T, BS, HA>& _config);

}//TP

//------------------------------

template<typename T, template<typename> class CF, typename H_T, int32_t BS, typename HA>
TP::Detail::Promise<T, H_T, BS, HA> TP::solve(const TP::Config<T, CF, H_T, BS, HA>& _conf) {

    using namespace TP;
    static_assert(CostFunction::isValidCF<T, CF>, "Costfunction has not a valid signature!");

    Detail::Promise<T, H_T, BS, HA> out;
    out.context_ = std::make_unique<SolverContext<T, H_T, BS, HA>>();
    auto c = out.context_.get();

    c->mt_ = std::thread(Detail::run_impl, _conf, c);

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

    const int32_t ee0 = conf->Bounds.Y + 2;
    const int32_t ee1 = conf->Bounds.X + 2;

    const int32_t N_ = Util::get_power_of_2(std::max(ee0, ee1), Tree::BUCKET_SIZE);

    const bool useRandomBox = conf.BoxType == BoxGenerationType::RANDOM;

    std::vector<H_T> map;
    map.resize(N_ * N_);
    std::fill(map.begin(), map.end(), 0.);

    for (int32_t n0 = 0; n0 < ee0; ++n0) {
        const int32_t i1 = n0 * N_;
        const int32_t i2 = (ee1 - 1) + n0 * N_;
        map[i1] = conf.Height;
        map[i2] = conf.Height;
    }

    for (int32_t n1 = 0; n1 < ee1; ++n1) {
        const int32_t i1 = n1;
        const int32_t i2 = n1 + (ee0 - 1) * N_;
        map[i1] = conf.Height;
        map[i2] = conf.Height;
    }

    Tree tree (map, N);

    //--------------------------

    const int32_t bc = (N_ / Tree::BUCKET_SIZE);
    std::vector<bool> mm;
    mm.resize(bc * bc);
    std::fill(mm.begin(), mm.end(), true);

    int32_t et = 0;

    std::random_device rd;
    int64_t Seed = conf.UseRandomSeed ? rd() : conf.Seed;
    Rand g(Seed);
    DistD dist = DistD(conf.MinBoxVolume, conf.MaxBoxVolume);

    const auto rbox = [&](const T _vol, const bool _c) -> glm::vec<3, T> {
        const T s = std::pow(_vol, 1. / 3.);
        if (_c) return glm::vec<3, T>(s) * 0.5;
        const T ss = 3 * s;
        const T d1 = DistD(0.2 * ss, ss)(g);
        const T d2 = DistD(0.2 * (ss - d1), ss - d1)(g);
        const T d3 = ss - d1 - d2;
        return { d1, d2, d3 };
    };

    const auto b = conf->Box;
    const auto box = conf->Box->GetDefaultObject<ARandomBox>();

    std::vector<glm::vec<3, T>> next_set(conf->SetSize);
    std::deque<glm::vec<3, T>> next_q;

    if (conf->BoxType == BoxGenerationType::LIST) {
        UBoxList* list = conf->List[conf->ListIndex]->GetDefaultObject<UBoxList>();
        std::vector<glm::vec<3, T>> tmp;
        for (const FPPair& p : list->List) {
            for (int32 i = 0; i < p.Count; ++i)
                tmp.push_back(p.Size);
        }
        if (list->ShuffleBoxes)
            std::shuffle(tmp.begin(), tmp.end(), g);
        for (const auto& f : tmp)
            next_q.push_back(f);
    }

    while (true) {

        if (conf->BoxType == BoxGenerationType::LIST && next_q.empty()) break;

        std::mutex mut;
        std::vector<Detail::Result<T>> res;
        std::atomic<double> minc = std::numeric_limits<double>::infinity();

        std::atomic<int32_t> c = 0;

        next_set.clear();
        
        switch (conf->BoxType) {
            case BoxGenerationType::RANDOM:
            {
                for (int32_t i = 0; i < conf.SetSize; ++i)
                    next_set.push_back(rbox(dist(g), conf->CubeRandomBoxes));
            }
            break;
            case BoxGenerationType::LIST:
            {
                for (int32_t i = 0; i < conf.SetSize; ++i) {
                    if (next_q.empty()) break;
                    next_set.push_back(next_q.front() * 0.5);
                    next_q.pop_front();
                }
            }
            break;
        }
        
        for (size_t i = 0; i < next_set.size(); ++i) {

            const glm::vec<3, T>& nextSize = next_set[i];

            const FBox aabb = box->get_aabb(nextSize);
            c += 6;

            //Z_XY
            dispatch_impl(res, mut, minc, c, i,
                std::rint(aabb.GetExtent().X),
                std::rint(aabb.GetExtent().Y),
                nextSize,
                std::rint(aabb.GetExtent().Z),
                EAxisPerm::Z_XY_0, b, co);

            //Z_YX
            dispatch_impl(res, mut, minc, c, i,
                std::rint(aabb.GetExtent().Y),
                std::rint(aabb.GetExtent().X),
                nextSize,
                std::rint(aabb.GetExtent().Z),
                EAxisPerm::Z_XY_1, b, co);

            //Y_XZ
            dispatch_impl(res, mut, minc, c, i,
                std::rint(aabb.GetExtent().X),
                std::rint(aabb.GetExtent().Z),
                nextSize,
                std::rint(aabb.GetExtent().Y),
                EAxisPerm::Y_XZ_0, b, co);

            //Y_ZX
            dispatch_impl(res, mut, minc, c, i,
                std::rint(aabb.GetExtent().Z),
                std::rint(aabb.GetExtent().X),
                nextSize,
                std::rint(aabb.GetExtent().Y),
                EAxisPerm::Y_XZ_1, b, co);

            //X_YZ
            dispatch_impl(res, mut, minc, c, i,
                std::rint(aabb.GetExtent().Y),
                std::rint(aabb.GetExtent().Z),
                nextSize,
                std::rint(aabb.GetExtent().X),
                EAxisPerm::X_YZ_0, b, co);

            //X_ZY
            dispatch_impl(res, mut, minc, c, i,
                std::rint(aabb.GetExtent().Z),
                std::rint(aabb.GetExtent().Y),
                nextSize,
                std::rint(aabb.GetExtent().X),
                EAxisPerm::X_YZ_1, b, co);

            while (c != 0) {}

            if (!isPacking_) break;

        }

        if (!isPacking_) break;

        if (res.empty()) {
            et++;
            mcc_++;
            if (useRandomBox && mcc_ > conf->MaxEmptryTries) break;
            if (useRandomBox && et > conf->EmptryTries) break;
            if (!useRandomBox && conf->EnforceMisses && mcc_ > conf->MaxEmptryTries) break;
            if (!useRandomBox && conf->EnforceMisses && et > conf->EmptryTries) break;
            continue;
        }

        et = 0;

        std::sort(res.begin(), res.end(), [](const ::Detail::Result& _e1, const ::Detail::Result& _e2) {
            return _e1.weight < _e2.weight;
        });

        auto& r = res[0];

        for (size_t i = 0; i < next_set.size(); ++i) {
            if (i == r.set_index) continue;
            next_q.push_front(next_set[i] * 2.);
        }

        const FBox tar = FBox(
            FVector(r.n0 - r.ext.X, r.n1 - r.ext.Y, r.h),
            FVector(r.n0 + r.ext.X, r.n1 + r.ext.Y, r.h + 2 * r.ext.Z));

        r.trans = ::Detail::make_transform(
            r.perm,
            tar,
            box->get_aabb(r.ext).GetCenter(),
            box->get_relative_location()
        );

        {
            std::lock_guard<std::mutex> lock(*m_);
            toSpawn_.push_back(r);
        }

        std::fill(mm.begin(), mm.end(), false);
        for (int32_t n0 = int32_t(tar.Min.X) / Tree::BUCKET_SIZE; n0 <= int32_t(tar.Max.X) / Tree::BUCKET_SIZE; ++n0) {
            for (int32_t n1 = int32_t(tar.Min.Y) / Tree::BUCKET_SIZE; n1 <= int32_t(tar.Max.Y) / Tree::BUCKET_SIZE; ++n1) {
                const int32_t iid = n0 + n1 * bc;
                mm[iid] = true;
            }
        }

        for (int32_t n0 = int32_t(tar.Min.X); n0 <= int32_t(tar.Max.X); ++n0) {
            for (int32_t n1 = int32_t(tar.Min.Y); n1 <= int32_t(tar.Max.Y); ++n1) {
                const int32_t i = n1 + n0 * N_;
                map_[i] = tar.Max.Z;
            }
        }

        tree_->recompute(mm);

    }
    
}

//------------------------------

template<class T>
glm::mat<4, 4, T> TP::Detail::make_transform (
    const EAxisPerm _perm,
    const FBox<T>& _target,
    const glm::vec<3, T>& _pivot_offset
) {

    const glm::vec<3, T>& ctr = _target.GetCenter();

    switch (_perm) {
        case EAxisPerm::Z_XY_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr - _pivot_offset);
            return tr;
        }
        case EAxisPerm::Z_XY_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Y, _pivot_offset.X, -_pivot_offset.Z });
            return tr;
        }
        case EAxisPerm::Z_XY_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(2 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.X, _pivot_offset.Y, -_pivot_offset.Z });
            return tr;
        }
        case EAxisPerm::Z_XY_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(3 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Y, _pivot_offset.X, -_pivot_offset.Z });
            return tr;
        }
        //---------------------
        case EAxisPerm::Y_XZ_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.X, -_pivot_offset.Z, _pivot_offset.Y });
            return tr;
        }
        case EAxisPerm::Y_XZ_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Z, -_pivot_offset.X, _pivot_offset.Y });
            return tr;
        }
        case EAxisPerm::Y_XZ_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(2 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.X, _pivot_offset.Z, _pivot_offset.Y });
            return tr;
        }
        case EAxisPerm::Y_XZ_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(3 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Z, _pivot_offset.X, _pivot_offset.Y });
            return tr;
        }
        //---------------------
        case EAxisPerm::X_YZ_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Y, -_pivot_offset.Z, -_pivot_offset.X });
            return tr;
        }
        case EAxisPerm::X_YZ_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(2 * 90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Z, _pivot_offset.Y, -_pivot_offset.X });
            return tr;
        }
        case EAxisPerm::X_YZ_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(3 * 90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Y, _pivot_offset.Z, -_pivot_offset.X });
            return tr;
        }
        case EAxisPerm::X_YZ_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Z, -_pivot_offset.Y, -_pivot_offset.X });
            return tr;
        }
        //-------------------------
        //-------------------------
        case EAxisPerm::Z_n_XY_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.X, _pivot_offset.Y, -_pivot_offset.Z });
            return tr;
        }
        case EAxisPerm::Z_n_XY_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Y, _pivot_offset.X, _pivot_offset.Z });
            return tr;
        }
        case EAxisPerm::Z_n_XY_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(2 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.X, _pivot_offset.Y, _pivot_offset.Z });
            return tr;
        }
        case EAxisPerm::Z_n_XY_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(3 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Y, -_pivot_offset.X, _pivot_offset.Z });
            return tr;
        }
        //-------------------------
        case EAxisPerm::Y_n_XZ_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.X, -_pivot_offset.Z, _pivot_offset.Y });
            return tr;
        }
        case EAxisPerm::Y_n_XZ_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Z, _pivot_offset.X, _pivot_offset.Y });
            return tr;
        }
        case EAxisPerm::Y_n_XZ_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(2 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.X, _pivot_offset.Z, _pivot_offset.Y });
            return tr;
        }
        case EAxisPerm::Y_n_XZ_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180.), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(3 * 90.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Z, -_pivot_offset.X, _pivot_offset.Y });
            return tr;
        }
        //-------------------------
        case EAxisPerm::X_n_YZ_0:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180. + 90), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0.), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Y, -_pivot_offset.Z, _pivot_offset.X });
            return tr;
        }
        case EAxisPerm::X_n_YZ_1:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(-90), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(-90), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { -_pivot_offset.Z, _pivot_offset.Y, _pivot_offset.X });
            return tr;
        }
        case EAxisPerm::X_n_YZ_2:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(3 * 90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(180 + 90), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(0), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Y, _pivot_offset.Z, _pivot_offset.X });
            return tr;
        }
        case EAxisPerm::X_n_YZ_3:
        {
            glm::mat<4, 4, T> tr = glm::mat4(1.);
            tr = glm::rotate(tr, glm::radians(-90.), glm::vec3(1, 0, 0));
            tr = glm::rotate(tr, glm::radians(270), glm::vec3(0, 1, 0));
            tr = glm::rotate(tr, glm::radians(270), glm::vec3(0, 0, 1));
            tr = glm::translate(tr, ctr + { _pivot_offset.Z, -_pivot_offset.Y, _pivot_offset.X });
            return tr;
        }
    }
}//TP::Detail::make_transform