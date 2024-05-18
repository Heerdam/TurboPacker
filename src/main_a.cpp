
#include <TP.hpp>
#include <Util.hpp>
#include <fstream>
#include <algorithm>
#include <numeric>

template<class T>
struct CF_Annealer {

    inline static CF_Annealer<T>* inst_;

    TP::Detail::Annealer<T, true> an_;

    CF_Annealer() = default;
    void init(const uint64_t _seed, const size_t _weights_count) {
        an_ = TP::Detail::Annealer<T, true>(_seed, _weights_count);
    }

    static float eval(const TP::Detail::Result<T>& _r) { 
        const auto nrml = _r.normalize();
        // constexpr size_t f = 9;
        // const size_t idx = (_r.n1 + _r.n0 * _r.ext1_) * f;  

        // return inst_->an_[idx] + inst_->an_[idx + 1] + 
        //     inst_->an_[idx + 2] * nrml.h + 
        //     inst_->an_[idx + 3] * nrml.l + 
        //     inst_->an_[idx + 4] * nrml.m + 
        //     inst_->an_[idx + 5] * nrml.h + 
        //     inst_->an_[idx + 6] * nrml.b_l + 
        //     inst_->an_[idx + 7] * nrml.b_m + 
        //     inst_->an_[idx + 8] * nrml.b_h;

        return
            inst_->an_[0] * (nrml.bin) +
            inst_->an_[1] * (nrml.n0) +
            inst_->an_[2] * (nrml.n1) +
            inst_->an_[3] * (nrml.h) +
            inst_->an_[4] * (nrml.l) +
            inst_->an_[5] * (nrml.m) +
            inst_->an_[6] * (nrml.b_l) +
            inst_->an_[7] * (nrml.b_m) +
            inst_->an_[9] * (nrml.b_h);
    };
};

//------------------------

template<class T>
struct CF_Custom_1 {
    static T eval(const TP::Detail::Result<T>& _r) {
        constexpr float w[] = {0.01426f, -0.00987833f, -0.0163784f, -0.0160077f, 0.f, 0.00836435f, -0.00417864f, 0.000154695f, -0.000847495f, -0.0046661f};
        const auto nrml = _r.normalize();
        return 
            w[0] * (nrml.bin) +
            w[1] * (nrml.n0) +
            w[2] * (nrml.n1) +
            w[3] * (nrml.h) +
            w[4] * (nrml.l) +
            w[5] * (nrml.m) +
            w[6] * (nrml.b_l) +
            w[7] * (nrml.b_m) +
            w[9] * (nrml.b_h);
    };
};//CF_Krass

//------------------------

int main() {

    using namespace TP;

    constexpr int32_t w_count = 10; //100*100*9;
    constexpr int32_t w_step = w_count;
    constexpr float w_size = 0.05f;
    constexpr float w_min = 0.3f;
    constexpr float w_max = 0.25;

    CF_Annealer<float>::inst_ = new CF_Annealer<float>();
    CF_Annealer<float>::inst_->init(1234567890, w_count);

    const auto sseed = std::random_device()();

    auto conf = Disk::create_default<float, CF_Annealer, uint16_t, uint32_t, 15, std::allocator<uint16_t>>();
    conf.BoxType = Detail::BoxGenerationType::LIST;
    conf.CubeRandomBoxes = false;
    conf.Bins[0].Bounds = { 100, 100 };
    conf.Bins[0].Height = 100;
    conf.NumThreads = 8;
    conf.LookAheadSize = 25;
    conf.MinBoxVolume = std::pow(5, 3);
    conf.MaxBoxVolume = std::pow(15, 3);

    auto conf_ref = Disk::create_default<float, TP::CostFunction::CF_Krass, uint16_t, uint32_t, 15, std::allocator<uint16_t>>();
    conf_ref.BoxType = conf.BoxType;
    conf_ref.CubeRandomBoxes = conf.CubeRandomBoxes;
    conf_ref.Bins[0].Bounds = conf.Bins[0].Bounds;
    conf_ref.Bins[0].Height = conf.Bins[0].Height;
    conf_ref.NumThreads = conf.NumThreads;
    conf_ref.LookAheadSize = conf.LookAheadSize;
    conf_ref.MinBoxVolume = conf.MinBoxVolume;
    conf_ref.MaxBoxVolume = conf.MaxBoxVolume;

    //---------------------------------------

    std::cout << "---------------------------" << std::endl;

    std::ofstream file("../res.txt", std::ios::app);
    file << std::format("Type: {}", (int32_t(conf.BoxType) ? "List" : "Random")) << std::endl;
    file << std::format("Bounds: {} x {} x {}", conf.Bins[0].Bounds.x, conf.Bins[0].Bounds.y, conf.Bins[0].Height) << std::endl;
    file << std::format("Look ahead: {}", conf.LookAheadSize) << std::endl;
    if(conf.BoxType == Detail::BoxGenerationType::RANDOM) {
        file << std::boolalpha << std::format("Cube: {}", conf.CubeRandomBoxes) << std::endl;
        file << std::format("Volume: {}, {}", conf.MinBoxVolume, conf.MaxBoxVolume) << std::endl;
    }
    file << std::format("Validation Seed: {}", sseed) << std::endl;
    file << std::format("Weights: {}, Step: {}, Type: {}, Size: {}, Min: {}, Max: {}", CF_Annealer<float>::inst_->an_.get_size(), w_step, "exp", w_size, w_min, w_max) << std::endl << std::endl;

    if constexpr(true){ //Mean: 0.30294323, Variance: 0.09177472
        std::cout << "----- CF_Krass -----" << std::endl;

        std::mt19937_64 rand (sseed);

        float mean = 0.f;
        float mmin = 1.f;
        float mmax = 0.f;
        std::vector<float> vls(150);
        for(int32_t ss = 0; ss < 150; ++ss) {
            conf_ref.Seed = rand();
            auto prms = solve(conf_ref);
            prms.wait();
            vls.push_back(prms.getTotalPackDensity());
            mean += prms.getTotalPackDensity();
            mmin = std::min(prms.getTotalPackDensity(), mmin);
            mmax = std::max(prms.getTotalPackDensity(), mmax);
            std::cout << std::format("Step {}: Time: {}s, Density: {}, Count: {}                     \r", ss, prms.getTime(), prms.getTotalPackDensity(), prms.getTotalBoxCount());
        }
        mean /= 150;
        float var = 0.f;
        for(int32_t ss = 0; ss < 150; ++ss) {
            var += std::pow(vls[ss] - mean, 2);
        }
        var /= 150;
        std::cout << std::endl << std::format("Result: {}, Variance: {}", mean, var) << std::endl;

        file << std::format("CF_Krass: Its: {}, Mean: {}, Variance: {} ({}), Bounds: {}, {}", 150, mean, var, std::sqrt(var), mmin, mmax) << std::endl;
    }

    //---------------------------------------
    if constexpr(true){

        std::cout << "----- Annealing -----" << std::endl;

        const uint64_t aseed = std::random_device()();
        std::mt19937_64 arand (aseed);

        std::vector<float> s_weights = CF_Annealer<float>::inst_->an_.cpy();
        float s_var = 1.f;
        float s_mean = 0.f;

        std::vector<float> last_weights = CF_Annealer<float>::inst_->an_.cpy();

        //stochastic gradient descent
        for(int32_t steps = 0; steps < 10; steps++){


            //sample region
            std::vector<float> b_weights = CF_Annealer<float>::inst_->an_.cpy();
            float b_mean = 0.f;
            float b_var = 1.f;
            for(int32_t r = 0; r < 5; r++){

                CF_Annealer<float>::inst_->an_.step(w_step, w_size);

                //sample current weights
                constexpr int32_t smpls = 10;
                std::vector<float> vls(smpls * 3);
                for(int32_t ss = 0; ss < smpls; ++ss) {

                    conf.Seed = arand();
                    auto prms1 = solve(conf);
                    conf.Seed = arand();
                    auto prms2 = solve(conf);
                    conf.Seed = arand();
                    auto prms3 = solve(conf);

                    prms1.wait();
                    prms2.wait();
                    prms3.wait();

                    vls.push_back(prms1.getTotalPackDensity());
                    vls.push_back(prms2.getTotalPackDensity());
                    vls.push_back(prms3.getTotalPackDensity());

                    std::cout << r << " - " << ss << "                                                                                \r";

                }

                const float mean = (float)std::accumulate(vls.begin(), vls.end(), 0.f) / float(smpls * 3);
                const float var = (float)std::accumulate(vls.begin(), vls.end(), 0.f, [&](float _a, float _b) -> float{
                        return _a + std::pow(_b - mean, 2);
                    }) / float(smpls * 3);

                //const float d_mean = std::abs(mean - b_mean);
                //const float d_var = std::abs(var - b_var);

                if(mean > b_mean && var < b_var){
                    b_weights = CF_Annealer<float>::inst_->an_.cpy();
                    b_mean = mean;
                    b_var = var;
                }

                std::cout << std::format("{}, {}, Best: {}, {}                        ", mean, var, b_mean, b_var) << std::endl;
   
                CF_Annealer<float>::inst_->an_.reverse();
            }

            std::cout << std::format("Step: {}, Mean: {}, Var: {}, Best: {}, {}                               ", steps, b_mean, b_var, s_mean, s_var) << std::endl;

            if(b_mean > s_mean && b_var < s_var){
                s_weights = b_weights;
                s_mean = b_mean;
                s_var = b_var;
            }

            //gradient to next region
            std::vector<float> gd;
            gd.resize(w_count);

            for(size_t i = 0; i < w_count; ++i){
                gd[i] = (b_weights[i] - last_weights[i]);
            }

            const float frac = std::accumulate(gd.begin(), gd.end(), 0.f);
            const float delta = std::uniform_real_distribution<float>(w_min, w_max)(arand);
            for(size_t i = 0; i < w_count; ++i){
                gd[i] /= frac;
                gd[i] *= delta;
                gd[i] += b_weights[i];
            }

            CF_Annealer<float>::inst_->an_.set(gd);
            last_weights = gd;

        }

        std::cout << std::format("Annealing: Seed: {}, Mean: {}, Var: {}", aseed, s_mean, s_var) << std::endl;
        file << std::format("Annealing: Seed: {}, Mean: {}, Var: {}", aseed, s_mean, s_var) << std::endl;
        file << "Weights: ";
        for(float w : s_weights)
            file << w << " ";
        file << std::endl;
        
        //----------------------------------

        std::cout << std::endl << "----- Validating -----" << std::endl;
        CF_Annealer<float>::inst_->an_.set(s_weights);
        std::mt19937_64 rand (sseed);

        float mean = 0.f;
        float mmin = 1.f;
        float mmax = 0.f;
        std::vector<float> vls(150);
        for(int32_t ss = 0; ss < 150; ++ss) {
            conf.Seed = rand();
            auto prms = solve(conf);
            prms.wait();
            mean += prms.getTotalPackDensity();
            mmin = std::min(prms.getTotalPackDensity(), mmin);
            mmax = std::max(prms.getTotalPackDensity(), mmax);
            vls.push_back(prms.getTotalPackDensity());
            std::cout << std::format("Step {}: Time: {}s, Density: {}, Count: {}                          \r", ss, prms.getTime(), prms.getTotalPackDensity(), prms.getTotalBoxCount());
        }
        mean /= 150;
        float var = 0.f;
        for(int32_t ss = 0; ss < 150; ++ss) {
            var += std::pow(vls[ss] - mean, 2);
        }
        var /= 150;
        std::cout << std::endl << std::format("Result: {}, Expected: {}, {}, Variance: {}", mean, s_mean, s_var, var) << std::endl;
        file << std::format("Validation: Its: {}, Mean: {}, Variance: {}, Bounds: {}, {}", 150, mean, var, mmin, mmax) << std::endl;

        std::cout << std::endl;

    }

    file << "-------------------------------" << std::endl;

    std::cout << "---------------------------" << std::endl;

    return EXIT_SUCCESS;
}