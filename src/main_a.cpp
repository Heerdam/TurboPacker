
#include <TP.hpp>
#include <Util.hpp>
#include <fstream>

template<class T>
struct CF_Annealer {

    inline static CF_Annealer<T>* inst_;

    TP::Detail::Annealer<T, true> an_;

    CF_Annealer() = default;
    void init(const uint64_t _seed, const size_t _weights_count) {
        an_ = TP::Detail::Annealer<T, true>(_seed, _weights_count);
    }

    static float eval(const TP::Detail::Result<T>& _r) { 
        //const auto nrml = _r.normalize();
        float out = 0.f;
        for(size_t i = 0; i < inst_->an_.get_size(); ++i){
            out += inst_->an_[i];
        }
        return out;
        // return out + 
        //     inst_->an_[0].first * std::pow(nrml.bin, inst_->an_[0].second) +
        //     inst_->an_[1].first * std::pow(nrml.n0, inst_->an_[1].second) +
        //     inst_->an_[2].first * std::pow(nrml.n1, inst_->an_[2].second) +
        //     inst_->an_[3].first * std::pow(nrml.h, inst_->an_[3].second) +
        //     inst_->an_[4].first * std::pow(nrml.l, inst_->an_[4].second) +
        //     inst_->an_[5].first * std::pow(nrml.m, inst_->an_[5].second) +
        //     inst_->an_[6].first * std::pow(nrml.b_l, inst_->an_[6].second) +
        //     inst_->an_[7].first * std::pow(nrml.b_m, inst_->an_[7].second) +
        //     inst_->an_[9].first * std::pow(nrml.b_h, inst_->an_[8].second);
    };
};



//------------------------

int main() {

    using namespace TP;

    constexpr int32_t w_count = 500;
    constexpr int32_t w_step = 500;
    constexpr float w_size = 0.1f;
    constexpr float w_min = 0.2f;
    constexpr float w_max = 0.3f;

    CF_Annealer<float>::inst_ = new CF_Annealer<float>();
    CF_Annealer<float>::inst_->init(1234567890, w_count);

    const auto sseed = std::random_device()();

    auto conf = Disk::create_default<float, CF_Annealer, uint16_t, uint32_t, 15, std::allocator<uint16_t>>();
    conf.BoxType = Detail::BoxGenerationType::RANDOM;
    conf.CubeRandomBoxes = false;
    conf.Bins[0].Bounds = { 100, 100 };
    conf.Bins[0].Height = 100;
    conf.NumThreads = 30;
    conf.LookAheadSize = 1;
    conf.MinBoxVolume = std::pow(5, 3);
    conf.MaxBoxVolume = std::pow(15, 3);

    auto conf_ref = Disk::create_default<float, CostFunction::CF_Krass, uint16_t, uint32_t, 15, std::allocator<uint16_t>>();
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

    if constexpr(true){
        std::cout << "----- CF Krass -----" << std::endl;

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

        file << std::format("CF Krass: Its: {}, Mean: {}, Variance: {}, Bounds: {}, {}", 150, mean, var, mmin, mmax) << std::endl;
    }

    //---------------------------------------
    {
        std::vector<float> dens;
        float last_dens = 0.f;
        int32_t stuck = 0;
        std::vector<float> bweights;

        std::cout << "----- Annealing -----" << std::endl;

        const uint64_t aseed = std::random_device()();
        std::mt19937_64 arand (aseed);

        for(int32_t ss = 0; ss < 15; ++ss) {
            conf.Seed = arand();

            for(int32_t s = 0; s < 25; ++s) {
            
                auto prms = solve(conf);
                prms.wait();

                std::cout << std::format("Step {}, {}: Time: {}s, Density: {}, Count: {}, Best: {}", ss, s, prms.getTime(), prms.getTotalPackDensity(), prms.getTotalBoxCount(), last_dens);
                
                if(stuck > 10){ //punch through
                    stuck = 0;
                    CF_Annealer<float>::inst_->an_.punch_through(w_min, w_max);
                    std::cout << " PT                 \r";
                } else if(last_dens > prms.getTotalPackDensity()){ //try different direction
                    stuck++;
                    CF_Annealer<float>::inst_->an_.reverse();
                    CF_Annealer<float>::inst_->an_.step(w_step, w_size);
                    std::cout << " S                \r";
                } else {
                    dens.push_back(last_dens);
                    last_dens = prms.getTotalPackDensity();
                    bweights = CF_Annealer<float>::inst_->an_.cpy();
                    CF_Annealer<float>::inst_->an_.step(w_step, w_size);
                    std::cout << " N                \r";
                }

            }
        }

        std::cout << std::endl << std::format("Gain: {}, {}", dens[1], dens.back()) << std::endl;
        file << std::format("Annealing: Seed: {}, Gain: {}, {}", aseed, dens[1], dens.back()) << std::endl;
        std::cout << std::endl << "----- Validating -----" << std::endl;
        CF_Annealer<float>::inst_->an_.set(bweights);

        std::mt19937_64 rand (sseed);

        float mean = 0.f;
        float mmin = 1.f;
        float mmax = 0.f;
        std::vector<float> vls(100);
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
        std::cout << std::endl << std::format("Result: {}, Expected: {}, Variance: {}", mean, last_dens, var) << std::endl;
        file << std::format("Validation: Its: {}, Mean: {}, Variance: {}, Bounds: {}, {}", 150, mean, var, mmin, mmax) << std::endl;

    }

    file << "-------------------------------" << std::endl;

    std::cout << "---------------------------" << std::endl;

    return EXIT_SUCCESS;
}