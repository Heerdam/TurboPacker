
#include <TP.hpp>
#include <Util.hpp>

template<class T>
struct CF_Annealer {

    inline static CF_Annealer<T>* inst_;

    TP::Detail::Annealer<T, true, false> an_;

    CF_Annealer() = default;
    void init(const uint64_t _seed, const size_t _weights_count) {
        an_ = TP::Detail::Annealer<T, true, false>(_seed, _weights_count);
    }

    static float eval(const TP::Detail::Result<T>& _r) { 
        const auto nrml = _r.normalize();
        return 
            inst_->an_[0].first * std::pow(nrml.bin, inst_->an_[0].second) +
            inst_->an_[1].first * std::pow(nrml.n0, inst_->an_[1].second) +
            inst_->an_[2].first * std::pow(nrml.n1, inst_->an_[2].second) +
            inst_->an_[3].first * std::pow(nrml.h, inst_->an_[3].second) +
            inst_->an_[4].first * std::pow(nrml.l, inst_->an_[4].second) +
            inst_->an_[5].first * std::pow(nrml.m, inst_->an_[5].second) +
            inst_->an_[6].first * std::pow(nrml.b_l, inst_->an_[6].second) +
            inst_->an_[7].first * std::pow(nrml.b_m, inst_->an_[7].second) +
            inst_->an_[9].first * std::pow(nrml.b_h, inst_->an_[8].second);
    };
};



//------------------------

int main() {

    using namespace TP;

    CF_Annealer<float>::inst_ = new CF_Annealer<float>();
    CF_Annealer<float>::inst_->init(1234567890, 10);

    auto conf = Disk::create_default<float, CF_Annealer, uint16_t, uint32_t, 15, std::allocator<uint16_t>>();
    conf.BoxType = Detail::BoxGenerationType::LIST;
    conf.Bins[0].Bounds = { 100, 100 };
    conf.Bins[0].Height = 100;
    conf.NumThreads = 30;
    conf.LookAheadSize = 5;

    auto conf_ref = Disk::create_default<float, CostFunction::CF_Krass, uint16_t, uint32_t, 15, std::allocator<uint16_t>>();
    conf_ref.BoxType = conf.BoxType;
    conf_ref.Bins[0].Bounds = conf.Bins[0].Bounds;
    conf_ref.Bins[0].Height = conf.Bins[0].Height;
    conf_ref.NumThreads = conf.NumThreads;
    conf_ref.LookAheadSize = conf.LookAheadSize;

    //---------------------------------------

    {
        std::cout << "----- CF Krass -----" << std::endl;

        float mean = 0.f;
        float mmin = 1.f;
        float mmax = 0.f;
        for(int32_t ss = 0; ss < 100; ++ss) {
            conf_ref.Seed = std::random_device()();
            auto prms = solve(conf_ref);
            prms.wait();
            mean += prms.getTotalPackDensity();
            mmin = std::min(prms.getTotalPackDensity(), mmin);
            mmax = std::max(prms.getTotalPackDensity(), mmax);
            std::cout << std::format("Step {}: Time: {}s, Density: {}, Count: {}\r", ss, prms.getTime(), prms.getTotalPackDensity(), prms.getTotalBoxCount());
            std::cout.flush();
        }
        mean /= 100;
        std::cout.flush();
        std::cout << std::endl << std::format("Result: {}, Variance: {}", mean, mmax - mmin) << std::endl;
    }

    //---------------------------------------
    {
        std::vector<float> dens;
        float last_dens = 0.f;
        int32_t stuck = 0;
        std::vector<std::pair<float, float>> bweights;

        std::cout << "----- Annealing -----" << std::endl;

        for(int32_t ss = 0; ss < 30; ++ss) {
            conf.Seed = std::random_device()();

            for(int32_t s = 0; s < 25; ++s) {
            
                auto prms = solve(conf);
                prms.wait();

                std::cout << std::format("Step {}, {}: Time: {}s, Density: {}, Count: {}, Best: {}", ss, s, prms.getTime(), prms.getTotalPackDensity(), prms.getTotalBoxCount(), last_dens);
                
                if(stuck > 10){ //punch through
                    stuck = 0;
                    CF_Annealer<float>::inst_->an_.step(9, 0.2f);
                    std::cout << " PT \r";
                } else if(last_dens > prms.getTotalPackDensity()){ //try different direction
                    stuck++;
                    CF_Annealer<float>::inst_->an_.reverse();
                    CF_Annealer<float>::inst_->an_.step(9, 0.015f);
                    std::cout << " S \r";
                } else {
                    dens.push_back(last_dens);
                    last_dens = prms.getTotalPackDensity();
                    bweights = CF_Annealer<float>::inst_->an_.cpy();
                    CF_Annealer<float>::inst_->an_.step(9, 0.015f);
                    std::cout << " N \r";
                }

                std::cout.flush();

            }
        }

        std::cout.flush();
        std::cout << std::format("Gain: {}, {}", dens[1], dens.back()) << std::endl;
        std::cout << std::endl << "----- Validating -----" << std::endl;
        CF_Annealer<float>::inst_->an_.set(bweights);

        float mean = 0.f;
        float mmin = 1.f;
        float mmax = 0.f;
        for(int32_t ss = 0; ss < 100; ++ss) {
            conf.Seed = std::random_device()();
            auto prms = solve(conf);
            prms.wait();
            mean += prms.getTotalPackDensity();
            mmin = std::min(prms.getTotalPackDensity(), mmin);
            mmax = std::max(prms.getTotalPackDensity(), mmax);
            std::cout << std::format("Step {}: Time: {}s, Density: {}, Count: {}\r", ss, prms.getTime(), prms.getTotalPackDensity(), prms.getTotalBoxCount());
            std::cout.flush();
        }
        mean /= 100;
        std::cout.flush();
        std::cout << std::endl << std::format("Result: {}, Expected: {}, Variance: {}", mean, last_dens, mmax - mmin) << std::endl;

    }

    return EXIT_SUCCESS;
}