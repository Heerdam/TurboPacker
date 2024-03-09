
#include <TP.hpp>

void taskgraphtest() {

    using namespace TP;
    using namespace TP::Detail;

    std::mutex m;
    int32_t i = 0;

    const auto f = [&]() {
        std::lock_guard<std::mutex> lock(m);
        std::cout << i++ << std::endl;
    };

    TaskGraph tg (4);

    for(int32_t i = 0; i < 10; ++i)
        tg.dispatch(f);

    tg.wait();

}

int main() {

    //taskgraphtest();
    //return 0;

    using namespace TP;

    Config<float, CostFunction::CF_Basic> conf;
    conf.Bounds = {80., 100. };
    conf.Height = 100.;
    conf.EmptryTries = 1;
    conf.MinBoxVolume = 20 * 20 * 20;
    conf.MaxBoxVolume = 30 * 30 * 30;

    const auto pr = solve(conf);

    while(!pr.isDone()){
        std::cout << pr.getBoxCount() << " [" << pr.getPackDensity() << "]"  << "\r";
    }

    std::cout << "done" << std::endl;
    std::cout << pr.data().size() << std::endl;

    return 0;
}