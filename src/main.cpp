
#include <TP.hpp>

int main() {

    using namespace TP;

    Config<float, CostFunction::CF_Basic> conf;
    conf.Bounds = {80., 100. };

    const auto pr = solve(conf);

    while(!pr.isDone()){
        std::cout << pr.getPackDensity() << std::endl;
    }

    std::cout << "done" << std::endl;


    return 0;
}