# TurboPacker
The TurboPacker is a high level library built upon the the MedianQuadTree (MQT) backend. It allows for simple integration 
into existing rendering front ends or any other application. It facilitates a rich interface for to pack boxes of various kind into one or multiple bins. <br />
This library was part of my master thesis. For a exhaustive discussion of the TurboPacker I refere to here.

## Quick start
For use as a library add this project with cmake or simply copy the TP.hpp header. Clone the whole project for use as App.  <br />

### Minimal example for library integration
```
#include <TP.hpp>
using namespace TP;
int main() {
    BinInfo<float> bin;
    bin.Bounds = {100., 100.}; 
    bin.Height = 100.;
    
    Config<float, CostFunction::CF_Krass> config;
    config.Bins.push_back(bin);

    auto promise = solve(config);
    promise.wait();
    assert(promise.done());
    const auto result = promise.data();
    //...   
    return 0;
}
```

### Minimal example for custom costfunction
```
#include <TP.hpp>
using namespace TP;

template<class T>
struct CF_Example {
    static float eval(const Detail::Result<T>& _r) {
        const auto nrml = _r.normalize();
        return nrml.n0 * nrml.n1;
    }
}

int main() {
    //...
    Config<float, CF_Example> config;
    //...
}
```
