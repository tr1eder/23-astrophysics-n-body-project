// #include "Python.h"
#include <iostream>
#include <vector>
#include <cstdint>
typedef int64_t __int64;

#include "matplotlibcpp.h"

#include <cstdint>

// Change __int64 to int64_t
// int64_t myVariable;


// #include "mytest"

namespace plt = matplotlibcpp;

int main() {
    // double a = 3.915;
    // int b = a;
    // std::cout << b << std::endl;
    std::vector<double> x = {1, 2, 3, 4};
    std::vector<double> y = {1, 3, 2, 4};

    plt::plot(x, y);
    plt::show();

    return 0;
}




// ${workspaceFolder}/**
// C:/cygwin64/usr/include
// C:/Python38/include
// C:/Python38/Lib/site-packages/numpy/core/include
// C:/Python38/libs

// g++ tests2.cpp -o tests2.exe -I C:/Python38/include -I C:/Python38/Lib/site-packages/numpy/core/include -L C:/Python38/libs -lpython38 