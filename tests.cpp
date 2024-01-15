#include <iostream>
#include <cstring>
#include <cmath>
// #include "my"

int main() {
    char** a = new char*[2];
    a[0] = new char[3];
    a[1] = new char[3];
    a[0][0] = 'a';
    a[0][1] = 'b';
    a[0][2] = 'c';
    a[1][0] = 'd';
    a[1][1] = 'e';
    a[1][2] = 'f';

    for (int i=0; i<2; i++) {
        std::cout << a[i] << std::endl;
    }

    // std::cout << "Hello World!" << std::endl;
    // std::cout << "C++ Standard Version: " << __cplusplus << std::endl;
    return 0;
}
