#include <iostream>
#include <algorithm>
#include <vector>


int main() {
    std::cout << "C++ Standard Version: " << __cplusplus << std::endl;
    return 0;
}


// int main() {

//     AnIntWrap wrap;
//     wrap.val = 3;

//     MyObj obj(wrap);
//     wrap.val = 4;
//     obj.print();
//     // std::vector<int> numbers = {3, 2, 8, 5, 1, 4, 7, 6};

//     // Using std::nth_element to find the median (middle element)
//     // auto middle = numbers.begin() + numbers.size() / 2;
//     // int find = 3;
//     // auto middle = numbers.begin() + find;
//     // std::nth_element(numbers.begin(), numbers.begin()+find, numbers.end());

//     // std::cout << "The median is: " << *middle << std::endl;

//     // // The rest of the vector is not sorted
//     // for (const auto& num : numbers) {
//     //     std::cout << num << " ";
//     // }

//     return 0;
// }




/*
almost works (22 linking errors) in Developer Command Prompt for VS 2022
cl /EHsc /arch:AVX2 /I C:/Python38/include /I C:/Python38/Lib/site-packages/numpy/core/include .\tests2.cpp /link /LIBPATH:C:/Python38/libs python3.lib /SUBSYSTEM:CONSOLE




*/