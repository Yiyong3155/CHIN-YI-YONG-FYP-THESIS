#include <iostream>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;

std::vector<int> addv(std::vector<int> &in) {
    std::vector<int> out(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
        out[i] = in[i] + 1;
    }
    return out;
}

NB_MODULE(my_ex2, m) {
    m.def("addv", &addv);
}