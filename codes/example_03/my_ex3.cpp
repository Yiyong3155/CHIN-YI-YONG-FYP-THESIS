#include <iostream>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>

namespace nb = nanobind;

std::vector<int> addv(const std::vector<int> in){
    std::vector<int> out(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
        out[i] = in[i] + 1;
    }
    return out;
}

void addvr(std::vector<int> &v){
    for (int &value: v)
        value *= 2;
}

NB_MODULE(my_ex3, m) {
    nb::bind_vector<std::vector<int>>(m, "IntVector");
    m.def("addv", &addv);
    m.def("addvr", &addvr);
}
