#include <nanobind/nanobind.h>

int add(int a, int b) { return a + b; }

NB_MODULE(my_ex1, m) {
    m.def("add", &add);
}