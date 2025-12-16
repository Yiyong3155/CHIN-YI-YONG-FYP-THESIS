#include <iostream>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

struct Cat { };
const char *meow(Cat *cat) {
    return cat != nullptr ? "meow!" : "(no cat)";
}


NB_MODULE(my_ex4, m) {
    m.def("double", [](float x) { return 2.f * x; }); // lambda function
    nb::class_<Cat>(m, "Cat")
        .def(nb::init<>());
        m.def("meow", &meow, nb::arg("cat").none());
}