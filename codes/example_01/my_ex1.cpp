#include <iostream>
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;
using namespace nb::literals;

int add(int a, int b = 1) { return a + b; }
int mult(int a, int b) { return a * b; }

class FinalYearProject {
public:
    FinalYearProject(
        std::string student_name
    )
    : name(student_name) {
    }

    void print_name() const {
        std::cout << "The student name is " << name << "\n";
    }

private:
    std::string name;
};

NB_MODULE(my_ex1, m) {
    m.def("add", &add, "a"_a, "b"_a = 1,
          "This function adds two numbers and increments if only one is provided.");
    m.def("mult", &mult);
    
    nb::class_<FinalYearProject>(m, "FinalYearProject")
        .def(nb::init<std::string &>())
        .def("print_name", &FinalYearProject::print_name);

    m.doc() = "A simple example of nanobind";
}