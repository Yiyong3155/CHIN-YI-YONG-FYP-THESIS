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
        std::string student_name,
        std::string student_matric_no,
        int student_year
    )
    : name(student_name),
      matric_no(student_matric_no),
      year(student_year) {
    }

    void print_info() const {
        std::cout << "The student name is " << name << "\n";
        std::cout << "The student matric no is " << matric_no << "\n";
        std::cout << "The student year is " << year << "\n";
    }

private:
    std::string name;
    std::string matric_no;
    int year;
};

NB_MODULE(my_ex1, m) {
    m.def("add", &add, "a"_a, "b"_a = 1,
          "This function adds two numbers and increments if only one is provided.");
    m.def("mult", &mult);
    
    nb::class_<FinalYearProject>(m, "FinalYearProject")
        .def(nb::init<std::string&, std::string&, int>())
        .def("print_info", &FinalYearProject::print_info);

    m.doc() = "A simple example of nanobind";
}