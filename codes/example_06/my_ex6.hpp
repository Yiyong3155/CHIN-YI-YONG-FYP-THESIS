#pragma once

#include <iostream>

template <typename T>
class Compute {
public:
    Compute(T a_);

    void run();

private:
    T a;
};