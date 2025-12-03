#pragma once

// #include "forms.h"
#include "spectral_op.hpp"

#include <fstream>
#include <memory>
#include <string>

#include <dolfinx.h>
#include <dolfinx/geometry/utils.h>
#include <dolfinx/la/Vector.h>

using namespace dolfinx;

template <typename T, int P>
class LinearSpectral2D {
public:
    LinearSpectral2D(
        basix::FiniteElement<T> element,
        mesh::Mesh<T> mesh);
};