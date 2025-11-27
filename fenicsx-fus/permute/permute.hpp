// Copyright (C) 2022 Adeeb Arif Kor
// SPDX-License-Identifier:    MIT

#pragma once

#include <basix/finite-element.h>
#include <span>
#include <vector>

/// Reorder dofmap into tensor product order
/// @param[in] in_arr Input dofmap
/// @param[in] celltype Cell type
/// @param[in] p Degree of basis function
/// @param[out] out_arr Output dofmap
void reorder_dofmap(
  std::span<std::int32_t> out_arr,
  MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan<const std::int32_t, MDSPAN_IMPL_STANDARD_NAMESPACE::dextents<std::size_t, 2>> in_arr, 
  basix::cell::type celltype, int p);