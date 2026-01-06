#pragma once

#include "misc/util/abc_global.h"

#include "66lut_dsd.hpp"
#include "node_global.hpp"

ABC_NAMESPACE_CXX_HEADER_START

namespace acd
{

inline bool stp66_find_mx_my( const TT& root_tt, Lut66DsdResult& result )
{
  result = run_66lut_dsd_by_mx_subset( root_tt.f01, root_tt.order, /*depth_for_print=*/0 );
  return result.found;
}

} // namespace acd

ABC_NAMESPACE_CXX_HEADER_END