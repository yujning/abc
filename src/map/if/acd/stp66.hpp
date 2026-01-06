#pragma once

#include "misc/util/abc_global.h"

#include "66lut_dsd.hpp"
#include "node_global.hpp"
#include "66lut_bidec.hpp"

ABC_NAMESPACE_CXX_HEADER_START

namespace acd
{

inline bool stp66_find_mx_my( const TT& root_tt, Lut66DsdResult& result )
{
  result = run_66lut_dsd_by_mx_subset( root_tt.f01, root_tt.order, /*depth_for_print=*/0 );
   if ( result.found )
    return true;

  // 若 DSD 失败，尝试强双分解算法
  const int n = static_cast<int>( root_tt.order.size() );
  if ( ( size_t( 1 ) << n ) != root_tt.f01.size() )
    return false;

  const std::vector<int> original_order = root_tt.order;

  for ( auto [x, y, z] : enumerate_xyz( n ) )
  {
    auto parts = enumerate_variable_partitions( original_order, x, y, z );
    for ( const auto& [A, B, C] : parts )
    {
      std::vector<int> new_order;
      new_order.insert( new_order.end(), A.begin(), A.end() );
      new_order.insert( new_order.end(), B.begin(), B.end() );
      new_order.insert( new_order.end(), C.begin(), C.end() );

      std::string MFp = reorder_tt_by_var_order( root_tt.f01, n, new_order, original_order );
      std::string MXY = compute_MXYX_from_MF( MFp, x, y, z );

      std::string MX, MY;
      if ( !solve_MX_MY_from_MXY( MXY, x, y, z, MX, MY ) )
        continue;

      result.found = true;
      result.Mx = MX;
      result.My = MY;
      result.block0.clear();
      result.block1.clear();
      result.reordered_tt.clear();

      result.mx_vars_msb2lsb.clear();
      result.mx_vars_msb2lsb.insert( result.mx_vars_msb2lsb.end(), B.begin(), B.end() );
      result.mx_vars_msb2lsb.insert( result.mx_vars_msb2lsb.end(), C.begin(), C.end() );

      result.my_vars_msb2lsb.clear();
      result.my_vars_msb2lsb.insert( result.my_vars_msb2lsb.end(), A.begin(), A.end() );
      result.my_vars_msb2lsb.insert( result.my_vars_msb2lsb.end(), B.begin(), B.end() );

      result.L = size_t( 1 ) << ( y + z );

      auto fill_pos = [&]( const std::vector<int>& vars, std::vector<int>& pos_out ) {
        pos_out.clear();
        for ( int var : vars )
        {
          auto it = std::find( original_order.begin(), original_order.end(), var );
          if ( it != original_order.end() )
            pos_out.push_back( static_cast<int>( std::distance( original_order.begin(), it ) ) );
        }
      };

      fill_pos( result.mx_vars_msb2lsb, result.mx_pos );
      fill_pos( result.my_vars_msb2lsb, result.my_pos );

      return true;
    }
  }

  return false;
}

} // namespace acd

ABC_NAMESPACE_CXX_HEADER_END