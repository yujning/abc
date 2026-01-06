/**C++File**************************************************************

  FileName    [ac_wrapper.cpp]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Ashenhurst-Curtis decomposition.]

  Synopsis    [Interface with the FPGA mapping package.]

  Author      [Alessandro Tempia Calvino]
  
  Affiliation [EPFL]

  Date        [Ver. 1.0. Started - November 20, 2023.]

***********************************************************************/

#include "ac_wrapper.h"
#include "ac_decomposition.hpp"
#include "acd66.hpp"
#include "acdXX.hpp"
#include "stp66.hpp"
#include <cstdio>
ABC_NAMESPACE_IMPL_START

int acd_evaluate( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
{
  using namespace acd;

  ac_decomposition_params ps;
  ps.lut_size = lutSize;
  ps.use_first = false;
  ps.try_no_late_arrival = static_cast<bool>( try_no_late_arrival );
  ac_decomposition_stats st;

  ac_decomposition_impl acd( nVars, ps, &st );
  int val = acd.run( pTruth, *pdelay );

  if ( val < 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();
  *cost = st.num_luts;

  return val;
}

int acd_decompose( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned char *decomposition )
{
  using namespace acd;

  ac_decomposition_params ps;
  ps.lut_size = lutSize;
  ps.use_first = true;
  ac_decomposition_stats st;

  ac_decomposition_impl acd( nVars, ps, &st );
  acd.run( pTruth, *pdelay );
  int val = acd.compute_decomposition();

  if ( val < 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();
  acd.get_decomposition( decomposition );
  return 0;
}

int acd2_evaluate( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
{
  using namespace acd;

  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );
  int val = acd.run( pTruth, *pdelay );

  if ( val == 0 )
  {
    *pdelay = 0;
    return -1;
  }

  acd.compute_decomposition();
  *pdelay = acd.get_profile();
  *cost = 2;

  return val;
}

int acd2_decompose( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned char *decomposition )
{
  using namespace acd;

  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );
  acd.run( pTruth, *pdelay );
  int val = acd.compute_decomposition();

  if ( val != 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();

  acd.get_decomposition( decomposition );
  return 0;
}

inline int acd66_evaluate( word * pTruth, unsigned nVars )
{
  using namespace acd;

  acd66_impl acd( nVars, true, false );

  if ( acd.run( pTruth ) == 0 )
    return 0;

  return 1;
}

int acdXX_evaluate( word * pTruth, unsigned lutSize, unsigned nVars )
{
  using namespace acd;

  if ( lutSize == 6 )
  {
    return acd66_evaluate( pTruth, nVars );
  }
  
  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );

  if ( acd.run( pTruth ) == 0 )
    return 0;

  return 1;
}

int acdXX_decompose( word * pTruth, unsigned lutSize, unsigned nVars, unsigned char *decomposition )
{
  using namespace acd;

  acdXX_params ps;
  ps.lut_size = lutSize;

  for ( int i = 0; i <= lutSize - 2; ++i )
  {
    ps.max_shared_vars = i;
    ps.min_shared_vars = i;
    acdXX_impl acd( nVars, ps );

    if ( acd.run( pTruth ) == 0 )
      continue;
    acd.compute_decomposition();
    acd.get_decomposition( decomposition );
    return 0;
  }

  return 1;
}

static inline std::string truth_to_string( word* pTruth, unsigned nVars )
{
  const size_t size = size_t( 1 ) << nVars;
  std::string tt( size, '0' );

  for ( size_t i = 0; i < size; ++i )
  {
    size_t word_id = i >> 6;
    size_t bit_id = i & 63;
    tt[i] = ( ( pTruth[word_id] >> bit_id ) & 1ull ) ? '1' : '0';
  }

  return tt;
}

static inline uint64_t truth_string_to_uint64( const std::string& tt )
{
  uint64_t value = 0;
  for ( size_t i = 0; i < tt.size(); ++i )
  {
    if ( tt[i] == '1' )
      value |= ( 1ull << i );
  }
  return value;
}

int stpxx_decompose( word * pTruth, unsigned nVars, unsigned lutSize, unsigned *pdelay, unsigned char *decomposition )
{
  using namespace acd;

  (void)pdelay;

  if ( lutSize != 6 )
  {
    std::printf( "[STP66] Skip: unsupported LUT size %u, falling back to default decomposition\n", lutSize );
    return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
  }

  if ( nVars == 0 || nVars > 11 )
  {
    std::printf( "[STP66] Skip: cut size %u is outside supported range [1, 11], falling back to default decomposition\n", nVars );
    return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
  }

  TT root_tt;
  root_tt.f01 = truth_to_string( pTruth, nVars );
  root_tt.order.reserve( nVars );
  for ( int v = static_cast<int>( nVars ) - 1; v >= 0; --v )
    root_tt.order.push_back( v );

  Lut66DsdResult res;
  if ( !stp66_find_mx_my( root_tt, res ) )
  {
    std::printf( "[STP66] No 6+6 decomposition found, falling back to default decomposition\n" );
    return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
  }

}

ABC_NAMESPACE_IMPL_END
