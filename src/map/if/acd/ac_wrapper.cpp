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
#include <string>

ABC_NAMESPACE_IMPL_START
std::string bin_to_hex(const std::string& bin)
{
    std::string b = bin;

    // 左侧补 0，使长度是 4 的倍数
    size_t pad = (4 - b.size() % 4) % 4;
    b.insert(0, pad, '0');

    std::string hex;
    hex.reserve(b.size() / 4);

    for (size_t i = 0; i < b.size(); i += 4) {
        int v = (b[i]   - '0') << 3 |
                (b[i+1] - '0') << 2 |
                (b[i+2] - '0') << 1 |
                (b[i+3] - '0');
        hex.push_back("0123456789ABCDEF"[v]);
    }

    return hex;
}
static void stp66_print_delay_profile( unsigned nVars, uint32_t delay_profile )
{
  std::printf( "[STP66] delay_profile=0x%08x, late_vars(abc_idx):",
               delay_profile );
  bool printed = false;
  for ( unsigned v = 0; v < nVars; ++v )
  {
    if ( delay_profile & ( 1u << v ) )
    {
      std::printf( " %u", v );
      printed = true;
    }
  }
  if ( !printed )
  {
    std::printf( " none" );
  }
  std::printf( "\n" );
}

// =====================================================
// 原始 ACD 接口（保持不变）
// =====================================================

int acd_evaluate( word * pTruth, unsigned nVars, int lutSize,
                  unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
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

int acd_decompose( word * pTruth, unsigned nVars, int lutSize,
                   unsigned *pdelay, unsigned char *decomposition )
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

// =====================================================
// XX ACD（保持不变）
// =====================================================

int acd2_evaluate( word * pTruth, unsigned nVars, int lutSize,
                   unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
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

int acd2_decompose( word * pTruth, unsigned nVars, int lutSize,
                    unsigned *pdelay, unsigned char *decomposition )
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

// =====================================================
// 66 only evaluate
// =====================================================

inline int acd66_evaluate( word * pTruth, unsigned nVars )
{
  using namespace acd;
  acd66_impl acd( nVars, true, false );
  return acd.run( pTruth ) == 0 ? 0 : 1;
}

int acdXX_evaluate( word * pTruth, unsigned lutSize, unsigned nVars )
{
  using namespace acd;

  if ( lutSize == 6 )
    return acd66_evaluate( pTruth, nVars );

  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );
  return acd.run( pTruth ) == 0 ? 0 : 1;
}

int acdXX_decompose( word * pTruth, unsigned lutSize, unsigned nVars,
                     unsigned char *decomposition )
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

// =====================================================
// 工具：truth table <-> string
// =====================================================

static inline std::string truth_to_string( word* pTruth, unsigned nVars )
{
  const size_t size = size_t(1) << nVars;
  std::string tt( size, '0' );

  for ( size_t i = 0; i < size; ++i )
  {
    size_t word_id = i >> 6;
    size_t bit_id  = i & 63;
    tt[i] = ( ( pTruth[word_id] >> bit_id ) & 1ull ) ? '1' : '0';
  }
  return tt;
}

// =====================================================
// ★ STP66 → ABC 编码辅助函数
// =====================================================

static inline void write_lut_tt_bytes(
    unsigned char*& p,
    const std::string& tt_bits )
{
  uint64_t v = 0;
  for ( size_t i = 0; i < tt_bits.size(); ++i )
    if ( tt_bits[i] == '1' )
      v |= ( 1ull << i );

  unsigned m = 0;
  while ( ( 1u << m ) < tt_bits.size() ) ++m;

  unsigned num_bytes = ( m <= 3 ) ? 1u : ( 1u << ( m - 3 ) );
  if ( num_bytes > 8 ) num_bytes = 8;

  for ( unsigned j = 0; j < num_bytes; ++j )
    *p++ = (unsigned char)( ( v >> ( 8 * j ) ) & 0xFF );
}

static inline std::string build_mux_tt_from_blocks(
    const std::string& block0,
    const std::string& block1 )
{
  const size_t L = block0.size();
  std::string out( 2 * L, '0' );

  for ( size_t i = 0; i < 2 * L; ++i )
  {
    const int sel = (int)( i >= L );
    const size_t mx = i & ( L - 1 );
    out[i] = sel ? block0[mx] : block1[mx];
  }
  return out;
}

static inline int encode_2lut_decomposition_abc(
    unsigned nVars,
    const Lut66DsdResult& res,
    unsigned char* decompArray )
{
  const unsigned m1 = (unsigned)res.my_vars_msb2lsb.size();
  const unsigned k  = (unsigned)res.mx_vars_msb2lsb.size();
  const unsigned m2 = 1u + k;

  if ( m1 == 0 || m1 > 6 || m2 == 0 || m2 > 6 )
    return -1;

  if ( res.My.size() != (size_t(1) << m1) )
    return -1;

  const size_t L = size_t(1) << k;
  if ( res.Mx.size() != 2 * L )
    return -1;

  const std::string block0 = res.Mx.substr( 0, L );
  const std::string block1 = res.Mx.substr( L, L );
  const std::string mx_tt  = build_mux_tt_from_blocks( block0, block1 );

  unsigned char* p = decompArray;
  unsigned bytes = 2;

  p++;
  *p++ = 2; // two LUTs

  // ---- LUT0: MY ----
  *p++ = (unsigned char)m1; bytes++;
  for ( int i = (int)m1 - 1; i >= 0; --i )
  {
    int abc_idx = res.my_vars_msb2lsb[i] - 1;
    if ( abc_idx < 0 || abc_idx >= (int)nVars ) return -1;
    *p++ = (unsigned char)abc_idx; bytes++;
  }
  write_lut_tt_bytes( p, res.My );
  bytes += ( m1 <= 3 ) ? 1u : ( 1u << ( m1 - 3 ) );

  const int my_node_idx = (int)nVars;

  // ---- LUT1: MX ----
  *p++ = (unsigned char)m2; bytes++;
  for ( int i = (int)k - 1; i >= 0; --i )
  {
    int abc_idx = res.mx_vars_msb2lsb[i] - 1;
    if ( abc_idx < 0 || abc_idx >= (int)nVars ) return -1;
    *p++ = (unsigned char)abc_idx; bytes++;
  }
  *p++ = (unsigned char)my_node_idx; bytes++;

  write_lut_tt_bytes( p, mx_tt );
  bytes += ( m2 <= 3 ) ? 1u : ( 1u << ( m2 - 3 ) );

  decompArray[0] = (unsigned char)bytes;
  return 0;
}

// =====================================================
// ★ stpxx_decompose：STP66 主入口
// =====================================================

int stpxx_decompose(
    word * pTruth,
    unsigned nVars,
    unsigned lutSize,
    unsigned *pdelay,
    unsigned char *decomposition )
{
  using namespace acd;

  if ( lutSize != 6 )
    return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );

  if ( nVars == 0 || nVars > 11 )
    return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );

  TT root_tt;
  root_tt.f01 = truth_to_string( pTruth, nVars );
  root_tt.order.reserve( nVars );
  for ( int v = (int)nVars - 1; v >= 0; --v )
    root_tt.order.push_back( v + 1 ); // 1-based

  const uint32_t delay_profile = *pdelay;
  stp66_print_delay_profile( nVars, delay_profile );
  std::printf( "[STP66] truth table (original order): %s\n", root_tt.f01.c_str() );
  std::string hex = bin_to_hex(root_tt.f01);

  std::printf(
      "[STP66] truth table (hex): 0x%s\n",
      hex.c_str()
  );
  
  Lut66DsdResult res;
  if ( stp66_find_mx_my( root_tt, delay_profile, res ) && res.found )
  {
    if ( encode_2lut_decomposition_abc( nVars, res, decomposition ) == 0 )
    {
       *pdelay = delay_profile;
      std::printf(
        "[STP66] Exported 2-LUT decomposition: |MY|=%zu, |MX|=%zu\n",
        res.my_vars_msb2lsb.size(),
        1 + res.mx_vars_msb2lsb.size() );
      return 0;
    }
  }

  return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
}

ABC_NAMESPACE_IMPL_END
