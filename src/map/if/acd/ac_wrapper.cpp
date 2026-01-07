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
extern "C" int stpxx_decompose(
    word * pTruth,
    unsigned nVars,
    int lutSize,
    unsigned *pdelay,
    unsigned char *decomposition );

extern "C" int stpxx_evaluate( word * pTruth, unsigned nVars, int lutSize,
                    unsigned *pdelay, unsigned *cost, int try_no_late_arrival );
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
// =====================================================
// ★ STP66 evaluate（ACD 同签名/同语义）
//   - 成功：返回 level = 1(有late) or 2(无late)
//   - 失败：返回 -1（cut 判死，避免 direct TT node 产生 >6 LUT）
//   - try_no_late_arrival == 1 时：失败会再尝试把 delay_profile 置 0 再跑一次（模仿 ACD）
// =====================================================

int stpxx_evaluate( word * pTruth, unsigned nVars, int lutSize,
                    unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
{
  using namespace acd;

  if ( lutSize != 6 )
    return -1;
  if ( nVars == 0 || nVars > 11 )
    return -1;

  // IF 传进来的 *pdelay 是 “late-arrival mask”(uLeafMask)
  uint32_t dp = *pdelay;

  auto try_run = [&](uint32_t delay_profile)->bool {
    TT root_tt;
    root_tt.f01 = truth_to_string( pTruth, nVars );
    root_tt.order.clear();
    root_tt.order.reserve( nVars );
    for ( int v = (int)nVars - 1; v >= 0; --v )
      root_tt.order.push_back( v );

    Lut66DsdResult res;
    if ( !stp66_find_mx_my( root_tt, delay_profile, res ) || !res.found )
      return false;

    // 你现在 STP 输出固定 2-LUT 结构
    *cost = 2;
    *pdelay = delay_profile; // 关键：告诉 IF 本次用的 delay_profile（可能被置 0）
    return true;
  };

  if ( try_run(dp) )
  {
    // ACD 的 run：delay_profile==0 -> 2 levels，否则 1 level
    return (dp == 0) ? 2 : 1;
  }

  // 模仿 ACD：如果带 dp 不行，且允许 try_no_late_arrival，则再试 dp=0
  if ( try_no_late_arrival )
  {
    if ( try_run(0) )
      return 2; // dp=0 -> 2 levels
  }

  // 不可分解：cut 判死（非常关键，否则就会走 direct truth table node，出现 10-LUT）
  *pdelay = 0;
  return -1;
}



static inline bool write_lut_tt_bytes_safe(
    unsigned char*& p,
    unsigned char* end,
    const std::string& tt_bits,
    unsigned m /* #fanin */)
{
  // tt_bits size must be 2^m
  if (tt_bits.size() != (size_t(1) << m)) return false;

  unsigned num_bytes = (m <= 3) ? 1u : (1u << (m - 3));
  if (p + num_bytes > end) return false;

  // pack bits into bytes (LSB-first)
  for (unsigned byte = 0; byte < num_bytes; ++byte)
  {
    unsigned char v = 0;
    for (unsigned b = 0; b < 8; ++b)
    {
      unsigned idx = byte * 8 + b;
      if (idx < tt_bits.size() && tt_bits[idx] == '1')
        v |= (unsigned char)(1u << b);
    }
    *p++ = v;
  }
  return true;
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

static inline int encode_2lut_decomposition_abc_safe(
    unsigned nVars,
    const Lut66DsdResult& res,
    unsigned char* decompArray,
    unsigned decompCapBytes = 92 )
{
  if (!decompArray || decompCapBytes < 2) return -1;

  const unsigned m1 = (unsigned)res.my_vars_msb2lsb.size(); // fanin of MY
  const unsigned k  = (unsigned)res.mx_vars_msb2lsb.size(); // mx vars (without MY)
  const unsigned m2 = 1u + k;                               // fanin of MX

  if (m1 < 1 || m1 > 6) return -1;
  if (m2 < 1 || m2 > 6) return -1;

  if (res.My.size() != (size_t(1) << m1)) return -1;

  const size_t L = size_t(1) << k;
  if (res.Mx.size() != 2 * L) return -1;

  const std::string block0 = res.Mx.substr(0, L);
  const std::string block1 = res.Mx.substr(L, L);
  const std::string mx_tt  = build_mux_tt_from_blocks(block0, block1);
  if (mx_tt.size() != (size_t(1) << m2)) return -1;

  unsigned char* p   = decompArray;
  unsigned char* end = decompArray + decompCapBytes;

  // header: [numBytes][numLuts]
  if (end - p < 2) return -1;
  *p++ = 0;       // placeholder numBytes
  *p++ = 2;       // 2 LUTs

  // ---- LUT0: MY ----
  if (p >= end) return -1;
  *p++ = (unsigned char)m1;

  // fanins list must be LSB->MSB of the LUT truth table.
  // Your res.my_vars_msb2lsb is MSB->LSB, so we output reversed.
  for (int i = (int)m1 - 1; i >= 0; --i)
  {
    int abc_idx = res.my_vars_msb2lsb[i];   // ✅ now 0-based already
    if (abc_idx < 0 || abc_idx >= (int)nVars) return -1;
    if (p >= end) return -1;
    *p++ = (unsigned char)abc_idx;
  }

  if (!write_lut_tt_bytes_safe(p, end, res.My, m1)) return -1;

  const int my_node_idx = (int)nVars; // internal node index in ABC record

  // ---- LUT1: MX ----
  if (p >= end) return -1;
  *p++ = (unsigned char)m2;

  for (int i = (int)k - 1; i >= 0; --i)
  {
    int abc_idx = res.mx_vars_msb2lsb[i];   // ✅ now 0-based already
    if (abc_idx < 0 || abc_idx >= (int)nVars) return -1;
    if (p >= end) return -1;
    *p++ = (unsigned char)abc_idx;
  }

  if (p >= end) return -1;
  *p++ = (unsigned char)my_node_idx;

  if (!write_lut_tt_bytes_safe(p, end, mx_tt, m2)) return -1;

  // finalize numBytes
  unsigned bytes = (unsigned)(p - decompArray);
  if (bytes > decompCapBytes) return -1;
  decompArray[0] = (unsigned char)bytes;
  return 0;
}

// =====================================================
// ★ stpxx_decompose：STP66 主入口
// =====================================================

int stpxx_decompose(
    word * pTruth,
    unsigned nVars,
    int lutSize,
    unsigned *pdelay,
    unsigned char *decomposition )

{
  using namespace acd;

  if ( lutSize != 6 )
    //return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
      {
    *pdelay = 0;
    return -1;
  }

  if ( nVars == 0 || nVars > 11 )
//return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
  {
    *pdelay = 0;
    return -1;
  }


  TT root_tt;
  root_tt.f01 = truth_to_string( pTruth, nVars );
  root_tt.order.reserve( nVars );
  for ( int v = (int)nVars - 1; v >= 0; --v )
    root_tt.order.push_back( v ); // 1-based

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
    if ( encode_2lut_decomposition_abc_safe( nVars, res, decomposition, 92 ) == 0 )

    {
       *pdelay = delay_profile;
      std::printf(
        "[STP66] Exported 2-LUT decomposition: |MY|=%zu, |MX|=%zu\n",
        res.my_vars_msb2lsb.size(),
        1 + res.mx_vars_msb2lsb.size() );
            
      return 0;
    }
  }

  //return acd_decompose( pTruth, nVars, lutSize, pdelay, decomposition );
   *pdelay = 0;
  return -1;
}
int stpXX_evaluate_simple( word * pTruth, unsigned lutSize, unsigned nVars )
{
    if ( lutSize != 6 )
        return 0;

    if ( nVars == 0 || nVars > 11 )
        return 0;

    TT root_tt;
    root_tt.f01 = truth_to_string( pTruth, nVars );
    root_tt.order.clear();
    root_tt.order.reserve( nVars );
    for ( int v = (int)nVars - 1; v >= 0; --v )
        root_tt.order.push_back( v );

    Lut66DsdResult res;
    if ( acd::stp66_find_mx_my( root_tt, /*delay_profile=*/0, res ) && res.found )
        return 1;

    return 0;
}

ABC_NAMESPACE_IMPL_END
