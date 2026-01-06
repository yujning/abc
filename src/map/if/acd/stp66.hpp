#pragma once

#include "misc/util/abc_global.h"

#include "node_global.hpp"   // 提供 TT
#include "66lut_dsd.hpp"     // 66-LUT DSD
#include "66lut_bidec.hpp"   // ★ 现在可以安全 include 了
#include <cstdint>
#include <iostream>

ABC_NAMESPACE_CXX_HEADER_START

namespace acd
{

inline bool stp66_find_mx_my(
    const TT& root_tt,
    uint32_t delay_profile,
    Lut66DsdResult& result )
{
    LUT66_DSD_DEBUG_PRINT   = true;
    LUT66_BIDEC_DEBUG_PRINT = true;

    std::cout << "[STP66] Try 66-LUT DSD\n";

    result = run_66lut_dsd_by_mx_subset(
                root_tt.f01, root_tt.order, delay_profile );

    bool meaningful_dsd =
        result.found &&
        !result.mx_vars_msb2lsb.empty() &&
        !result.my_vars_msb2lsb.empty() &&
        result.mx_vars_msb2lsb.size() < root_tt.order.size() &&
        result.my_vars_msb2lsb.size() < root_tt.order.size();

    if ( meaningful_dsd )
    {
        std::cout << "[STP66] meaningful 66-LUT DSD succeeded\n";
        return true;
    }

    std::cout << "[STP66] DSD trivial or failed, try strong bi-dec\n";

    StrongBiDecResult bi =
            run_strong_bi_dec(root_tt.f01, root_tt.order, delay_profile);

    if ( !bi.found )
    {
        std::cout << "[STP66] strong bi-dec failed\n";
        return false;
    }

    // map bi-dec → unified result
    result = Lut66DsdResult{};
    result.found = true;
    result.Mx = bi.MX;
    result.My = bi.MY;

    result.my_vars_msb2lsb = bi.A;
    result.my_vars_msb2lsb.insert(
        result.my_vars_msb2lsb.end(),
        bi.B.begin(), bi.B.end());

    result.mx_vars_msb2lsb = bi.B;
    result.mx_vars_msb2lsb.insert(
        result.mx_vars_msb2lsb.end(),
        bi.C.begin(), bi.C.end());

    std::cout << "[STP66] strong bi-dec succeeded\n";
    return true;
}

} // namespace acd

ABC_NAMESPACE_CXX_HEADER_END
