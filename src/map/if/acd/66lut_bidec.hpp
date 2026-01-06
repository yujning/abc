#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <iostream>

// =====================================================
// Debug
// =====================================================
inline bool LUT66_BIDEC_DEBUG_PRINT = false;
inline int  LUT66_BIDEC_DEBUG_INDENT = 0;

#define BIDEC_DBG(msg) \
    do { if (LUT66_BIDEC_DEBUG_PRINT) { \
        for (int _i = 0; _i < LUT66_BIDEC_DEBUG_INDENT; ++_i) std::cout << "  "; \
        std::cout << msg << std::endl; \
    }} while (0)

#define BIDEC_INDENT_INC() ++LUT66_BIDEC_DEBUG_INDENT
#define BIDEC_INDENT_DEC() --LUT66_BIDEC_DEBUG_INDENT

inline uint64_t bidec_pow2(int k) { return 1ull << k; }

// =====================================================
// Result
// =====================================================
struct StrongBiDecResult
{
    bool found = false;
    std::vector<int> A;
    std::vector<int> B;
    std::vector<int> C;
    std::string MY;
    std::string MX;
};

// =====================================================
// reorder tt by variable order (MSB->LSB)
// =====================================================
inline std::string reorder_tt_by_var_order(
    const std::string& tt,
    int n,
    const std::vector<int>& new_order,
    const std::vector<int>& old_order)
{
    std::string out(tt.size(), '0');

    std::unordered_map<int,int> pos;
    for (int i = 0; i < n; ++i)
        pos[old_order[i]] = i;

    for (size_t idx = 0; idx < tt.size(); ++idx)
    {
        uint64_t old_idx = 0;
        for (int i = 0; i < n; ++i)
        {
            int bit = (idx >> (n-1-i)) & 1;
            old_idx |= uint64_t(bit) << (n-1-pos[new_order[i]]);
        }
        out[idx] = tt[old_idx];
    }
    return out;
}

// =====================================================
// index (A|B|C) → MF
// =====================================================
inline uint64_t mf_index(uint64_t a, uint64_t b, uint64_t c, int y, int z)
{
    return (a << (y + z)) | (b << z) | c;
}

// =====================================================
// Step 1: MF → MXY
// =====================================================
inline std::string compute_MXYX_from_MF(
    const std::string& MF,
    int x, int y, int z)
{
    uint64_t HN = bidec_pow2(x);
    uint64_t MN = bidec_pow2(y);
    uint64_t LN = bidec_pow2(z);

    std::string MXY(bidec_pow2(x + 2*y + z), 'x');

    for (uint64_t a = 0; a < HN; ++a)
        for (uint64_t b = 0; b < MN; ++b)
            for (uint64_t c = 0; c < LN; ++c)
            {
                uint64_t idx = (a*MN*MN + b*MN + b)*LN + c;
                MXY[idx] = MF[mf_index(a,b,c,y,z)];
            }

    return MXY;
}

// =====================================================
// Step 2: solve MX / MY from MXY
// =====================================================
inline bool solve_MX_MY_from_MXY(
    const std::string& MXY,
    int x, int y, int z,
    std::string& MX,
    std::string& MY)
{
    uint64_t blocks = bidec_pow2(x + y);
    uint64_t blk_sz = bidec_pow2(y + z);

    std::vector<std::string> classes;
    MY.assign(blocks, '0');

    for (uint64_t i = 0; i < blocks; ++i)
    {
        std::string blk = MXY.substr(i * blk_sz, blk_sz);

        bool hit = false;
        for (size_t k = 0; k < classes.size(); ++k)
        {
            bool ok = true;
            for (uint64_t j = 0; j < blk_sz; ++j)
                if (blk[j] != 'x' && classes[k][j] != 'x' && blk[j] != classes[k][j])
                { ok = false; break; }

            if (ok)
            {
                for (uint64_t j = 0; j < blk_sz; ++j)
                    if (classes[k][j] == 'x') classes[k][j] = blk[j];
                MY[i] = (k == 0 ? '1' : '0');
                hit = true;
                break;
            }
        }

        if (!hit)
        {
            if (classes.size() == 2) return false;
            classes.push_back(blk);
            MY[i] = (classes.size() == 1 ? '1' : '0');
        }
    }

    if (classes.size() != 2) return false;

    MX.clear();
    for (auto& c : classes)
        for (char b : c)
            MX.push_back(b == 'x' ? '0' : b);

    return true;
}

// =====================================================
// enumerate (x,y,z)
// =====================================================
inline std::vector<std::tuple<int,int,int>> enumerate_xyz(int n)
{
    std::vector<std::tuple<int,int,int>> out;
    for (int y = 1; y <= 4; ++y)
        for (int x = 0; x <= 6; ++x)
        {
            int z = n - x - y;
            if (z < 0) continue;
            if (x + y > 6) continue;
            if (y + z + 1 > 6) continue;
            out.emplace_back(x,y,z);
        }
    return out;
}

// =====================================================
// enumerate partitions A,B,C
// =====================================================
inline std::vector<
    std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>>
enumerate_variable_partitions(
    const std::vector<int>& vars,
    int x,int y,int z)
{
    std::vector<
        std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>> out;

    int n = vars.size();
    std::vector<bool> sel_b(n,false);
    std::fill(sel_b.begin(), sel_b.begin()+y, true);

    do {
        std::vector<int> B, rest;
        for (int i=0;i<n;i++) (sel_b[i]?B:rest).push_back(vars[i]);

        std::vector<bool> sel_a(rest.size(),false);
        std::fill(sel_a.begin(), sel_a.begin()+x, true);

        do {
            std::vector<int> A,C;
            for (size_t i=0;i<rest.size();i++)
                (sel_a[i]?A:C).push_back(rest[i]);
            out.emplace_back(A,B,C);
        } while (std::prev_permutation(sel_a.begin(), sel_a.end()));

    } while (std::prev_permutation(sel_b.begin(), sel_b.end()));

    return out;
}

// =====================================================
// main entry (delay-aware)
// =====================================================
inline StrongBiDecResult
run_strong_bi_dec(const std::string& tt,
                  const std::vector<int>& order,
                  uint32_t delay_profile = 0)
{
    StrongBiDecResult res;
    int n = order.size();
    if ((size_t(1)<<n) != tt.size()) return res;

    auto is_late = [&](int v){ return (delay_profile >> v) & 1; };

    int total_late = 0;
    for (int v : order) if (is_late(v)) total_late++;

    for (auto [x,y,z] : enumerate_xyz(n))
    {
        auto parts = enumerate_variable_partitions(order,x,y,z);
        for (auto const& tpl : parts)
        {
            const auto& A = std::get<0>(tpl);
            const auto& B = std::get<1>(tpl);
            const auto& C = std::get<2>(tpl);

            bool bad = false;
            for (int v : A)
                if (is_late(v)) { bad = true; break; }
            if (bad) continue;

            int lateC = 0;
            for (int v : C) if (is_late(v)) lateC++;
            if (lateC < std::min((int)C.size(), total_late)) continue;

            std::vector<int> new_order;
            new_order.insert(new_order.end(), A.begin(), A.end());
            new_order.insert(new_order.end(), B.begin(), B.end());
            new_order.insert(new_order.end(), C.begin(), C.end());

            std::string MFp = reorder_tt_by_var_order(tt,n,new_order,order);
            std::string MXY = compute_MXYX_from_MF(MFp,x,y,z);

            std::string MX, MY;
            if (!solve_MX_MY_from_MXY(MXY,x,y,z,MX,MY)) continue;

            res.found = true;
            res.A=A; res.B=B; res.C=C;
            res.MX=MX; res.MY=MY;
            return res;
        }
    }
    return res;
}
