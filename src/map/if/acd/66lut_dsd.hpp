#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <iostream>

// =====================================================
// Debug switch
// =====================================================
inline bool LUT66_DSD_DEBUG_PRINT = false;

inline void print_tt_with_order_66(
    const std::string& title,
    const std::string& tt,
    const std::vector<int>& order,
    int depth = 0)
{
    if (!LUT66_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "📌 " << title << "\n";
    std::cout << indent << "   TT    = " << tt << "\n";
    std::cout << indent << "   order = { ";
    for (int v : order) std::cout << v << " ";
    std::cout << "}\n";
}

// =====================================================
// 66-LUT DSD result
// =====================================================
struct Lut66DsdResult {
    bool found = false;
    size_t L = 0;              // = 2^{|Mx|}
    std::string Mx;            // = block0 + block1, length 2*L
    std::string My;            // length 2^{|My|}

    // 0-based positions in current order
    std::vector<int> mx_pos;
    std::vector<int> my_pos;

    // var IDs in MSB->LSB order (by position)
    std::vector<int> mx_vars_msb2lsb;
    std::vector<int> my_vars_msb2lsb;

    // optional debug
    std::string block0;
    std::string block1;
    std::string reordered_tt;  // concatenated blocks by My assignment (My|Mx)
};

// =====================================================
// next_combination: comb is strictly increasing indices in [0..m-1]
// =====================================================
inline bool next_combination_66(std::vector<int>& comb, int m)
{
    int k = (int)comb.size();
    for (int i = k - 1; i >= 0; --i) {
        if (comb[i] < m - k + i) {
            ++comb[i];
            for (int j = i + 1; j < k; ++j) {
                comb[j] = comb[j - 1] + 1;
            }
            return true;
        }
    }
    return false;
}

// =====================================================
// Extract block (Mx subfunction) for a fixed My assignment
// - order positions: 0..n-1 (0 is MSB position)
// - mx_pos_sorted / my_pos_sorted must be sorted ascending by position (MSB->LSB)
// - my_assignment encoded in My's MSB->LSB order
// - mx_index encoded in Mx's MSB->LSB order
// =====================================================
inline std::string extract_block_for_mx_66(
    const std::string& mf,
    int n,
    const std::vector<int>& mx_pos_sorted,
    const std::vector<int>& my_pos_sorted,
    uint64_t my_assignment)
{
    int k = (int)mx_pos_sorted.size();
    int m = (int)my_pos_sorted.size();

    const size_t sub_size = 1ull << k;
    std::string sub(sub_size, '0');

    for (size_t mx_index = 0; mx_index < sub_size; ++mx_index)
    {
        uint64_t full_index = 0;

        // 1) place My bits (MSB->LSB)
        for (int t = 0; t < m; ++t)
        {
            int pos = my_pos_sorted[t];
            int bit = (int)((my_assignment >> (m - 1 - t)) & 1ull);
            full_index |= (uint64_t(bit) << (n - 1 - pos));
        }

        // 2) place Mx bits (MSB->LSB)
        for (int t = 0; t < k; ++t)
        {
            int pos = mx_pos_sorted[t];
            int bit = (int)((mx_index >> (k - 1 - t)) & 1ull);
            full_index |= (uint64_t(bit) << (n - 1 - pos));
        }

        sub[mx_index] = mf[(size_t)full_index];
    }

    return sub;
}

inline void print_candidate_info_66(
    int depth,
    int k,
    int m,
    const std::vector<int>& mx_vars_msb2lsb,
    const std::vector<int>& my_vars_msb2lsb)
{
    if (!LUT66_DSD_DEBUG_PRINT) return;

    std::string indent((size_t)depth * 2, ' ');
    std::cout << indent << "🔎 66-DSD 尝试 |My|=" << m << " (<=6), |Mx|=" << k << " (<=5)\n";
    std::cout << indent << "   My={ ";
    for (int v : my_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}  Mx={ ";
    for (int v : mx_vars_msb2lsb) std::cout << v << " ";
    std::cout << "}\n";
}

// =====================================================
// ★ 66-LUT Strong DSD (subset enumeration, single split)
// Constraint: |My| <= 6, |Mx| <= 5, and EXACTLY 2 distinct blocks
// Search: iterate k within feasible range; enumerate Mx subset by var_id ascending.
// =====================================================
inline Lut66DsdResult run_66lut_dsd_by_mx_subset(
    const std::string& mf,
    const std::vector<int>& order,
    int depth_for_print = 0)
{
    Lut66DsdResult out;

    int n = (int)order.size();
    if (n <= 1) return out;
    if ((size_t(1) << n) != mf.size()) return out;

    struct VarPos { int var; int pos; };
    std::vector<VarPos> vp;
    vp.reserve(n);
    for (int pos = 0; pos < n; ++pos) vp.push_back({order[pos], pos});
    std::sort(vp.begin(), vp.end(), [](const VarPos& x, const VarPos& y){
        return x.var < y.var;
    });

    int k_min = std::max(1, n - 6);
    int k_max = std::min(5, n - 1);
    if (k_min > k_max) return out;

    for (int k = k_max; k >= k_min; --k)
    {
        std::vector<int> comb(k);
        for (int i = 0; i < k; ++i) comb[i] = i;

        while (true)
        {
            std::vector<int> mx_pos;
            mx_pos.reserve(k);
            for (int idx : comb) mx_pos.push_back(vp[idx].pos);

            std::vector<int> my_pos;
            my_pos.reserve(n - k);
            {
                std::vector<char> is_mx(n, 0);
                for (int p : mx_pos) is_mx[p] = 1;
                for (int p = 0; p < n; ++p) if (!is_mx[p]) my_pos.push_back(p);
            }

            std::sort(mx_pos.begin(), mx_pos.end());
            std::sort(my_pos.begin(), my_pos.end());

            std::vector<int> mx_vars_msb2lsb;
            std::vector<int> my_vars_msb2lsb;
            mx_vars_msb2lsb.reserve(mx_pos.size());
            my_vars_msb2lsb.reserve(my_pos.size());
            for (int p : mx_pos) mx_vars_msb2lsb.push_back(order[p]);
            for (int p : my_pos) my_vars_msb2lsb.push_back(order[p]);

            int m = n - k;
            if (m > 6 || k > 5) {
                if (!next_combination_66(comb, (int)vp.size())) break;
                continue;
            }

            print_candidate_info_66(depth_for_print, k, m, mx_vars_msb2lsb, my_vars_msb2lsb);

            const uint64_t my_cnt = 1ull << m;
            const uint64_t L = 1ull << k;

            if (LUT66_DSD_DEBUG_PRINT)
            {
                std::string reordered;
                reordered.reserve((size_t)my_cnt * (size_t)L);

                for (uint64_t y = 0; y < my_cnt; ++y)
                    reordered += extract_block_for_mx_66(mf, n, mx_pos, my_pos, y);

                std::vector<int> reordered_order;
                reordered_order.reserve(n);
                reordered_order.insert(reordered_order.end(), my_vars_msb2lsb.begin(), my_vars_msb2lsb.end());
                reordered_order.insert(reordered_order.end(), mx_vars_msb2lsb.begin(), mx_vars_msb2lsb.end());

                print_tt_with_order_66("候选 split 的重排 TT (My|Mx)", reordered, reordered_order, depth_for_print);
            }

            std::unordered_map<std::string, int> block_map;
            std::vector<std::string> blocks;
            blocks.reserve(2);

            std::string My;
            My.reserve((size_t)my_cnt);

            bool too_many = false;

            for (uint64_t y = 0; y < my_cnt; ++y)
            {
                std::string block = extract_block_for_mx_66(mf, n, mx_pos, my_pos, y);

                auto it = block_map.find(block);
                if (it == block_map.end())
                {
                    if (blocks.size() >= 2) { too_many = true; break; }
                    int id = (int)blocks.size();
                    block_map.emplace(block, id);
                    blocks.push_back(block);
                    My.push_back(id == 0 ? '1' : '0');
                }
                else
                {
                    My.push_back(it->second == 0 ? '1' : '0');
                }
            }

            if (!too_many && blocks.size() == 2)
            {
                out.found = true;
                out.L = (size_t)L;
                out.Mx = blocks[0] + blocks[1];
                out.My = My;

                out.mx_pos = mx_pos;
                out.my_pos = my_pos;
                out.mx_vars_msb2lsb = mx_vars_msb2lsb;
                out.my_vars_msb2lsb = my_vars_msb2lsb;

                out.block0 = blocks[0];
                out.block1 = blocks[1];

                std::string reordered;
                reordered.reserve((size_t)my_cnt * (size_t)L);
                for (uint64_t y = 0; y < my_cnt; ++y)
                    reordered += extract_block_for_mx_66(mf, n, mx_pos, my_pos, y);
                out.reordered_tt = reordered;

                if (LUT66_DSD_DEBUG_PRINT)
                {
                    std::string indent((size_t)depth_for_print * 2, ' ');
                    std::cout << indent << "✅ 命中 66-LUT Strong DSD split\n";
                    std::cout << indent << "   block0 = " << blocks[0] << "\n";
                    std::cout << indent << "   block1 = " << blocks[1] << "\n";
                    print_tt_with_order_66("当前 split 的 My", My, my_vars_msb2lsb, depth_for_print);
                }

                return out;
            }

            if (!next_combination_66(comb, (int)vp.size())) break;
        }
    }

    return out;
}
