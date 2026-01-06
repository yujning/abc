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

// =====================================================
// Result
// =====================================================
struct Lut66DsdResult {
    bool found = false;
    size_t L = 0;
    std::string Mx;
    std::string My;

    std::vector<int> mx_pos;
    std::vector<int> my_pos;

    std::vector<int> mx_vars_msb2lsb;
    std::vector<int> my_vars_msb2lsb;
};

// =====================================================
// next combination
// =====================================================
inline bool next_combination_66(std::vector<int>& comb, int m)
{
    int k = comb.size();
    for (int i=k-1;i>=0;--i)
        if (comb[i] < m-k+i)
        {
            ++comb[i];
            for (int j=i+1;j<k;++j) comb[j]=comb[j-1]+1;
            return true;
        }
    return false;
}

// =====================================================
// extract block
// =====================================================
inline std::string extract_block_for_mx_66(
    const std::string& mf,
    int n,
    const std::vector<int>& mx_pos,
    const std::vector<int>& my_pos,
    uint64_t my_assignment)
{
    int k = mx_pos.size();
    int m = my_pos.size();
    std::string sub(1ull<<k,'0');

    for (uint64_t mx=0;mx<(1ull<<k);++mx)
    {
        uint64_t idx=0;
        for (int i=0;i<m;i++)
            idx |= ((my_assignment>>(m-1-i))&1ull) << (n-1-my_pos[i]);
        for (int i=0;i<k;i++)
            idx |= ((mx>>(k-1-i))&1ull) << (n-1-mx_pos[i]);
        sub[mx]=mf[idx];
    }
    return sub;
}

// =====================================================
// main entry (delay-aware)
// =====================================================
inline Lut66DsdResult
run_66lut_dsd_by_mx_subset(
    const std::string& mf,
    const std::vector<int>& order,
    uint32_t delay_profile,
    int depth_for_print = 0)
{
    Lut66DsdResult out;
    int n = order.size();
    if ((size_t(1)<<n) != mf.size()) return out;

    auto is_late = [&](int v){ return (delay_profile>>v)&1; };

    struct VarPos{int var;int pos;};
    std::vector<VarPos> vp;
    for (int i=0;i<n;i++) vp.push_back({order[i],i});

    // ★ slow vars first → MX
    std::sort(vp.begin(),vp.end(),
        [&](auto&a,auto&b){
            bool la=is_late(a.var), lb=is_late(b.var);
            if (la!=lb) return la;
            return a.var<b.var;
        });

    int k_min = std::max(1,n-6);
    int k_max = std::min(5,n-1);

    for (int k=k_max;k>=k_min;--k)
    {
        std::vector<int> comb(k);
        for (int i=0;i<k;i++) comb[i]=i;

        while (true)
        {
            std::vector<int> mx_pos, my_pos;
            std::vector<char> used(n,0);

            for (int id:comb)
            {
                mx_pos.push_back(vp[id].pos);
                used[vp[id].pos]=1;
            }
            for (int i=0;i<n;i++) if (!used[i]) my_pos.push_back(i);

            // forbid slow vars in MY
            bool bad=false;
            for (int p:my_pos)
                if (is_late(order[p])) { bad=true; break; }
            if (bad)
            {
                if (!next_combination_66(comb,vp.size())) break;
                continue;
            }

            int m=n-k;
            if (m>6)
            {
                if (!next_combination_66(comb,vp.size())) break;
                continue;
            }

            std::unordered_map<std::string,int> cls;
            std::vector<std::string> blocks;
            std::string My;
            bool fail=false;

            for (uint64_t y=0;y<(1ull<<m);++y)
            {
                std::string blk=extract_block_for_mx_66(mf,n,mx_pos,my_pos,y);
                auto it=cls.find(blk);
                if (it==cls.end())
                {
                    if (blocks.size()==2){fail=true;break;}
                    int id=blocks.size();
                    cls[blk]=id;
                    blocks.push_back(blk);
                    My.push_back(id==0?'1':'0');
                }
                else My.push_back(it->second==0?'1':'0');
            }

            if (!fail && blocks.size()==2)
            {
                out.found=true;
                out.L=1ull<<k;
                out.Mx=blocks[0]+blocks[1];
                out.My=My;

                out.mx_pos=mx_pos;
                out.my_pos=my_pos;
                for (int p:mx_pos) out.mx_vars_msb2lsb.push_back(order[p]);
                for (int p:my_pos) out.my_vars_msb2lsb.push_back(order[p]);
                return out;
            }

            if (!next_combination_66(comb,vp.size())) break;
        }
    }
    return out;
}
