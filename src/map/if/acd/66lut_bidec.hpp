#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <iostream>

// =====================================================
// Debug infrastructure (aligned with 66lut_dsd style)
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

inline std::string vec2str(const std::vector<int>& v)
{
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i)
    {
        s += std::to_string(v[i]);
        if (i + 1 < v.size()) s += ",";
    }
    s += "}";
    return s;
}

// =====================================================
// 工具
// =====================================================
inline uint64_t pow2(int k) { return 1ull << k; }

// =====================================================
// 结果结构体
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
// 变量重排
// =====================================================
inline std::string reorder_tt_by_var_order(
    const std::string& tt,
    int n,
    const std::vector<int>& new_order_msb2lsb,
    const std::vector<int>& old_order_msb2lsb)
{
    std::unordered_map<int,int> pos;
    for (int i = 0; i < n; ++i)
        pos.emplace(old_order_msb2lsb[i], i);

    std::string out(tt.size(), '0');

    for (size_t idx = 0; idx < tt.size(); ++idx)
    {
        uint64_t old_idx = 0;
        for (int i = 0; i < n; ++i)
        {
            int bit = (idx >> (n - 1 - i)) & 1;
            int var = new_order_msb2lsb[i];
            int p = pos[var];
            old_idx |= (uint64_t(bit) << (n - 1 - p));
        }
        out[idx] = tt[old_idx];
    }
    return out;
}

// =====================================================
// 枚举 (x,y,z)
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
// 枚举变量划分
// =====================================================
inline std::vector<
    std::tuple<std::vector<int>,std::vector<int>,std::vector<int>>>
enumerate_variable_partitions(
    const std::vector<int>& vars,
    int x,int y,int z)
{
    std::vector<
        std::tuple<std::vector<int>,std::vector<int>,std::vector<int>>> out;

    int n = vars.size();
    std::vector<bool> sel_b(n,false);
    std::fill(sel_b.begin(), sel_b.begin()+y, true);

    do {
        std::vector<int> B, rest;
        for (int i=0;i<n;i++)
            (sel_b[i]?B:rest).push_back(vars[i]);

        if (x==0)
        {
            out.emplace_back(std::vector<int>{}, B, rest);
            continue;
        }

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
// index
// =====================================================
inline uint64_t mf_index(uint64_t a,uint64_t b,uint64_t c,int y,int z)
{
    return (a<<(y+z)) | (b<<z) | c;
}

// =====================================================
// Step 1: MF' → MXY
// =====================================================
inline std::string compute_MXYX_from_MF(
    const std::string& MF,int x,int y,int z)
{
    uint64_t H=pow2(x), M=pow2(y), L=pow2(z);
    std::string MXY(pow2(x+2*y+z),'x');

    for (uint64_t a=0;a<H;a++)
        for (uint64_t b=0;b<M;b++)
        {
            uint64_t blk=a*M*M+b*M+b;
            for (uint64_t c=0;c<L;c++)
                MXY[blk*L+c]=MF[mf_index(a,b,c,y,z)];
        }
    return MXY;
}

// =====================================================
// Step 2: solve MX / MY
// =====================================================
inline bool solve_MX_MY_from_MXY(
    const std::string& MXY,
    int x,int y,int z,
    std::string& MX,
    std::string& MY)
{
    BIDEC_DBG("[solve_MX_MY]");
    BIDEC_INDENT_INC();

    uint64_t MYN = pow2(x+y);
    uint64_t BS  = pow2(y+z);

    std::vector<std::string> blocks;
    MY.assign(MYN,'0');

    for (uint64_t i=0;i<MYN;i++)
    {
        std::string blk=MXY.substr(i*BS,BS);
        BIDEC_DBG("Block " << i << " = " << blk);

        bool hit=false;
        for (size_t k=0;k<blocks.size();k++)
        {
            bool ok=true;
            for (uint64_t j=0;j<BS;j++)
            {
                char a=blk[j], b=blocks[k][j];
                if (a!='x' && b!='x' && a!=b)
                {
                    ok=false;
                    BIDEC_DBG("  conflict at bit " << j);
                    break;
                }
            }
            if (ok)
            {
                for (uint64_t j=0;j<BS;j++)
                    if (blocks[k][j]=='x')
                        blocks[k][j]=blk[j];
                MY[i]=(k==0?'1':'0');
                hit=true;
                break;
            }
        }

        if (!hit)
        {
            if (blocks.size()==2)
            {
                BIDEC_DBG("  FAIL: >2 classes");
                BIDEC_INDENT_DEC();
                return false;
            }
            blocks.push_back(blk);
            MY[i]=(blocks.size()==1?'1':'0');
            BIDEC_DBG("  new class");
        }
    }

    if (blocks.size()!=2)
    {
        BIDEC_DBG("FAIL: class count !=2");
        BIDEC_INDENT_DEC();
        return false;
    }

    MX.clear();
    for (auto& b:blocks)
        for (char c:b)
            MX.push_back(c=='x'?'0':c);

    BIDEC_DBG("MY = " << MY);
    BIDEC_DBG("MX = " << MX);
    BIDEC_INDENT_DEC();
    return true;
}

// =====================================================
// ★ 主入口
// =====================================================
inline StrongBiDecResult
run_strong_bi_dec(const std::string& tt,
                  const std::vector<int>& order)
{
    StrongBiDecResult res;
    int n=order.size();
    if ((size_t(1)<<n)!=tt.size()) return res;

    BIDEC_DBG("[BiDec] start");
    BIDEC_DBG("order = " << vec2str(order));

    for (auto [x,y,z]:enumerate_xyz(n))
    {
        BIDEC_DBG("Try (x,y,z)=("<<x<<","<<y<<","<<z<<")");
        BIDEC_INDENT_INC();

        auto parts=enumerate_variable_partitions(order,x,y,z);
        for (auto& p:parts)
        {
            auto& A=std::get<0>(p);
            auto& B=std::get<1>(p);
            auto& C=std::get<2>(p);

            BIDEC_DBG("A="<<vec2str(A)<<" B="<<vec2str(B)<<" C="<<vec2str(C));

            std::vector<int> new_order=A;
            new_order.insert(new_order.end(),B.begin(),B.end());
            new_order.insert(new_order.end(),C.begin(),C.end());

            std::string MFp=reorder_tt_by_var_order(tt,n,new_order,order);
            std::string MXY=compute_MXYX_from_MF(MFp,x,y,z);

            std::string MX,MY;
            if (!solve_MX_MY_from_MXY(MXY,x,y,z,MX,MY))
                continue;

            res.found=true;
            res.A=A; res.B=B; res.C=C;
            res.MX=MX; res.MY=MY;

            BIDEC_DBG("[BiDec] SUCCESS");
            BIDEC_INDENT_DEC();
            return res;
        }
        BIDEC_INDENT_DEC();
    }
    return res;
}
