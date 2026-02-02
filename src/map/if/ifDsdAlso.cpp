#include "map/if/if.h"
#include "misc/util/abc_global.h"

#ifdef ABC_USE_ALSO
// 不要 include 任何 kitty/mockturtle/percy/stp 的头
// 只声明实现函数（定义在另一个 .cpp 里）
int If_CutDsdBalanceEvalAlso_Impl( If_Man_t * p, If_Cut_t * pCut, Vec_Int_t * vAig );
#endif

int If_CutDsdBalanceEvalAlso( If_Man_t * p, If_Cut_t * pCut, Vec_Int_t * vAig )
{
#ifndef ABC_USE_ALSO
    (void)p;
    (void)pCut;
    (void)vAig;
    return -1;
#else
    return If_CutDsdBalanceEvalAlso_Impl( p, pCut, vAig );
#endif
}
