#include "map/if/if.h"

#include "misc/util/abc_global.h"
#include "misc/util/utilTruth.h"

#include <algorithm>
#include <vector>

// ===== ALSO / STP / mockturtle dependencies =====
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/topo_view.hpp>

#include "networks/aoig/stp_dsd_resynthesis.hpp"

namespace {

using AigNetwork = mockturtle::aig_network;
using AigSignal  = AigNetwork::signal;

// ----------------------------------------------------------------------
// Build kitty truth table from ABC cut
// ----------------------------------------------------------------------
bool If_AlsoBuildTruthTable( If_Man_t * p,
                             If_Cut_t * pCut,
                             kitty::dynamic_truth_table & tt )
{
    if ( !p || !pCut )
        return false;
    if ( !p->pPars->fTruth )
        return false;

    word * pTruth = If_CutTruthW( p, pCut );
    if ( !pTruth )
        return false;

    const int nLeaves = If_CutLeaveNum( pCut );
    tt = kitty::dynamic_truth_table( nLeaves );

    const int nBits = 1 << nLeaves;
    for ( int i = 0; i < nBits; ++i )
    {
        if ( Abc_TtGetBit( pTruth, i ) )
            kitty::set_bit( tt, i );
    }

    return true;
}

// ----------------------------------------------------------------------
// Convert mockturtle AIG into ABC mini-AIG representation
// ----------------------------------------------------------------------
int If_AlsoBuildMiniAig( AigNetwork const & ntk,
                         AigSignal const & root,
                         Vec_Int_t * vAig,
                         int const * pTimes,
                         int nLeaves,
                         int * pArea )
{
    if ( vAig )
        Vec_IntClear( vAig );

    const auto root_node = ntk.get_node( root );

    // constant
    if ( ntk.is_constant( root_node ) )
    {
        if ( vAig )
            Vec_IntPush( vAig, ntk.is_complemented( root ) ? 1 : 0 );
        if ( pArea )
            *pArea = 0;
        return 0;
    }

    // PI as root -> unsupported
    if ( ntk.is_pi( root_node ) )
        return -1;

    std::vector<int> var_map( ntk.size(), -1 );
    std::vector<int> delay_map( ntk.size(), 0 );

    // map leaves
    int leaf_index = 0;
    ntk.foreach_pi( [&]( auto n ) {
        if ( leaf_index < nLeaves )
        {
            var_map[n]   = leaf_index;
            delay_map[n] = pTimes[leaf_index];
        }
        ++leaf_index;
    } );

    int gate_index = 0;
    bool ok = true;

    mockturtle::topo_view<AigNetwork> topo( ntk );
    topo.foreach_gate( [&]( auto n ) {
        AigSignal fanins[2];
        int fanin_idx = 0;

        ntk.foreach_fanin( n, [&]( AigSignal const & f ) {
            fanins[fanin_idx++] = f;
        } );

        const int var0 = var_map[ntk.get_node( fanins[0] )];
        const int var1 = var_map[ntk.get_node( fanins[1] )];

        if ( var0 < 0 || var1 < 0 )
        {
            ok = false;
            return;
        }

        const int lit0 = Abc_Var2Lit( var0, ntk.is_complemented( fanins[0] ) );
        const int lit1 = Abc_Var2Lit( var1, ntk.is_complemented( fanins[1] ) );

        if ( vAig )
        {
            Vec_IntPush( vAig, lit0 );
            Vec_IntPush( vAig, lit1 );
        }

        const int var = nLeaves + gate_index++;
        var_map[n]   = var;
        delay_map[n] = 1 + std::max(
            delay_map[ntk.get_node( fanins[0] )],
            delay_map[ntk.get_node( fanins[1] )]
        );
    } );

    if ( !ok )
        return -1;

    const int root_var = var_map[root_node];
    if ( root_var < 0 )
        return -1;

    if ( vAig )
        Vec_IntPush( vAig,
                     Abc_Var2Lit( root_var,
                                  ntk.is_complemented( root ) ) );

    if ( pArea )
        *pArea = gate_index;

    return delay_map[root_node];
}

} // namespace

// ======================================================================
// Public entry called from ifDsdAlso.cpp
// ======================================================================
int If_CutDsdBalanceEvalAlso_Impl( If_Man_t * p,
                                  If_Cut_t * pCut,
                                  Vec_Int_t * vAig )
{
    if ( !p || !pCut || pCut->nLeaves <= 1 )
        return -1;

    kitty::dynamic_truth_table tt;
    if ( !If_AlsoBuildTruthTable( p, pCut, tt ) )
        return -1;

    const int nLeaves = If_CutLeaveNum( pCut );

    std::vector<AigSignal> leaves;
    leaves.reserve( nLeaves );

    AigNetwork ntk;
    for ( int i = 0; i < nLeaves; ++i )
        leaves.push_back( ntk.create_pi() );

    AigSignal root;
    bool has_root = false;

    also::stp_dsd_lut_resynthesis<AigNetwork> resyn;
    resyn( ntk, tt, leaves.begin(), leaves.end(),
           [&]( AigSignal const & signal ) {
               root     = signal;
               has_root = true;
           } );

    if ( !has_root )
        return -1;

    int pTimes[IF_MAX_FUNC_LUTSIZE];
    for ( int i = 0; i < nLeaves; ++i )
        pTimes[i] =
            (int)If_ObjCutBest( If_CutLeaf( p, pCut, i ) )->Delay;

    int area  = 0;
    int delay = If_AlsoBuildMiniAig(
        ntk, root, vAig, pTimes, nLeaves, &area
    );

    if ( delay >= 0 )
        pCut->Cost = area;

    return delay;
}
