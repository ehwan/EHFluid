#include "engine.h"
#include "fluidgroup.h"

#include <iostream>


FluidEngine::FluidEngine()
{

}
void FluidEngine::Load( const vector2& minbound , const vector2& maxbound , unsigned int maxc )
{
    group.reserve( GROUP_SIZE );

    pairlist.reserve(
            ( maxc*(maxc-1) ) >> 3
            );

    aabb.min = minbound - METER_PER_GRID*1.1f;
    aabb.max = maxbound + METER_PER_GRID*1.1f;
    aabb.mingrid = aabb.min * GRID_PER_METER;

    //const vector2 gsize = (aabb.max - aabb.min) * GRID_PER_METER;
    const vector2 gsize = EH::Matrix::fma( aabb.max , GRID_PER_METER , -aabb.mingrid );
    grid_attrib.Load( (unsigned int)std::ceil( gsize[0] ) , (unsigned int)std::ceil( gsize[1] ) , maxc );
}

void FluidEngine::AddParticle( const vector2& pos , const vector2& vel , FluidGroup *grp )
{
    const vector2 gridf = EH::Matrix::fma( pos , GRID_PER_METER , -aabb.mingrid );
    //const vector2 gridf = ( pos - aabb.min )*FluidEngine::GRID_PER_METER;
    const unsigned int gx = (unsigned int)std::floor( gridf[0] );
    const unsigned int gy = (unsigned int)std::floor( gridf[1] );

    if( gx < grid_attrib.w && gy < grid_attrib.h )
    {
        grp->particle.Add( pos , vel , GetGridIndex( gx , gy ) );
    }
}
void FluidEngine::AddGroup( FluidGroup *grp )
{
    const auto iter = std::upper_bound( group.begin() , group.end() , grp ,
            []( const FluidGroup *g1 , const FluidGroup *g2 )->bool
            { return g1->type < g2->type; }
            );
    group.insert( iter , grp );
    //group.push_back( grp );
    //std::sort( group.begin() , group.end() , []( const FluidGroup *g1 , const FluidGroup *g2 )->bool { return g1->type<g2->type; } );
}



void FluidEngine::GridSort()
{
using std::fill;
using std::partial_sum;
    fill( grid_attrib.gridcount , grid_attrib.gridcount + grid_attrib.count , 0 );

    for( FluidGroup* grp : group )
    {
        if( grp->type != FluidGroup::STATIC )
        {
            unsigned int i = grp->particle.count;
            //while( i<c )
            while( i )
            {
                --i;
                const vector2 gridf = EH::Matrix::fma( grp->particle.position[i] , GRID_PER_METER , -aabb.mingrid );
                const unsigned int gx = (unsigned int)std::floor(gridf[0]);
                const unsigned int gy = (unsigned int)std::floor(gridf[1]);

                /*
                 *      out of bound check  i>=min && i<=max
                 *      try this
                 *      (unsigned)( i-min ) <= (unsigned)( max-min )
                 */
                if( gx < grid_attrib.w && gy < grid_attrib.h )
                {
                    grp->particle.grid[i] = GetGridIndex( gx , gy );
                }
                else if( grp->autoremove )
                {
                    grp->particle.DeleteBack( i );
                }else
                {
                    grp->particle.grid[i] = grid_attrib.count;
                }
            }
        }
    }
    for( FluidGroup* grp : group )
    {
using std::transform;
        unsigned int *grdcnt = grid_attrib.gridcount;
        transform( grp->particle.grid , grp->particle.grid + grp->particle.count , grp->particle.gridindex ,
                [grdcnt]( const unsigned int grd_id ) -> unsigned int
                {
                    return grdcnt[ grd_id ]++;
                }
                );

    }
    partial_sum( grid_attrib.gridcount , grid_attrib.gridcount + grid_attrib.count , grid_attrib.gridoffset + 1 );
    for( FluidGroup *grp : group )
    {
        const unsigned int *grd = grp->particle.grid;
        const unsigned int *grdind = grp->particle.gridindex;
        const unsigned int c = grp->particle.count;
        for( unsigned int i=0; i<c; ++i )
        {
            grid_attrib.grid[
                grid_attrib.gridoffset[ grd[i] ] + grdind[i]
                ].Set( grp , i );
        }
    }
}
void FluidEngine::Pressure()
{
    pairlist.clear();

    for( FluidGroup *grp : group )
    {
        grp->particle.Init();
    }

    const unsigned int left_top_gap = GetGridIndex( 1 , 1 );
    unsigned int i,j;

    //std::sort( onlinegrid.begin() , onlinegrid.end() );
    const auto& grdoffs = grid_attrib.gridoffset;

    for( unsigned int gy=1; gy<grid_attrib.h-1; ++gy )
    {
        for( unsigned int gx=1; gx<grid_attrib.w-1; ++gx)
        {
            const unsigned int g1 = GetGridIndex( gx , gy );
            const unsigned int s1 = grdoffs[ g1 - left_top_gap ];       // left-top
            const unsigned int s2 = grdoffs[ g1 - left_top_gap + 3 ];   // right-top
            const unsigned int s3 = grdoffs[ g1 - 1 ];                  // left
                                                        //s4            // cur
            const unsigned int s5 = grdoffs[ g1 ];
            const unsigned int s6 = grdoffs[ g1+1 ];
            for( i=s5; i<s6; ++i )
            {
                const particle_index_type& p1 = grid_attrib.grid[ i ];

                for( j=s1; j<s2; ++j )
                {
                    AddPair( p1 , grid_attrib.grid[j] );
                }
                for( j=s3; j<i; ++j )
                {
                    AddPair( p1 , grid_attrib.grid[j] );
                }
            }
        }
    }

    for( FluidGroup *grp : group )
    {
        grp->CalculatePressure();
    }
}
void FluidEngine::AddPair( const particle_index_type& p1 , const particle_index_type& p2 )
{
    const bool same_group = p1.group == p2.group;
    if( same_group && p1.group->pressureconnected==false ){ return; }
    vector2 rel = p2.group->particle.position[p2.index] - p1.group->particle.position[p1.index];
    const float lsq = rel.LengthSquared();
    if( lsq < H2 )
    {
        const float h = 1.0f - lsq*INV_H2;
        const float m1 = p1.group->mass.m;
        const float m2 = p2.group->mass.m;

        p1.group->particle.rho[p1.index] += m2*h;
        p2.group->particle.rho[p2.index] += m1*h;

        if( same_group && p1.group->collideconnected==false ){ return; }
        if( p1.group->type==FluidGroup::STATIC && p2.group->type==FluidGroup::STATIC ){ return; }

        //pairlist.resize( pairlist.size() + 1 );
        pairlist.emplace_back();
        ParticlePair& pair = pairlist.back();
        pair.length = rel.Normalize();
        pair.q = 1.0f - pair.length*INV_H;
        pair.q2 = pair.q * pair.q;
        pair.relq = rel * pair.q2;
        pair.force1 = p1.group->particle.force + p1.index;
        pair.force2 = p2.group->particle.force + p2.index;
        pair.press1 = p1.group->particle.pressure + p1.index;
        pair.press2 = p2.group->particle.pressure + p2.index;
        pair.veldif1 = p1.group->particle.veldiffpressure + p1.index;
        pair.veldif2 = p2.group->particle.veldiffpressure + p2.index;
        pair.omega1 = p1.group->particle.angularpressure + p1.index;
        pair.omega2 = p2.group->particle.angularpressure + p2.index;
        pair.vortice1 = p1.group->particle.vorticity + p1.index;
        pair.vortice2 = p2.group->particle.vorticity + p2.index;

        const vector2 dv = ( p2.group->particle.velocity[p2.index] - p1.group->particle.velocity[p1.index] );
        const vector2 dvf = dv * pair.q2 * p1.group->viscosity * p2.group->viscosity;
        *pair.force1 += dvf;
        *pair.force2 -= dvf;

        p1.group->particle.surfacenormal[ p1.index ] -= pair.relq;
        p2.group->particle.surfacenormal[ p2.index ] += pair.relq;

        //const float dvr = -dv.dot( rel ) * h;
        const float dvr = -EH::Matrix::ATxB( dv , rel ) * h;
        *pair.veldif1 += dvr*m2;
        *pair.veldif2 += dvr*m1;

        //const float tempw =  rel.cross( dv ) * h;
        const float tempw = EH::Matrix::Cross( rel , dv ) * h;
        *pair.omega1 += tempw;
        *pair.omega2 += tempw;
    }
}
void FluidEngine::Force()
{
using std::for_each;

    for_each( pairlist.cbegin() , pairlist.cend() ,
            []( const ParticlePair&  pair )
            {
                const vector2 f =
                    (
                     ( *pair.press1 + *pair.press2 ) * pair.q
                     + ( *pair.veldif1 + *pair.veldif2 )
                    ) * pair.relq;

                const vector2 vort =
                    ( std::fabs(*pair.omega1) - std::fabs(*pair.omega2) ) * pair.relq;
                *pair.vortice1 -= vort;
                *pair.vortice2 += vort;

                *pair.force1 -= f;
                *pair.force2 += f;

            }
            );
}
void FluidEngine::Advect( float dt )
{
    for( FluidGroup *grp : group )
    {
        grp->Advect( dt );
    }
}

void FluidEngine::Step( float dt )
{
    GridSort();
    Pressure();
    Force();
    Advect( dt );
}

#include <random>
#include <chrono>
int main()
{
    FluidEngine engine;
    engine.Load(
            vector2{{ -3.0f , -3.0f}},
            vector2{{ 3.0f , 3.0f }},
            1000u
            );
    std::default_random_engine randeng;
    std::uniform_real_distribution<float> randdist( -3.0f , 3.0f );

    FluidGroup grp1;
    grp1.particle.Load( 1000 );
    engine.AddGroup( &grp1 );

    for( auto i=0; i<1000; ++i )
    {
        engine.AddParticle(
                { randdist(randeng) , randdist(randeng) },
                { 0.0f , 0.0f },
                &grp1
                );
    }
    std::chrono::time_point< std::chrono::high_resolution_clock > start = std::chrono::high_resolution_clock::now();
    engine.Step( 1.0f/60.0f );
    std::chrono::time_point< std::chrono::high_resolution_clock > end = std::chrono::high_resolution_clock::now();

    std::chrono::duration< float , std::ratio<1,1>  > gap = end-start;
    std::cout << gap.count() << " Seconds elapsed; " << std::endl;


    std::cout << "Pair Count : " << engine.pairlist.size() << " / " << engine.pairlist.capacity() << std::endl;
}
