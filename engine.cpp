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

    const vec2< unsigned int > gsize = EH::Matrix::ceil( (aabb.max - aabb.min)*GRID_PER_METER ).Convert< unsigned int >();
    grid_attrib.Load( gsize.x , gsize.y , maxc );
}

void FluidEngine::AddParticle( const vector2& pos , const vector2& vel , FluidGroup *grp )
{
    const vec2< unsigned int > gr = EH::Matrix::floor( ( pos - aabb.min ) * GRID_PER_METER ).Convert< unsigned int >();

    if( gr.x < grid_attrib.w && gr.y < grid_attrib.h )
    {
        grp->particle.Add( pos , vel , GetGridIndex( gr.x , gr.y ) );
    }
}
void FluidEngine::AddGroup( FluidGroup *grp )
{
    const auto iter = std::upper_bound( group.begin() , group.end() , grp ,
            []( const FluidGroup *g1 , const FluidGroup *g2 )->bool
            { return g1->type < g2->type; }
            );
    group.insert( iter , grp );
}



void FluidEngine::GridSort()
{
    for( FluidGroup* grp : group )
    {
        if( grp->type != FluidGroup::STATIC )
        {
            unsigned int i = grp->particle.count;
            //while( i<c )
            while( i )
            {
                --i;
                const vec2< unsigned int > gr = EH::Matrix::floor( (grp->particle.position[i]-aabb.min) * GRID_PER_METER ).Convert< unsigned int >();

                /*
                 *      out of bound check  i>=min && i<=max
                 *      try this
                 *      (unsigned)( i-min ) <= (unsigned)( max-min )
                 */
                if( gr.x < grid_attrib.w && gr.y < grid_attrib.h )
                {
                    grp->particle.grid[i] = GetGridIndex( gr.x , gr.y );
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
    std::fill( grid_attrib.gridcount.get() , grid_attrib.gridcount.get() + grid_attrib.count+1 , 0 );
    for( FluidGroup* grp : group )
    {
        for( unsigned int i=0; i<grp->particle.count; ++i )
        {
            grp->particle.gridindex[i] = grid_attrib.gridcount[ grp->particle.grid[i] ]++;
        }
    }
    std::partial_sum( grid_attrib.gridcount.get() , grid_attrib.gridcount.get() + grid_attrib.count , grid_attrib.gridoffset.get() + 1 );
    for( FluidGroup *grp : group )
    {
        for( unsigned int i=0; i<grp->particle.count; ++i )
        {
            grid_attrib.grid[
                grid_attrib.gridoffset[ grp->particle.grid[i] ] + grp->particle.gridindex[i]
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

    for( unsigned int gy=1; gy<grid_attrib.h-1; ++gy )
    {
        for( unsigned int gx=1; gx<grid_attrib.w-1; ++gx)
        {
            const unsigned int g1 = GetGridIndex( gx , gy );
            const unsigned int s1 = grid_attrib.gridoffset[ g1 - left_top_gap ];       // left-top
            const unsigned int s2 = grid_attrib.gridoffset[ g1 - left_top_gap + 3 ];   // right-top
            const unsigned int s3 = grid_attrib.gridoffset[ g1 - 1 ];                  // left
                                                        //s4            // cur
            const unsigned int s5 = grid_attrib.gridoffset[ g1 ];
            const unsigned int s6 = grid_attrib.gridoffset[ g1+1 ];
            for( i=s5; i<s6; ++i )
            {
                const particle_index_type p1 = grid_attrib.grid[ i ];

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
        pair.force1 = p1.group->particle.force.get() + p1.index;
        pair.force2 = p2.group->particle.force.get() + p2.index;
        pair.press1 = p1.group->particle.pressure.get() + p1.index;
        pair.press2 = p2.group->particle.pressure.get() + p2.index;
        pair.veldif1 = p1.group->particle.veldiffpressure.get() + p1.index;
        pair.veldif2 = p2.group->particle.veldiffpressure.get() + p2.index;

        const vector2 dv = ( p2.group->particle.velocity[p2.index] - p1.group->particle.velocity[p1.index] );
        const vector2 dvf = dv * pair.q2 * p1.group->viscosity * p2.group->viscosity;
        *pair.force1 += dvf;
        *pair.force2 -= dvf;

        p1.group->particle.surfacenormal[ p1.index ] -= pair.relq;
        p2.group->particle.surfacenormal[ p2.index ] += pair.relq;

        //const float dvr = -dv.dot( rel ) * h;
        const float dvr = -h * ( dv.Transpose()*rel );
        *pair.veldif1 += dvr*m2;
        *pair.veldif2 += dvr*m1;
    }
}
void FluidEngine::Force()
{
    for( ParticlePair& pair : pairlist )
    {
        const vector2 f =
            (
             ( *pair.press1 + *pair.press2 ) * pair.q
             + ( *pair.veldif1 + *pair.veldif2 )
            ) * pair.relq;

        *pair.force1 -= f;
        *pair.force2 += f;

    }
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
/*
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
                0.0f ,
                &grp1
                );
    }
    std::chrono::time_point< std::chrono::high_resolution_clock > start = std::chrono::high_resolution_clock::now();
    engine.Step( 1.0f/60.0f );
    std::chrono::time_point< std::chrono::high_resolution_clock > end = std::chrono::high_resolution_clock::now();

    std::chrono::duration< float , std::ratio<1,1>  > gap = end-start;
    std::cout << gap.count() << " Seconds elapsed; " << std::endl;


    std::cout << "Pair Count : " << engine.pairlist.size() << " / " << engine.pairlist.capacity() << std::endl;
}*/
