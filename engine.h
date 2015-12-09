#pragma once

#include <utility>
#include <numeric>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
//#include <valarray>
#include "../Matrix/EHMatrix.h"
//#include "vector2.h"

#ifndef VECTOR2_TYPE
#define VECTOR2_TYPE

typedef EH::Matrix::Matrix< float , 2 , 1 > vector2;
template < typename T >
using vec2 = EH::Matrix::Matrix< T , 2 , 1 >;
//#include "vector2.h"

#endif

struct FluidGroup;

class FluidEngine
{
public:
    constexpr static float H = 0.3f;
    constexpr static float H2 = H*H;
    constexpr static float INV_H = 1.0f/H;
    constexpr static float INV_H2 = 1.0f/(H*H);

    constexpr static float METER_PER_GRID = H;
    constexpr static float GRID_PER_METER = 1.0f/METER_PER_GRID;

    constexpr static std::size_t GROUP_SIZE = 10;

    struct particle_index_type
    {
        FluidGroup *group;
        unsigned int index;
        particle_index_type(){}
        void Set( FluidGroup *grp , unsigned int ind )
        {
            group = grp;
            index = ind;
        }
    };
    // p2's relative to p1
    struct ParticlePair
    {
        vector2 relq;
        float length;
        float q,q2;
        vector2 *force1 , *force2;
        float *press1 , *press2;
        float *veldif1 , *veldif2;
        float *omega1 , *omega2;
        vector2 *vortice1 , *vortice2;

        ParticlePair(){}
    };

    struct GridContainer
    {
        unsigned int w , h;
        unsigned int count;
        particle_index_type *grid;
        unsigned int *gridcount;
        unsigned int *gridoffset;

        GridContainer()
        {
            std::cout << "Grid Construct" << "\n";
            grid = 0;
            gridcount = 0;
            gridoffset = 0;
        }
        ~GridContainer()
        {
            std::cout << "Grid Destruct" << "\n";
            if( grid ){ delete[] grid; }
            if( gridcount ){ delete[] gridcount; }
            if( gridoffset ){ delete[] gridoffset; }
        }
        void Load( unsigned int _w , unsigned int _h , unsigned int maxc )
        {
            std::cout << "Grid Load" << "\n";
            w = _w;
            h = _h;
            count = _w * _h;
            grid = new particle_index_type[ maxc ];
            gridcount = new unsigned int[ count ];
            gridoffset = new unsigned int[ count + 1 ];
            gridoffset[0] = 0;
        }
    } grid_attrib;

    std::vector< FluidGroup* > group;
    std::vector< ParticlePair > pairlist;

    struct
    {
        vector2 min;
        vector2 max;
        vector2 mingrid;
    } aabb;

    //int gridsize_bitfactor;
    inline unsigned int GetGridIndex( unsigned int x , unsigned int y )
    {
        return y*grid_attrib.w + x;
        //return x + (y<<gridsize_bitfactor);
    }

    FluidEngine();

    void Load( const vector2& minbound , const vector2& maxbound , unsigned int maxc );

    void AddParticle( const vector2& pos , const vector2& vel , FluidGroup *grp );
    void AddGroup( FluidGroup * grp );
    inline std::vector< FluidGroup* >& GetGroups() {return group;}

    void GridSort();
    void AddPair( const particle_index_type& p1 , const particle_index_type& p2 );
    void Pressure();
    void Force();
    void Advect( float dt );
    void Step( float dt );

protected:
};
