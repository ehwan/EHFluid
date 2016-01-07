#pragma once

#include <utility>
#include <numeric>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <memory>

#include "fluid_global.h"


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
        const float *press1 , *press2;
        const float *veldif1 , *veldif2;
    };

    struct GridContainer
    {
        unsigned int w , h;
        unsigned int count;
        std::unique_ptr< particle_index_type[] > grid;
        std::unique_ptr< unsigned int[] > gridcount , gridoffset;

        void Load( unsigned int _w , unsigned int _h , unsigned int maxc )
        {
            std::cout << "Grid Load" << "\n";
            w = _w;
            h = _h;
            count = _w * _h;
            grid = std::template make_unique< particle_index_type[] >( maxc );
            gridcount = std::template make_unique< unsigned int[] >( count+1 );
            gridoffset = std::template make_unique< unsigned int[] >( count + 1 );
            gridoffset[0] = 0;
        }
    } grid_attrib;

    std::vector< FluidGroup* > group;
    std::vector< ParticlePair > pairlist;

    struct
    {
        vector2 min;
        vector2 max;
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
