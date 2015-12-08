#pragma once

#include <algorithm>
#include <numeric>
#include <vector>
#include "../Matrix/EHMatrix.h"

#ifndef VECTOR2_TYPE
#define VECTOR2_TYPE

typedef EH::EHMatrix< float , 2 , 1 > vector2;
template < typename T >
using vec2 = EH::EHMatrix< T , 2 , 1 >;
//#include "vector2.h"
//
#endif

struct FluidGroup
{
    enum GroupType
    {
        STATIC = 0,
        FLUID = 1,
        BODY = 2,
        ROPE = 3
    } type;
    bool autoremove;
    bool pressureconnected;
    bool collideconnected;

    struct _Container
    {
        unsigned int count;
        vector2 *position;
        vector2 *velocity;
        vector2 *force;
        vector2 *relative;
        vector2 *surfacenormal;
        unsigned int *grid;
        float *rho;
        float *pressure;
        float *veldiffpressure;
        float *angularpressure;
        vector2 *vorticity;
        unsigned int *gridindex;

        _Container()
        {
            position = 0;
            velocity = 0;
            force = 0;
            relative = 0;
            surfacenormal = 0;
            grid = 0;
            rho = 0;
            pressure = 0;
            veldiffpressure = 0;
            angularpressure = 0;
            vorticity = 0;
            gridindex = 0;
        }
        ~_Container()
        {
            if( position )  { delete[] position; }
            if( velocity )  { delete[] velocity; }
            if( force )     { delete[] force; }
            if( relative )  { delete[] relative; }
            if( surfacenormal ) { delete[] surfacenormal; }
            if( grid )      { delete[] grid; }
            if( rho )       { delete[] rho; }
            if( pressure )  { delete[] pressure; }
            if( veldiffpressure ) { delete[] veldiffpressure; }
            if( angularpressure ) { delete[] angularpressure; }
            if( vorticity ) { delete[] vorticity; }
            if( gridindex ) { delete[] gridindex; }
        }
        void Init()
        {
            using std::fill;
            fill( rho , rho + count , 0.0f );
            fill( veldiffpressure , veldiffpressure + count , 0.0f );
            fill( angularpressure , angularpressure + count , 0.0f );
            fill( surfacenormal , surfacenormal + count , vector2(0.0f) );
            fill( vorticity , vorticity + count , vector2(0.0f) );
            //fill( pressure.begin() , pressure.end() , 0.0f );
        }
        void Load( std::size_t size )
        {
            count = 0;

            position                        = new vector2[ size ];
            velocity                        = new vector2[ size ];
            force                           = new vector2[ size ];
            relative                        = new vector2[ size ];
            surfacenormal                   = new vector2[ size ];
            grid                            = new unsigned int[ size ];
            rho                             = new float[ size ];
            pressure                        = new float[ size ];
            veldiffpressure                 = new float[ size ];
            angularpressure                 = new float[ size ];
            vorticity                       = new vector2[ size ];
            gridindex                       = new unsigned int[ size ];
        }
        void DeleteBack( unsigned int pos )
        {
            --count;
            if( pos == count ){ return; }
            position[ pos ] = position[ count ];
            velocity[ pos ] = velocity[ count ];
            grid[ pos ] = grid[ count ];
            force[ pos ] = force[ count ];
        }
        void Add( const vector2& pos , const vector2& vel , unsigned int grd )
        {
            position[ count ] = pos;
            velocity[ count ] = vel;
            grid[ count ] = grd;
            force[ count ] = vector2(0.0f);

            ++count;
        }
    } particle;

    vector2 gravity;

    struct
    {
        float m;
        float invm;
        float ascale;

        void Set( float mas , float invmas )
        {
            m = mas;
            invm = invmas;
            ascale = mas * invmas;
        }
        void Set( float mas )
        {
            m = mas;
            invm = 1.0f/mas;
            ascale = 1.0f;
        }
        void operator =( float mas )
        {
            Set( mas );
        }

    } mass;

    float rest_density;
    float pressurek;
    float veldiffk;
    float vorticityk;

    float viscosity;

    FluidGroup()
    {
        pressureconnected = true;
        collideconnected = true;
        autoremove = true;

        type = FLUID;
        gravity = { 0.0f , 0.0f };

        mass = 1.0f;

        rest_density = 3.0f;
        pressurek = 1.2f;
        veldiffk = 10.0f;
        vorticityk = 3.0f;

        viscosity = 0.3f;
    }
    void CalculatePressure()
    {
using std::transform;
        const float veldifk = veldiffk;
        transform( particle.veldiffpressure , particle.veldiffpressure + particle.count , particle.rho , particle.veldiffpressure ,
                [ veldifk ]( float dv , float p )->float
                {
                    return dv/p * veldifk;
                }
                );
        const float restdens = rest_density;
        const float pk = pressurek;
        const float mas = mass.m;
        transform( particle.rho , particle.rho + particle.count , particle.pressure ,
                [ mas , restdens , pk ]( float p )->float
                {
                    return ( p + mas - restdens ) * pk;
                }
                );
        const float vortk = vorticityk;
        transform( particle.angularpressure , particle.angularpressure + particle.count , particle.angularpressure ,
                [ vortk ]( float p )->float
                {
                    return p * vortk;
                }
                );
    }

    void ClampVelocity()
    {
constexpr static float MAX_VELOCITY_MAG = 10.0f;
using std::for_each;
        for_each( particle.velocity , particle.velocity + particle.count ,
                []( vector2& vel )
                {
                    if( vel.LengthSquared() > MAX_VELOCITY_MAG*MAX_VELOCITY_MAG )
                    {
                        vel.Normalize();
                        vel *= MAX_VELOCITY_MAG;
                    }
                    ////or non-branching code.
                    //const float l = std::fmin( vel.Normalize() , MAX_VELOCITY_MAG );
                    //vel *= l;
                }
                );
    }
    virtual void Advect( float dt )
    {
using std::transform;
using std::fill;
        const vector2 gravdt = gravity * (mass.ascale*dt);
        const float invmdt = mass.invm*dt;
        for( decltype( particle.count ) i=0; i<particle.count; ++i )
        {
            particle.vorticity[i].Normalize();
            particle.force[i] =
                gravdt + invmdt * ( particle.force[i]
                        + EH::Cross( particle.vorticity[i] , particle.angularpressure[i] )
                        );
        }
        transform( particle.velocity , particle.velocity + particle.count , particle.force , particle.velocity ,
                []( const vector2& vel , const vector2& force )->vector2
                {
                    return vel + force;
                }
                );
        for( decltype( particle.count ) i=0; i<particle.count; ++i )
        {
            particle.position[i] += dt * ( particle.velocity[i] + 0.5f*particle.force[i] );
        }
        fill( particle.force , particle.force + particle.count , vector2(0.0f) );

        ClampVelocity();
    }
};

struct FluidBody : public FluidGroup
{
    vector2 position;
    vector2 velocity;
    float omega;

    float invI;
    float invN;

    FluidBody():FluidGroup()
    {
        type = FluidGroup::BODY;
    }

    void LoadBody( float invIfactor=1.0f )
    {
using std::transform;
using std::accumulate;
        invN = 1.0f/(float)particle.count;

        const vector2 centerpos = accumulate( particle.position , particle.position + particle.count ,
                vector2( 0.0f ) ) * invN;

        position = centerpos;
        transform( particle.position , particle.position + particle.count , particle.relative ,
                [&centerpos]( const vector2& partpos )->vector2
                {
                    return partpos - centerpos;
                }
                );
        const float inertia = accumulate( particle.relative , particle.relative + particle.count , 0.0f ,
                []( const float cur , const vector2& rel )->float
                {
                    return cur + rel.LengthSquared();
                }
                ) * mass.m;
        invI = 1.0f/inertia * invIfactor;
    }

    void Advect( float dt ) override
    {
using std::accumulate;
using std::inner_product;
using std::transform;
using std::for_each;
        const vector2 adt = ( accumulate( particle.force , particle.force + particle.count , vector2(0.0f) )*mass.invm*invN + gravity*mass.ascale ) * dt;
        velocity += adt;
        position += ( velocity + 0.5f*adt ) * dt;

        const float alphadt = inner_product( particle.force , particle.force + particle.count , particle.relative , 0.0f ,
                std::plus<float>() ,
                []( const vector2& force , const vector2& rel )->float
                {
                    return EH::Cross( rel , force );
                }
                )
            * invI * dt;
        omega += alphadt;

        const float dtheta = ( omega + 0.5f*alphadt ) * dt;
        const vector2 dthetap = EH::Complex::Complex( dtheta );


        for_each( particle.relative , particle.relative ,
                [&dthetap]( vector2& rel )
                {
                    rel = EH::Complex::Multiply( rel , dthetap );
                    //rel.rotate( dthetap.x , dthetap.y );
                }
                );

        const vector2 pos = position;
        transform( particle.relative , particle.relative + particle.count , particle.position ,
                [&pos]( const vector2& rel )->vector2
                {
                    return pos + rel;
                }
                );
        const vector2 vel = velocity;
        const float omg = omega;
        transform( particle.relative , particle.relative + particle.count , particle.velocity ,
                [&vel , omg]( const vector2& rel )->vector2
                {
                    return EH::Cross( omg , rel );
                }
                );
    }
};

struct FluidRope : public FluidGroup
{
constexpr static float FIX_FACTOR = 0.9f;
constexpr static float VELOCITY_FACTOR = 0.75f;
constexpr static unsigned int ITERATION = 5;

    struct RopePair
    {
        vector2 *pos1 , *pos2;
        vector2 *vel1 , *vel2;
    };
    std::vector< RopePair > ropepair;
    float rope_length;


    void LoadRope()
    {
        ropepair.resize( particle.count - 1 );
        for( unsigned int i=1; i<particle.count; ++i )
        {
            ropepair[i-1] =
            {
                &particle.position[i-1] , &particle.position[i] ,
                &particle.velocity[i-1] , &particle.velocity[i]
            };
        }
    }
    void Iterate()
    {
        const float rl = rope_length;
        std::for_each( ropepair.cbegin() , ropepair.cend() ,
                [rl]( const RopePair& pair )
                {
                    vector2 rel = *pair.pos2 - *pair.pos1;
                    const float gap = rel.Normalize() - rl;
                    rel *= gap * FIX_FACTOR * 0.5f;
                    *pair.pos1 += rel;
                    *pair.pos2 -= rel;
                }
                );
    }
    void Advect( float dt ) override
    {
using std::copy;
using std::transform;
        FluidGroup::Advect( dt );
        copy( particle.position , particle.position + particle.count , particle.relative );
        for( unsigned int i=0; i<ITERATION; ++i )
        {
            Iterate();
        }
        const float invdtfactor = VELOCITY_FACTOR/dt;
        for( decltype( particle.count ) i=0; i<particle.count; ++i )
        {
            particle.velocity[i] += ( particle.position[i] - particle.relative[i] ) * ( invdtfactor ) ;
        }
    }
};
