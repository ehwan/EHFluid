#pragma once

#include <memory>
#include <algorithm>
#include <vector>
#include "fluid_global.h"

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
        unsigned int max;
        std::unique_ptr< vector2[] > position ,
                                   velocity ,
                                   force ,
                                   relative ,
                                   surfacenormal;
        std::unique_ptr< float[] >   rho ,
                                   pressure ,
                                   veldiffpressure;
        std::unique_ptr< unsigned int[] > grid ,
                                        gridindex;
        void Init()
        {
            std::fill( rho.get() , rho.get() + count , 0.0f );
            std::fill( veldiffpressure.get() , veldiffpressure.get() + count , 0.0f );
            std::fill( surfacenormal.get() , surfacenormal.get() + count , vector2(0.0f) );
            //fill( pressure.begin() , pressure.end() , 0.0f );
        }
        void Load( std::size_t size )
        {
            using std::make_unique;
            max = size;
            count = 0;

            position            = make_unique< vector2[] >( size );
            velocity            = make_unique< vector2[] >( size );
            force               = make_unique< vector2[] >( size );
            relative            = make_unique< vector2[] >( size );
            surfacenormal       = make_unique< vector2[] >( size );
            grid                = make_unique< unsigned int[] >( size );
            rho                 = make_unique< float[] >( size );
            pressure            = make_unique< float[] >( size );
            veldiffpressure     = make_unique< float[] >( size );
            gridindex           = make_unique< unsigned int[] >( size );
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
            if( count == max ){ return; }
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

        operator float () const
        {
            return m;
        }

    } mass;

    float rest_density;
    float pressurek;
    float veldiffk;

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

        viscosity = 0.3f;
    }
    void CalculatePressure()
    {
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.veldiffpressure[i] *= veldiffk/particle.rho[i];
        }
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.pressure[i] = ( particle.rho[i] + mass - rest_density ) * pressurek;
        }
        //for( unsigned int i=0; i<particle.count; ++i )
        //{
            //particle.angularpressure[i] *= vorticityk/particle.rho[i];
        //}
    }

    void ClampVelocity()
    {
constexpr static float MAX_VELOCITY_MAG = 10.0f;
        for( unsigned int i=0; i<particle.count; ++i )
        {
            vector2& vel = particle.velocity[i];
            if( vel.LengthSquared() > MAX_VELOCITY_MAG*MAX_VELOCITY_MAG )
            {
                vel.Normalize();
                vel *= MAX_VELOCITY_MAG;
            }
        }
    }
    virtual void Advect( float dt )
    {
        const vector2 gravdt = gravity * (mass.ascale*dt);
        const float invmdt = mass.invm*dt;
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.force[i] =
                gravdt + invmdt *
                ( particle.force[i] );
        }
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.velocity[i] += particle.force[i];
        }
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.position[i] += dt * ( particle.velocity[i] + 0.5f*particle.force[i] );
        }
        std::fill( particle.force.get() , particle.force.get() + particle.count , vector2(0.0f) );

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
        invN = 1.0f/(float)particle.count;

        const vector2 centerpos = std::accumulate( particle.position.get() , particle.position.get() + particle.count , vector2( 0.0f ) ) * invN;

        position = centerpos;
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.relative[i] = particle.position[i] - centerpos;
        }
        const float inertia = std::accumulate( particle.relative.get() , particle.relative.get() + particle.count , 0.0f ,
                []( const float cur , const vector2& rel )->float
                {
                    return cur + rel.LengthSquared();
                }
                ) * mass.m;
        invI = 1.0f/inertia * invIfactor;
    }

    void Advect( float dt ) override
    {
        const vector2 adt =
            ( std::accumulate( particle.force.get() , particle.force.get() + particle.count , vector2(0.0f) )*mass.invm*invN + gravity*mass.ascale )
            * dt;
        velocity += adt;
        position += ( velocity + 0.5f*adt ) * dt;

        const float alphadt = std::inner_product( particle.force.get() , particle.force.get() + particle.count , particle.relative.get() , 0.0f ,
                std::plus<float>() ,
                []( const vector2& force , const vector2& rel )->float
                {
                    return EH::Matrix::Cross( rel , force );
                }
                )
            * invI * dt;
        omega += alphadt;

        const float dtheta = ( omega + 0.5f*alphadt ) * dt;
        const vector2 dthetap = EH::Matrix::Complex::Complex( dtheta );

        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.relative[i] = EH::Matrix::Complex::Multiply( dthetap , particle.relative[i] );
            particle.position[i] = position + particle.relative[i];
            particle.velocity[i] = EH::Matrix::Cross( omega , particle.relative[i] );
        }
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
        for( RopePair& pair : ropepair )
        {
            vector2 rel = *pair.pos2 - *pair.pos1;
            const float gap = rel.Normalize() - rope_length;
            rel *= gap * FIX_FACTOR * 0.5f;
            *pair.pos1 += rel;
            *pair.pos2 -= rel;
        }
    }
    void Advect( float dt ) override
    {
        FluidGroup::Advect( dt );
        std::copy( particle.position.get() , particle.position.get() + particle.count , particle.relative.get() );
        for( unsigned int i=0; i<ITERATION; ++i )
        {
            Iterate();
        }
        const float invdtfactor = VELOCITY_FACTOR/dt;
        for( unsigned int i=0; i<particle.count; ++i )
        {
            particle.velocity[i] += ( particle.position[i] - particle.relative[i] ) * ( invdtfactor ) ;
        }
    }
};
