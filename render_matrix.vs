#version 130

#define RHO_FACTOR 0.1
#define SURFACE_FACTOR 0.3
#define VELOCITY_FACTOR 1.9
#define SKEW_THRESHOLD 0.7


in vec2 velocity;
in vec2 surfacenormal;
in float rho;

out mat2 outmat;

mat2 skew_matrix( vec2 axis , float factor )
{
    float lsq = dot( axis , axis );
    if( lsq <= 0.0 )
    {
        return mat2( 1.0 , 0.0 , 0.0 , 1.0 );
    }
    float k    = 1.0 + factor*sqrt( lsq );
    float invk = 1.0/k;
    float a    = ( k - invk ) / lsq;

    return
        mat2(
                a * axis.x*axis.x + invk , a * axis.x*axis.y , a * axis.x*axis.y , a * axis.y*axis.y + invk
            );
}

void main()
{
    float rho_scale = 1.0 + RHO_FACTOR * rho;

    outmat =  //vec4 ( 1.0 , 0.0 , 0.0 , 1.0 );
        0.2 *
        rho_scale *
        skew_matrix( velocity , VELOCITY_FACTOR ) *
        skew_matrix( vec2( surfacenormal.y , -surfacenormal.x ) , SURFACE_FACTOR );

    gl_Position = vec4( 0.0 );
}
