#version 130

in vec2 rel;
in float dv_white;

const vec3 ambient = vec3( 0.1 , 0.2 , 0.4 );
const vec3 dv_spec = vec3( 0.3 , 0.3 , 0.4 );

void main()
{
    float h = 1.0 - dot( rel , rel );
    gl_FragColor = vec4( ambient + dv_white * dv_spec , h );
    /*gl_FragColor = vec4( vec3(1.0) , 1.0-dot(rel,rel) );*/
    /*gl_FragColor = vec4( 1.0 );*/
}
