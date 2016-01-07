#version 130

uniform mat4 u_projection;

in vec2 position;
in vec4 modelview;
in float dv_white_in;

out vec2 rel;
out float dv_white;

void main()
{
    dv_white = max( dv_white_in , 0.0 );
    vec2 vert = vec2( (gl_VertexID&1)*2 - 1 , (gl_VertexID&2)-1 );
    rel = vert;
    gl_Position = u_projection * vec4( modelview.xy * vert.x + modelview.zw * vert.y + position , 0.0 , 1.0 );
}
