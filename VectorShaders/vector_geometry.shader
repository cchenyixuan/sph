#version 460 core

layout (points) in;
layout (line_strip, max_vertices = 2) out;
in GeometryOutput{
    vec4 v_pos;
    vec4 v_color;
}g_in[];
out vec4 v_color;

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

uniform mat4x4 projection;
uniform mat4x4 view;


void CreateVector(){
    // 2 vertices
    vec4 p1 = g_in[0].v_pos;
    vec4 p2 = vec4(g_in[0].v_pos.xyz+g_in[0].v_color.xyz, 1.0);
    // vertex color
    v_color = g_in[0].v_color;
    // 2 vertices emittion to generate a vector line
    gl_Position = projection*view*p1;
    EmitVertex();
    gl_Position = projection*view*p2;
    EmitVertex();
    // end of line-strip
    EndPrimitive();
}

void main() {
    CreateVector();

}
