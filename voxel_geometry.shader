#version 460 core

layout (points) in;
layout (triangle_strip, max_vertices = 14) out;
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


void CreateCube(){
    // center of voxel
    vec4 center = g_in[0].v_pos;
    // 8 vertices
    vec4 p1 = vec4(center.x+h/2, center.y-h/2, center.z-h/2, 1.0);
    vec4 p2 = vec4(p1.x, p1.y+h, p1.z, 1.0);
    vec4 p3 = vec4(p1.x-h, p1.y, p1.z, 1.0);
    vec4 p4 = vec4(p3.x, p3.y+h, p3.z, 1.0);
    vec4 p5 = vec4(p1.x, p1.y, p1.z+h, 1.0);
    vec4 p6 = vec4(p5.x, p5.y+h, p5.z, 1.0);
    vec4 p7 = vec4(p4.x, p4.y, p4.z+h, 1.0);
    vec4 p8 = vec4(p3.x, p3.y, p3.z+h, 1.0);

    // vertex color
    v_color = g_in[0].v_color;
    // 12 vertices emittion to generate a full cube in order 4-3-7-8-5-3-1-4-2-7-6-5-2-1
    gl_Position = projection*view*p4;
    EmitVertex();
    gl_Position = projection*view*p3;
    EmitVertex();
    gl_Position = projection*view*p7;
    EmitVertex();
    gl_Position = projection*view*p8;
    EmitVertex();
    gl_Position = projection*view*p5;
    EmitVertex();
    gl_Position = projection*view*p3;
    EmitVertex();
    gl_Position = projection*view*p1;
    EmitVertex();
    gl_Position = projection*view*p4;
    EmitVertex();
    gl_Position = projection*view*p2;
    EmitVertex();
    gl_Position = projection*view*p7;
    EmitVertex();
    gl_Position = projection*view*p6;
    EmitVertex();
    gl_Position = projection*view*p5;
    EmitVertex();
    gl_Position = projection*view*p2;
    EmitVertex();
    gl_Position = projection*view*p1;
    EmitVertex();

    // end of triangle-strip
    EndPrimitive();


}

void main() {
    CreateCube();

}
