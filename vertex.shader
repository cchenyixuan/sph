#version 460 core

layout(location=0) in int v_index; // vertex id
out vec4 v_color; // color output

layout(std430, binding=0) buffer Particles{
    // particle inside domain with x, y, z, voxel_id; vx, vy, vz, mass; rho0, p0, rho, p; r, g, b, a
    mat4x4 Particle[];
};
layout(std430, binding=1) buffer BoundaryParticles{
    // particle at boundary with x, y, z, voxel_id; vx, vy, vz, mass; rho0, p0, rho, p; r, g, b, a
    mat4x4 BoundaryParticle[];
};
layout(std430, binding=2) coherent buffer Voxels{
    // each voxel has 20 mat44 and first 2 matrices contains its id, x_offset of h, y_offset of h, z_offset of h; and neighborhood voxel ids
    // other 18 matrices containing current-indoor-particle-ids, particles getting out and particles stepping in
    // matrices are changed into integer arrays to apply atomic operations, first 32 integers for first 2 matrices and one voxel costs 320 integers
    int Voxel[];
};
layout(std430, binding=3) coherent buffer VoxelParticleNumbers{
    int VoxelParticleNumber[];
};

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

uniform mat4 projection;
uniform mat4 view;


void main() {
    gl_Position = projection*view*vec4(Particle[v_index][0].xyz, 1.0); // set vertex position, w=1.0
    int voxel_id = int(round(Particle[v_index][0].w));
    vec3 voxel_center = vec3(float(Voxel[(voxel_id-1)*320+1])*h, float(Voxel[(voxel_id-1)*320+2])*h, float(Voxel[(voxel_id-1)*320+3])*h);
    float l = length(Particle[v_index][3].xyz);
    v_color = vec4(abs(Particle[v_index][3].x)/l, abs(Particle[v_index][3].y)/l, abs(Particle[v_index][3].z)/l, 0.3); // set output color by its voxel id
    //v_color = vec4(abs(sin(float(voxel_id/2))), abs(cos(float(voxel_id/3))), abs(sin(float(voxel_id/5))), 0.3);
}