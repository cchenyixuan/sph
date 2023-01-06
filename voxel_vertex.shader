#version 460 core

layout(location=0) in int v_index; // vertex id
out GeometryOutput{
    vec4 v_pos;
    vec4 v_color;
}g_out;



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

uniform mat4 projection;
uniform mat4 view;

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius


void main() {
    g_out.v_pos = vec4(float(Voxel[v_index*320+1])*h, float(Voxel[v_index*320+2])*h, float(Voxel[v_index*320+3])*h, 1.0);
    vec4 color = vec4(0.0);
    if(Voxel[v_index*320+319]==0){color = vec4(0.5, 0.0, 0.0, 0.3);}
    if(Voxel[v_index*320+319]==1){color = vec4(0.0, 0.5, 1.0, 1.0);}
    if(Voxel[v_index*320+319]==2){color = vec4(0.0, 1.0, 0.5, 0.5);}
    g_out.v_color = color;

}
