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
    // each voxel has 182 mat44 and first 2 matrices contains its id, x_offset of h, y_offset of h, z_offset of h; and neighborhood voxel ids
    // other 180 matrices containing current-indoor-particle-ids, particles getting out and particles stepping in
    // matrices are changed into integer arrays to apply atomic operations, first 32 integers for first 2 matrices and one voxel costs 2912 integers
    int Voxel[];
};
layout(std430, binding=3) coherent buffer VoxelParticleNumbers{
    int VoxelParticleNumber[];
};
layout(std430, binding=4) coherent buffer VoxelParticleInNumbers{
    int VoxelParticleInNumber[];
};
layout(std430, binding=5) coherent buffer VoxelParticleOutNumbers{
    int VoxelParticleOutNumber[];
};


uniform mat4 projection;
uniform mat4 view;

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

const int voxel_memory_length = 2912;
const int voxel_block_size = 960;


void main() {
    g_out.v_pos = vec4(float(Voxel[v_index*voxel_memory_length+1])*h, float(Voxel[v_index*voxel_memory_length+2])*h, float(Voxel[v_index*voxel_memory_length+3])*h, 1.0);
    vec4 color = vec4(0.0);
    if(Voxel[v_index*voxel_memory_length+31]==0){color = vec4(1.0, 1.0, 0.0, 0.3);}
    else if(Voxel[v_index*voxel_memory_length+31]==1){color = vec4(0.0, 0.5, 1.0, 0.5);}
    else if(Voxel[v_index*voxel_memory_length+31]==2){color = vec4(0.0, 1.0, 0.5, 0.5);}
    else if(Voxel[v_index*voxel_memory_length+31]==3){color = vec4(1.0, 1.0, 0.2, 0.5);}
    else if(Voxel[v_index*voxel_memory_length+31]==4){
        // left
        color = vec4(0.5, 0.5, 0.0, 1.0);
    }
    else if(Voxel[v_index*voxel_memory_length+31]==8){
        // right
        color = vec4(1.0, 0.0, 0.0, 1.0);
    }

    g_out.v_color = color;

}
