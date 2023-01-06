#version 460 compatibility

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

layout(local_size_x=1, local_size_y=1, local_size_z=1) in;

uint gid = gl_GlobalInvocationID.x;
int particle_index = int(gid)+1;
float particle_index_float = float(particle_index);

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

uniform int id;

void main() {
    // set 0 everyone
    for(int i=0; i<n_voxel; ++i){
        int c = Voxel[i*320+319];
        atomicAdd(Voxel[i*320+319], -c);
    }
    barrier();
    // set 1 itself
    int c = Voxel[id*320+319];
    atomicAdd(Voxel[id*320+319], -c+1);
    barrier();
    // set 2 neighborhoods
    for(int i=4; i<32; ++i){
        int voxel_id = Voxel[id*320+i];
        if(voxel_id!=0){
            int c = Voxel[(voxel_id-1)*320+319];
            atomicAdd(Voxel[(voxel_id-1)*320+319], -c+2);
        }
    }

}
