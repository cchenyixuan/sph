#version 460 compatibility

layout(std430, binding=0) buffer Particles{
    // particle inside domain with x, y, z, voxel_id; vx, vy, vz, mass; wx, wy, wz, rho; ax, ay, az, P;
    mat4x4 Particle[];
};
layout(std430, binding=1) buffer BoundaryParticles{
    // particle at boundary with x, y, z, voxel_id; vx, vy, vz, mass; wx, wy, wz, rho; ax, ay, az, P;
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


layout(local_size_x=1, local_size_y=1, local_size_z=1) in;

uint gid = gl_GlobalInvocationID.x;
int particle_index = int(gid)+1;
float particle_index_float = float(particle_index);

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

uniform int id;


const int voxel_memory_length = 2912;
const int voxel_block_size = 960;

void set_voxel_data(){
    Voxel[id*voxel_memory_length+31] = 100;
    for(int i=4; i<30; ++i){
        if(Voxel[id*voxel_memory_length+i]!=0 && (i>=4 && i<=9)){
            Voxel[(Voxel[id*voxel_memory_length+i]-1)*voxel_memory_length+31] = i;
        }

    }
}

void main() {
    // set 0 everyone
    for(int i=0; i<n_voxel; ++i){
        int c = Voxel[i*voxel_memory_length+31];
        atomicAdd(Voxel[i*voxel_memory_length+31], -c);
    }
    barrier();
    // set 1 itself if particle inside
    for(int i=0; i<n_voxel; ++i){
        for(int j=0; j<voxel_block_size; ++j){
            if(Voxel[i*voxel_memory_length+32+j]>0){
                int c = Voxel[i*voxel_memory_length+31];
                atomicAdd(Voxel[i*voxel_memory_length+31], -c+1);
            }
        }

    }
    barrier();
    set_voxel_data();
}
