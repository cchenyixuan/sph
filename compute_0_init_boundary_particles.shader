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
layout(std430, binding=6) coherent buffer GlobalStatus{
    // simulation global settings and status such as max velocity etc.
    // [n_particle, n_boundary_particle, n_voxel, voxel_memory_length, voxel_block_size, h_p, h_q, r_p, r_q, max_velocity_n-times_than_r, rest_dense, eos_constant]
    int Status[];
};
layout(std430, binding=7) buffer ParticlesSubData{
    // particle inside domain has additional data: t_transfer.xyz, 0.0, 0.0...;
    mat4x4 ParticleSubData[];
};


layout(local_size_x=1, local_size_y=1, local_size_z=1) in;

uint gid = gl_GlobalInvocationID.x;
int particle_index = int(gid)+1;
float particle_index_float = float(particle_index);

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

const int voxel_memory_length = 2912;
const int voxel_block_size = 960;



void AllocateBoundaryParticles(){
    // position of current particle focused
    vec3 particle_pos = BoundaryParticle[particle_index-1][0].xyz;
    // for all voxels
    for(int i=0; i < n_voxel; ++i){
        // current voxel center position
        vec3 voxel_pos = vec3(float(Voxel[i*voxel_memory_length+1])*h, float(Voxel[i*voxel_memory_length+2])*h, float(Voxel[i*voxel_memory_length+3])*h);
        // current particle inside current voxel (vx-2/h<=px<vx+2/h)
        if(
            voxel_pos.x-h/2<=particle_pos.x && particle_pos.x<voxel_pos.x+h/2 &&
            voxel_pos.y-h/2<=particle_pos.y && particle_pos.y<voxel_pos.y+h/2 &&
            voxel_pos.z-h/2<=particle_pos.z && particle_pos.z<voxel_pos.z+h/2
            ){
                // one particle found inside current voxel, get its slot id (start from 0) and add 1 to next slot id
                int c = atomicAdd(VoxelParticleNumber[i], 1);
                barrier();
                // set slot with index value
                atomicAdd(Voxel[i*voxel_memory_length+32+c%voxel_block_size], -particle_index);  // starts from -1 (boundary particle)
                barrier();
                // set particle's voxel id
                BoundaryParticle[particle_index-1][0].w = float(i+1);  // starts from 1.0

        };
    }
}

void main() {
    AllocateBoundaryParticles();
}
