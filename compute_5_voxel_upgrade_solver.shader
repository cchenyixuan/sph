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
    // all zeros, length equals voxel numbers
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
int voxel_index = int(gid)+1;
float voxel_index_float = float(voxel_index);

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius



void UpgradeVoxel(){
    // create a counter to check all out buffer has been considered
    int out_counter = 0;
    for(int i=0; i<96; ++i){
        // read out buffer, out buffer is continuous
        int out_particle_id = Voxel[(voxel_index-1)*320+32+96+i%96];  // starts from 1
        // if zero found, iteration should stop
        if(out_particle_id==0){break;}
        for(int j=0; j<96; ++j){
            // read inside buffer and check if particle id match
            int inside_particle_id = Voxel[(voxel_index-1)*320+32+j%96];
            if(inside_particle_id==out_particle_id){
                // erase voxel particle slot and add counter
                Voxel[(voxel_index-1)*320+32+j%96] = 0;
                out_counter += 1;
                // break
                break;
            }
        }
    }
    // check if out_counter == VoxelParticleOutNumber[voxel_index-1]
    if(out_counter==VoxelParticleOutNumber[voxel_index-1]){
        // verified
    }
    else{
        // particle number not match! debug here
    }
    // create a counter to check all in buffer has been considered
    int in_counter = 0;
    for(int i=0; i<96; ++i){
        // read in buffer, in buffer is continuous
        int in_particle_id = Voxel[(voxel_index-1)*320+32+96+96+i%96];  // starts from 1
        // if zero found, iteration should stop
        if(in_particle_id==0){break;}
        for(int j=0; j<96; ++j){
            // read inside buffer and check if a slot is found
            if(Voxel[(voxel_index-1)*320+32+j%96]==0){
                // insert voxel particle slot and add counter
                Voxel[(voxel_index-1)*320+32+j%96] = in_particle_id;
                in_counter += 1;
                // break
                break;
            }
        }
    }
    // check if in_counter == VoxelParticleInNumber[voxel_index-1]
    if(in_counter==VoxelParticleInNumber[voxel_index-1]){
        // verified
    }
    else{
        // particle number not match! debug here
    }
    // code above could be modified to have O(n) time complexity
    // re-arrange inside buffer using 2 pointers
    int ptr1 = 0;  // slow
    int ptr2 = 0;  // fast
    while(ptr2<96){
        if(Voxel[(voxel_index-1)*320+32+ptr2%96]==0){
            // empty slot found
            ptr2 += 1;
        }
        else if(Voxel[(voxel_index-1)*320+32+ptr2%96]!=0 && ptr1!=ptr2){
            // filled slot found
            Voxel[(voxel_index-1)*320+32+ptr1%96] = Voxel[(voxel_index-1)*320+32+ptr2%96];
            Voxel[(voxel_index-1)*320+32+ptr2%96] = 0;
            ptr1 += 1;
        }
        else{
            ptr1 += 1;
            ptr2 += 1;
        }
    }
    // above code persists points order
    // following code breaks points order but could have better performance
    /*
    int ptr1 = 0;
    int ptr2 = 95;
    while(ptr1<ptr2){
        if(Voxel[(voxel_index-1)*320+32+ptr1%96]==0){
            if(Voxel[(voxel_index-1)*320+32+ptr2%96]!=0){
                Voxel[(voxel_index-1)*320+32+ptr1%96] = Voxel[(voxel_index-1)*320+32+ptr2%96];
                Voxel[(voxel_index-1)*320+32+ptr2%96] = 0;
            }
            else{
                ptr2 -= 1;
            }
        }
        else{
            ptr1 += 1;
        }
    }
    */

    // all out particles have been checked, clear out buffer
    for(int i=0; i<96; ++i){
        Voxel[(voxel_index-1)*320+32+96+i%96] = 0;
    }
    // all in particles have been checked, clear in buffer
    for(int i=0; i<96; ++i){
        Voxel[(voxel_index-1)*320+32+96+96+i%96] = 0;
    }
    // re-calculate VoxelParticleNumber and clear VoxelParticleInNumber and VoxelParticleOutNumber
    VoxelParticleNumber[voxel_index-1] += -VoxelParticleOutNumber[voxel_index-1]+VoxelParticleInNumber[voxel_index-1];
    VoxelParticleOutNumber[voxel_index-1] = 0;
    VoxelParticleInNumber[voxel_index-1] = 0;
}

void main() {
    UpgradeVoxel();
}