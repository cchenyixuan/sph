#version 460 compatibility

layout(std430, binding=0) buffer Particles{
    // particle inside domain with x, y, z, 0, vx, vy, vz, 0, ...
    mat4x4 Particle[];
};
layout(std430, binding=1) buffer BoundaryParticles{
    // particle at boundary with x, y, z, 0, vx, vy, vz, 0, ...
    mat4x4 BoundaryParticle[];
};
layout(std430, binding=2) buffer Voxels{
    // each voxel has 20 mat44 and first 2 matrices contains its id, loc_x, loc_y, loc_z and neighborhood voxel ids
    // other 18 matrices containing current-indoor-particle-ids, particles getting out and particles stepping in
    mat4x4 Voxel[];
};

layout(local_size_x=1, local_size_y=1, local_size_z=1) in;

uint gid = gl_GlobalInvocationID.x;
int index = int(gid);
float index_float = float(index);

uniform int n_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

void AllocateParticles(){
    // position of current particle focused
    vec3 particle_pos = Particle[index][0].xyz;
    // for all voxels
    for(int i=0; i < n_voxel; ++i){
        // current voxel center position
        vec3 voxel_pos = Voxel[i*20+0][0].yzw;
        // current particle inside current voxel
        if(
            abs(particle_pos.x-voxel_pos.x) < h/2 &&
            abs(particle_pos.y-voxel_pos.y) < h/2 &&
            abs(particle_pos.z-voxel_pos.z) < h/2
            ){
                int j=0;
                // a voxel has maximum 96 slots for particles, whicl will be upgrated later  TODO
                while(j < 96){
                    // an empty slot is found
                    if(Voxel[i*20+2+j/16][j%16/4][j%16%4]==0){
                        // write index to slot
                        Voxel[i*20+2+j/16][j%16/4][j%16%4] = index_float;
                        // check if the operation above has been proformed properly.
                        // if not, try again
                        // check will be proformed after 0.001s
                        float startTime = time;
                        while (time - startTime < 0.001) {
                          // do nothing
                        }
                        if(Voxel[i*20+2+j/16][j%16/4][j%16%4]==index_float){
                            // well-performed
                            break;
                        }
                        else{
                            // conflict found, try again with j not changed
                            continue;
                        }
                    }
                    // move on to the next slot
                    else{
                        j += 1;
                    }
                }
        };
    }
}

void main() {
    AllocateParticles();
}
