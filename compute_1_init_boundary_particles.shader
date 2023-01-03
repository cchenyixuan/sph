#version 460 compatibility

layout(std430, binding=0) buffer Particles{
    // particle inside domain with x, y, z, voxel_id; vx, vy, vz, mass; rho0, p0, rho, p; r, g, b, a
    mat4x4 Particle[];
};
layout(std430, binding=1) buffer BoundaryParticles{
    // particle at boundary with x, y, z, voxel_id; vx, vy, vz, mass; rho0, p0, rho, p; r, g, b, a
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

uniform int n_boundary_particle;  // particle number
uniform int n_voxel;  // voxel number
uniform float h;  // smooth radius

void AllocateBoundaryParticles(){
    // position of current particle focused
    vec3 particle_pos = BoundaryParticle[index][0].xyz;
    // for all voxels
    for(int i=0; i < n_voxel; ++i){
        // current voxel center position
        vec3 voxel_pos = Voxel[i*20+0][0].yzw;
        // current BoundaryParticle inside current voxel
        if(
            abs(particle_pos.x-voxel_pos.x) < h/2 &&
            abs(particle_pos.y-voxel_pos.y) < h/2 &&
            abs(particle_pos.z-voxel_pos.z) < h/2
            ){
                int j=0;
                // a voxel has maximum 96 slots for particles, whicl will be upgrated later  TODO
                while(j < 96){
                    // an empty slot is found
                    if(Voxel[i*20+2+j/16][(j%16)/4][(j%16)%4]==0){
                        // write index to slot
                        Voxel[i*20+2+j/16][(j%16)/4][(j%16)%4] = index_float;
                        // check if the operation above has been proformed properly.
                        // if not, try again
                        // check will be proformed after 100 operations
                        int counter = 0;
                        while (counter < 100) {
                          // do nothing
                            counter += 1;
                        }
                        if(Voxel[i*20+2+j/16][(j%16)/4][(j%16)%4]==index_float){
                            // well-performed
                            // set particle's voxel id to i+1 at BoundaryParticle[index][0].w
                            BoundaryParticle[index][0].w = float(i+1);
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
    AllocateBoundaryParticles();
}
