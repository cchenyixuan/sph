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

uniform int id;
uniform int n_voxel;  // voxel number

void main() {
    for(int i=0; i<n_voxel; ++i){
        Voxel[i*20+19][3].xyzw = vec4(0.5, 0.0, 0.0, 0.3);
    }
    barrier();
    // set blue itself
    Voxel[id*20+19][3].xyzw = vec4(0.0, 0.5, 1.0, 1.0);
    // set neighborhoods green
    for(int i=4; i<32; ++i){
        int voxel_id = int(round(Voxel[id*20+i/16][(i%16)/4][(i%16)%4]));
        if(voxel_id!=0){
            Voxel[(voxel_id-1)*20+19][3].xyzw = vec4(0.0, 1.0, 0.5, 0.5);
        }
    }

}
