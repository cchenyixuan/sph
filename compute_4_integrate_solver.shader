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



// function definitions

const float PI = 3.141592653589793;
const float REST_DENS = 1000.0;
const float EOS_CONST = 2000.0;
const float VISC = 10.0;
const float DELTA_T = 0.01;
const float EPS = 0.01;  // to control movement not too large


float h2 = h * h;

void EulerMethod(){
    // calculate future position
    //   move =             P_velocity           *   dt   +      P_acceleration  *     dt/2
    vec3 move = Particle[particle_index-1][1].xyz*DELTA_T + Particle[particle_index-1][3].xyz*DELTA_T*DELTA_T/2;
    // estimate future position
    //   future_pos =             P_position            + move
    vec3 future_pos = Particle[particle_index-1][0].xyz + move;
    // which voxel this particle would go
    // current voxel id
    int voxel_id = int(round(Particle[particle_index-1][0].w));  // starts from 1
    // current voxel center
    vec3 voxel_center = vec3(float(Voxel[(voxel_id-1)*320+1])*h, float(Voxel[(voxel_id-1)*320+2])*h, float(Voxel[(voxel_id-1)*320+3])*h);
    // x axis flag (-1: left, 0: current, 1: right)
    int x_axis_flag;
    if     (voxel_center.x+h/2<=future_pos.x){x_axis_flag=1;}
    else if(voxel_center.x-h/2<=future_pos.x && future_pos.x<voxel_center.x+h/2){x_axis_flag=0;}
    else if(future_pos.x<voxel_center.x-h/2){x_axis_flag=-1;}
    // y axis flag (-1: below, 0: current, 1: above)
    int y_axis_flag;
    if     (voxel_center.y+h/2<=future_pos.y){y_axis_flag=1;}
    else if(voxel_center.y-h/2<=future_pos.y && future_pos.y<voxel_center.y+h/2){y_axis_flag=0;}
    else if(future_pos.y<voxel_center.y-h/2){y_axis_flag=-1;}
    // z axis flag (-1: back, 0: current, 1: front)
    int z_axis_flag;
    if     (voxel_center.z+h/2<=future_pos.z){z_axis_flag=1;}
    else if(voxel_center.z-h/2<=future_pos.z && future_pos.z<voxel_center.z+h/2){z_axis_flag=0;}
    else if(future_pos.z<voxel_center.z-h/2){z_axis_flag=-1;}
    // identify which voxel the particle will go
    ivec3 flag = ivec3(x_axis_flag, y_axis_flag, z_axis_flag);
    // remain inside
    if(flag==ivec3(0, 0, 0)){
        // do nothing
    }
    // particle goes to other voxel
    else{
        // set voxel out buffer
        // o starts from 0 and max is 95, this will add 1 to o and return o before addition happends
        int o = atomicAdd(VoxelParticleOutNumber[voxel_id-1], 1);
        barrier();
        atomicAdd(Voxel[(voxel_id-1)*320+32+96+o%96], particle_index);  // starts from 1 (domain particle)
        barrier();
        // left
        if(flag==ivec3(-1, 0, 0)){
            // set voxel in buffer
            // particle goes to left and the voxel left to the current one is Voxel[(voxel_id-1)*320+4] by our standard
            int left_voxel_id = Voxel[(voxel_id-1)*320+4];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right
        else if(flag==ivec3(1, 0, 0)){
            // set voxel in buffer
            // particle goes to right and the voxel right to the current one is Voxel[(voxel_id-1)*320+5] by our standard
            int right_voxel_id = Voxel[(voxel_id-1)*320+5];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // down
        else if(flag==ivec3(0, -1, 0)){
            // set voxel in buffer
            // particle goes to down and the voxel down to the current one is Voxel[(voxel_id-1)*320+6] by our standard
            int down_voxel_id = Voxel[(voxel_id-1)*320+6];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[down_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(down_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // up
        else if(flag==ivec3(0, 1, 0)){
            // set voxel in buffer
            // particle goes to up and the voxel up to the current one is Voxel[(voxel_id-1)*320+7] by our standard
            int up_voxel_id = Voxel[(voxel_id-1)*320+7];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[up_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(up_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // back
        else if(flag==ivec3(0, 0, -1)){
            // set voxel in buffer
            // particle goes to back and the voxel back to the current one is Voxel[(voxel_id-1)*320+8] by our standard
            int back_voxel_id = Voxel[(voxel_id-1)*320+8];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // front
        else if(flag==ivec3(0, 0, 1)){
            // set voxel in buffer
            // particle goes to front and the voxel front to the current one is Voxel[(voxel_id-1)*320+9] by our standard
            int front_voxel_id = Voxel[(voxel_id-1)*320+9];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_down
        else if(flag==ivec3(-1, -1, 0)){
            // set voxel in buffer
            // particle goes to left_down and the voxel left_down to the current one is Voxel[(voxel_id-1)*320+10] by our standard
            int left_down_voxel_id = Voxel[(voxel_id-1)*320+10];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_down_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_down_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_up
        else if(flag==ivec3(-1, 1, 0)){
            // set voxel in buffer
            // particle goes to left_up and the voxel left_up to the current one is Voxel[(voxel_id-1)*320+11] by our standard
            int left_up_voxel_id = Voxel[(voxel_id-1)*320+11];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_up_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_up_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_down
        else if(flag==ivec3(1, -1, 0)){
            // set voxel in buffer
            // particle goes to right_down and the voxel right_down to the current one is Voxel[(voxel_id-1)*320+12] by our standard
            int right_down_voxel_id = Voxel[(voxel_id-1)*320+12];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_down_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_down_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_up
        else if(flag==ivec3(1, 1, 0)){
            // set voxel in buffer
            // particle goes to right_up and the voxel right_up to the current one is Voxel[(voxel_id-1)*320+13] by our standard
            int right_up_voxel_id = Voxel[(voxel_id-1)*320+13];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_up_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_up_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_back
        else if(flag==ivec3(-1, 0, -1)){
            // set voxel in buffer
            // particle goes to left_back and the voxel left_back to the current one is Voxel[(voxel_id-1)*320+14] by our standard
            int left_back_voxel_id = Voxel[(voxel_id-1)*320+14];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_front
        else if(flag==ivec3(-1, 0, 1)){
            // set voxel in buffer
            // particle goes to left_front and the voxel left_front to the current one is Voxel[(voxel_id-1)*320+15] by our standard
            int left_front_voxel_id = Voxel[(voxel_id-1)*320+15];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_back
        else if(flag==ivec3(1, 0, -1)){
            // set voxel in buffer
            // particle goes to right_back and the voxel right_back to the current one is Voxel[(voxel_id-1)*320+16] by our standard
            int right_back_voxel_id = Voxel[(voxel_id-1)*320+16];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_front
        else if(flag==ivec3(1, 0, 1)){
            // set voxel in buffer
            // particle goes to right_front and the voxel right_front to the current one is Voxel[(voxel_id-1)*320+17] by our standard
            int right_front_voxel_id = Voxel[(voxel_id-1)*320+17];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // down_back
        else if(flag==ivec3(0, -1, -1)){
            // set voxel in buffer
            // particle goes to down_back and the voxel down_back to the current one is Voxel[(voxel_id-1)*320+18] by our standard
            int down_back_voxel_id = Voxel[(voxel_id-1)*320+18];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[down_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(down_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // down_front
        else if(flag==ivec3(0, -1, 1)){
            // set voxel in buffer
            // particle goes to down_front and the voxel down_front to the current one is Voxel[(voxel_id-1)*320+19] by our standard
            int down_front_voxel_id = Voxel[(voxel_id-1)*320+19];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[down_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(down_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // up_back
        else if(flag==ivec3(0, 1, -1)){
            // set voxel in buffer
            // particle goes to up_back and the voxel up_back to the current one is Voxel[(voxel_id-1)*320+20] by our standard
            int up_back_voxel_id = Voxel[(voxel_id-1)*320+20];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[up_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(up_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // up_front
        else if(flag==ivec3(0, 1, 1)){
            // set voxel in buffer
            // particle goes to up_front and the voxel up_front to the current one is Voxel[(voxel_id-1)*320+21] by our standard
            int up_front_voxel_id = Voxel[(voxel_id-1)*320+21];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[up_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(up_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_down_back
        else if(flag==ivec3(-1, -1, -1)){
            // set voxel in buffer
            // particle goes to left_down_back and the voxel left_down_back to the current one is Voxel[(voxel_id-1)*320+22] by our standard
            int left_down_back_voxel_id = Voxel[(voxel_id-1)*320+22];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_down_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_down_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_down_front
        else if(flag==ivec3(-1, -1, 1)){
            // set voxel in buffer
            // particle goes to left_down_front and the voxel left_down_front to the current one is Voxel[(voxel_id-1)*320+23] by our standard
            int left_down_front_voxel_id = Voxel[(voxel_id-1)*320+23];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_down_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_down_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_up_back
        else if(flag==ivec3(-1, 1, -1)){
            // set voxel in buffer
            // particle goes to left_up_back and the voxel left_up_back to the current one is Voxel[(voxel_id-1)*320+24] by our standard
            int left_up_back_voxel_id = Voxel[(voxel_id-1)*320+24];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_up_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_up_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // left_up_front
        else if(flag==ivec3(-1, 1, 1)){
            // set voxel in buffer
            // particle goes to left_up_front and the voxel left_up_front to the current one is Voxel[(voxel_id-1)*320+25] by our standard
            int left_up_front_voxel_id = Voxel[(voxel_id-1)*320+25];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[left_up_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(left_up_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_down_back
        else if(flag==ivec3(1, -1, -1)){
            // set voxel in buffer
            // particle goes to right_down_back and the voxel right_down_back to the current one is Voxel[(voxel_id-1)*320+26] by our standard
            int right_down_back_voxel_id = Voxel[(voxel_id-1)*320+26];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_down_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_down_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_down_front
        else if(flag==ivec3(1, -1, 1)){
            // set voxel in buffer
            // particle goes to right_down_front and the voxel right_down_front to the current one is Voxel[(voxel_id-1)*320+27] by our standard
            int right_down_front_voxel_id = Voxel[(voxel_id-1)*320+27];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_down_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_down_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_up_back
        else if(flag==ivec3(1, 1, -1)){
            // set voxel in buffer
            // particle goes to right_up_back and the voxel right_up_back to the current one is Voxel[(voxel_id-1)*320+28] by our standard
            int right_up_back_voxel_id = Voxel[(voxel_id-1)*320+28];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_up_back_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_up_back_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
        // right_up_front
        else if(flag==ivec3(1, 1, 1)){
            // set voxel in buffer
            // particle goes to right_up_front and the voxel right_up_front to the current one is Voxel[(voxel_id-1)*320+29] by our standard
            int right_up_front_voxel_id = Voxel[(voxel_id-1)*320+29];  // starts from 1
            int i = atomicAdd(VoxelParticleInNumber[right_up_front_voxel_id-1], 1);
            barrier();
            atomicAdd(Voxel[(right_up_front_voxel_id-1)*320+32+96+96+i%96], particle_index);  // starts from 1 (domain particle)
            barrier();
        }
    }

    // particle position and velocity will be set and particle acceleration will be erased
    Particle[particle_index-1][0].xyz = future_pos;
    Particle[particle_index-1][1].xyz += Particle[particle_index-1][3].xyz*DELTA_T;
    Particle[particle_index-1][3].xyz = vec3(0.0);
}

void main() {
    EulerMethod();
}
