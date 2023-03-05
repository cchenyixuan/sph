#version 460 core

layout(location=0) in int v_index; // vertex id
out vec4 v_color; // color output

layout(std430, binding=0) buffer Particles{
    // particle inside domain with x, y, z, voxel_id; vx, vy, vz, mass; wx, wy, wz, rho; ax, ay, az, P;
    mat4x4 Particle[];
};
layout(std430, binding=1) buffer ParticlesSubData{
    // particle inside domain has additional data: t_transfer.xyz, 0.0, 0.0...;
    mat4x4 ParticleSubData[];
};
layout(std430, binding=2) buffer BoundaryParticles{
    // particle at boundary with x, y, z, voxel_id; vx, vy, vz, mass; wx, wy, wz, rho; ax, ay, az, P;
    mat4x4 BoundaryParticle[];
};
layout(std430, binding=3) coherent buffer Voxels{
    // each voxel has 182 mat44 and first 2 matrices contains its id, x_offset of h, y_offset of h, z_offset of h; and neighborhood voxel ids
    // other 180 matrices containing current-indoor-particle-ids, particles getting out and particles stepping in
    // matrices are changed into integer arrays to apply atomic operations, first 32 integers for first 2 matrices and one voxel costs 2912 integers
    int Voxel[];
};
layout(std430, binding=4) coherent buffer Voxels2{
    // each voxel has 182 mat44 and first 2 matrices contains its id, x_offset of h, y_offset of h, z_offset of h; and neighborhood voxel ids
    // other 180 matrices containing current-indoor-particle-ids, particles getting out and particles stepping in
    // matrices are changed into integer arrays to apply atomic operations, first 32 integers for first 2 matrices and one voxel costs 2912 integers
    int Voxel2[];
};
layout(std430, binding=5) coherent buffer VoxelParticleNumbers{
    int VoxelParticleNumber[];
};
layout(std430, binding=6) coherent buffer VoxelParticleInNumbers{
    int VoxelParticleInNumber[];
};
layout(std430, binding=7) coherent buffer VoxelParticleOutNumbers{
    int VoxelParticleOutNumber[];
};
layout(std430, binding=8) buffer GlobalStatus{
    // simulation global settings and status such as max velocity etc.
    // [n_particle, n_boundary_particle, n_voxel, voxel_memory_length, voxel_block_size, h_p, h_q, r_p, r_q, max_velocity_n-times_than_r, rest_dense, eos_constant, t_p, t_q, v_p, v_q, c_p, c_q, a_p, a_q]
    int StatusInt[];
};
layout(std430, binding=9) buffer GlobalStatus2{
    // simulation global settings and status such as max velocity etc.
    // [n_particle, n_boundary_particle, n_voxel, voxel_memory_length, voxel_block_size, h_p, h_q, r_p, r_q, max_velocity_n-times_than_r, rest_dense, eos_constant, t_p, t_q, v_p, v_q, c_p, c_q, a_p, a_q]
    float StatusFloat[];
};


/*
const float PI = 3.141592653589793;
const int n_boundary_particle = Status[0];
const int n_particle = Status[1];
const int n_voxel = Status[2];
const float h = float(Status[5])/float(Status[6]);
const float r = float(Status[7])/float(Status[8]);
const int voxel_memory_length = Status[3];
const int voxel_block_size = Status[4];
const float rest_dense = float(Status[10]);
const float eos_constant = float(Status[11]);
const float delta_t = float(Status[12])/float(Status[13]);
const float viscosity = float(Status[14])/float(Status[15]);
const float cohesion = float(Status[16])/float(Status[17]);
const float adhesion = float(Status[18])/float(Status[19]);
*/
const float PI = 3.141592653589793;
const int n_boundary_particle = StatusInt[0];
const int n_particle = StatusInt[1];
const int n_voxel = StatusInt[2];
const float h = StatusFloat[0];
const float r = StatusFloat[1];
const int voxel_memory_length = 2912;
const int voxel_block_size = 960;
const float rest_dense = 1000;
const float eos_constant = 276.571;
const float delta_t = StatusFloat[2];
const float viscosity = StatusFloat[3];
const float cohesion = StatusFloat[4];
const float adhesion = StatusFloat[5];

uniform mat4 projection;
uniform mat4 view;

uniform int color_type;



vec3 get_color_gradient(float ratio, float range){
    // ratio in range 0.9-1.1
    /*
        0.9 --> purple(0.5, 0.0, 1.0)
        0.93 --> blue(0.0, 0.0, 1.0)
        0.96 --> cyan(0.0, 1.0, 1.0)
        1.0 --> green(0.0, 1.0, 0.0)
        1.03 --> yellow(1.0, 1.0, 0.0)
        1.06 --> orange(1.0, 0.5, 0.0)
        1.1 --> red(1.0, 0.0, 0.0)
    */
    float red = 0.0;
    if(ratio<1.0){red=min(max(-15*(ratio-(1.0-range*2/3)), 0.0), 1.0);}
    else if(ratio>=1.0){red=min(max(3/range*(ratio-1.0), 0.0), 1.0);}
    float green = 0.0;
    if(ratio<1.0){green=min(max(30*(ratio-(1.0-range*2/3)), 0.0), 1.0);}
    else if(ratio>=1.0){green=min(max(-1.5/range*(ratio-(1.0+range)), 0.0), 1.0);}
    float blue = min(max(-3/range*(ratio-1.0), 0.0), 1.0);
    return vec3(red, green, blue);
}

void main() {
    gl_Position = projection*view*vec4(Particle[v_index][0].xyz, 1.0); // set vertex position, w=1.0
    // int voxel_id = int(round(Particle[v_index][0].w));
    // vec3 voxel_center = vec3(float(Voxel[(voxel_id-1)*voxel_memory_length+1])*h, float(Voxel[(voxel_id-1)*voxel_memory_length+2])*h, float(Voxel[(voxel_id-1)*voxel_memory_length+3])*h);
    // float l = length(Particle[v_index][3].xyz);
    //v_color = vec4(abs(Particle[v_index][3].xyz), 1.0); // set output color by its acc
    //v_color = vec4(abs(sin(float(voxel_id/2))), abs(cos(float(voxel_id/3))), abs(sin(float(voxel_id/5))), 0.3);
    switch(color_type){
        case 0:  // velocity
            v_color = vec4(abs(Particle[v_index][1].xyz), 1.0);
            break;
        case 1:  // acc
            v_color = vec4(abs(Particle[v_index][3].xyz), 1.0);
            break;
        case 2:  // pressure(density)
            v_color = vec4(get_color_gradient(Particle[v_index][2].w/1000.0, 0.5).xyz, 1.0);
            break;
        case 3:  // 2phase
            if(ParticleSubData[v_index][3].w==1.0){v_color = vec4(0.3, 0.9, 1.0, 1.0);}
            else if(ParticleSubData[v_index][3].w==2.0){v_color = vec4(0.8, 1.0, 0.3, 1.0);}
            break;
    }
}