#version 460 compatibility

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



layout(local_size_x=1, local_size_y=1, local_size_z=1) in;

uint gid = gl_GlobalInvocationID.x;
int particle_index = int(gid)+1;
float particle_index_float = float(particle_index);
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


float h2 = h * h;

// coefficients
struct Coefficient{
    float Poly6_2d;
    float Poly6_3d;
    float Spiky_2d;
    float Spiky_3d;
    float Viscosity_2d;
    float Viscosity_3d;
};

Coefficient coeff = Coefficient(
    4 / (PI * pow(h, 8)),
    315 / (64 * PI * pow(h, 9)),
    10 / (PI * pow(h, 5)),
    15 / (PI * pow(h, 6)),
    40 / (PI * h2),
    15 / (2 * PI * pow(h, 3))
);


// poly6
float poly6_2d(float rij, float h){
    return max(0.0, coeff.Poly6_2d * pow((h2 - rij * rij),3));
}
float poly6_3d(float rij, float h){
    return max(0.0, coeff.Poly6_3d * pow((h2 - rij * rij),3));
}
vec2 grad_poly6_2d(float x, float y, float rij, float h){
    if (rij > h){return vec2(0.0, 0.0);}
    float w_prime = - 6 * coeff.Poly6_2d * pow((h2 - rij * rij),2);
    return vec2(w_prime * x, w_prime * y);
}
vec3 grad_poly6_3d(float x, float y, float z, float rij, float h){
    if (rij > h){return vec3(0.0, 0.0, 0.0);}
    float w_prime = - 6 * coeff.Poly6_3d * pow((h2 - rij * rij),2);
    return vec3(w_prime * x, w_prime * y, w_prime * z);
}
float lap_poly6_2d(float rij, float h){
    if (rij > h){return 0;}
    return - 12 * coeff.Poly6_2d * (h2 - rij * rij) * (h2 - 3 * rij * rij);
}
float lap_poly6_3d(float rij, float h){
    if (rij > h){return 0;}
    return - 6 * coeff.Poly6_3d * (h2 - rij * rij) * (3 * h2 - 7 * rij * rij);
}

// spiky
float spiky_2d(float rij, float h){
    return max(0.0, coeff.Spiky_2d * pow((h - rij),3));
}
float spiky_3d(float rij, float h){
    return max(0.0, coeff.Spiky_3d * pow((h - rij),3));
}
vec2 grad_spiky_2d(float x, float y, float rij, float h){
    if (rij > h){return vec2(0.0, 0.0);}
    float w_prime = - 3 * coeff.Spiky_2d * pow((h - rij),2);
    return vec2(w_prime * x / rij, w_prime * y / rij);
}
vec3 grad_spiky_3d(float x, float y, float z, float rij, float h){
    if (rij > h){return vec3(0.0, 0.0, 0.0);}
    float w_prime = - 3 * coeff.Spiky_3d * pow((h - rij),2);
    return vec3(w_prime * x / rij, w_prime * y / rij, w_prime * z / rij);
}
float lap_spiky_2d(float rij, float h){
    if (rij > h){return 0;}
    return coeff.Spiky_2d * (- 3 * h2 / rij + 12 * h - 9 * rij);
}
float lap_spiky_3d(float rij, float h){
    if (rij > h){return 0;}
    return coeff.Spiky_3d * (- 6 * h2 / rij + 18 * h - 12 * rij);
}

// viscosity
float viscosity_2d(float rij, float h){
    return max(0.0, coeff.Viscosity_2d * (- rij * rij * rij / (2 * h2) + rij * rij / h2 + h / (2 * rij) -1));
}
float viscosity_3d(float rij, float h){
    return max(0.0, coeff.Viscosity_3d * (- rij * rij * rij / (9 * h2) + rij * rij / (4 * h2) + log(rij / h) / 6 - 5 / 36));
}
vec2 grad_viscosity_2d(float x, float y, float rij, float h){
    if (rij > h){return vec2(0.0, 0.0);}
    float w_prime = coeff.Viscosity_2d * (- rij * rij / (3 * h2 * h) + rij / (2 * h2) - 1/ (6 * rij));
    return vec2(w_prime * x / rij, w_prime * y / rij);
}
vec3 grad_viscosity_3d(float x, float y, float z, float rij, float h){
    if (rij > h){return vec3(0.0, 0.0, 0.0);}
    float w_prime = coeff.Viscosity_3d * (- 3 * rij * rij / (h2 * h) + 2 * rij / h2 - h / (2 * rij * rij));
    return vec3(w_prime * x / rij, w_prime * y / rij, w_prime * z / rij);
}
float lap_viscosity_2d(float rij, float h){
    if (rij > h){return 0;}
    return 6 * coeff.Viscosity_2d / (h * h2) * (h - rij);
}
float lap_viscosity_3d(float rij, float h){
    if (rij > h){return 0;}
    return coeff.Viscosity_3d / (h * h2) * (h - rij);
}

vec3 GetExternalForce(vec3 pos){
    // vec3 gravity = vec3(0.0, -9.81, 0.0);
    vec3 force = vec3(0.0, 0.0, 0.0);
    if(pos.y<0.8){
        vec3 center = vec3(0.61317, 0.67125, pos.z);

        vec3 po = center-pos;
        vec3 right = vec3(0.0, 0.0, 1.0);
        force = normalize(cross(right, po))*2.0;
    }
    else{
        force = vec3(0.0, 0.0, 0.0);
    }



    return force;
}



void ComputeParticleForce(){
    // position of current particle focused
    vec3 particle_pos = Particle[particle_index-1][0].xyz;
    // voxel_id of current particle
    int voxel_id = int(round(Particle[particle_index-1][0].w));  // starts from 1
    // empty f_pressure, f_viscosity, f_external
    vec3 f_pressure = vec3(0.0, 0.0, 0.0);
    vec3 f_viscosity = vec3(0.0, 0.0, 0.0);
    vec3 f_external = vec3(0.0, -9.81, 0.0);//GetExternalForce(particle_pos);  // gravity
    vec3 f_cohesion = vec3(0.0, 0.0, 0.0);  // surface tension of domain particles
    vec3 f_adhesion = vec3(0.0, 0.0, 0.0);  // surface tension of boundary particles
    // vec3 f_transfer = vec3(0.0, 0.0, 0.0);  // Vorticity transfer force
    // vec3 t_transfer = vec3(0.0, 0.0, 0.0);  // Vorticity transfer torque
    // neighbourhoods count
    vec2 neighbourhood_counter = vec2(0.0, 0.0);  // .x: domain_particle; .y: boundary_particle
    // find neighbourhood vertices, i.e., P_j
    // search in same voxel
    // calculate vertices inside
    if(voxel_id<=300000){
        for (int j=0; j<voxel_block_size; ++j){
            // vertex index
            int index_j = Voxel[(voxel_id-1)*voxel_memory_length+32+j];// starts from 1
            if (index_j==0){ break; }// empty slot
            if (particle_index==index_j){ continue; }
            // P_j is a domain particle
            if (index_j>0){
                // vector xij
                vec3 xij = particle_pos - Particle[index_j-1][0].xyz;
                // distance rij
                float rij = length(xij);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.x += 1.0;
                    // add f_pressure and f_viscosity
                    // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                    f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (Particle[index_j-1][1].w*(Particle[index_j-1][3].w/pow(Particle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                    // f_visco  += VISC * (         P_j_mass        /         P_j_rho         ) * (        P_j_velocity       -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                    // f_viscosity += VISC * (Particle[index_j-1][1].w / Particle[index_j-1][2].w) * (Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                    //f_viscosity += Particle[particle_index-1][1].w*VISC*(Particle[index_j-1][1].w/Particle[index_j-1][2].z)*(Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                    // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                    f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (Particle[index_j-1][1].w/Particle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-Particle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                    // f_cohesion -= COHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                    f_cohesion -= cohesion * Particle[index_j-1][1].w*xij*poly6_3d(rij, h);

                    // f_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v_angluar-P_j_v_angluar)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // VISC_TRANSFER refers to mu/rho
                    // f_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][2].xyz-Particle[index_j-1][2].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));
                    // t_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v-P_j_v)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // this term needs to -2*P_i_v_angluar afrterwards
                    // t_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));

                }
            }
            // P_j is a boundary particle
            else if (index_j<0){
                // reverse index_j
                index_j = abs(index_j);
                // vector xij
                vec3 xij = particle_pos - BoundaryParticle[index_j-1][0].xyz;
                // distance rij
                float rij = length(xij);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.y += 1.0;
                    // add f_pressure and f_viscosity
                    // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                    f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (BoundaryParticle[index_j-1][1].w*(BoundaryParticle[index_j-1][3].w/pow(BoundaryParticle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                    // f_visco  += VISC * (             P_j_mass            /             P_j_rho             ) * (            P_j_velocity           -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                    // f_viscosity += VISC * (BoundaryParticle[index_j-1][1].w / BoundaryParticle[index_j-1][2].w) * (BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                    //f_viscosity += Particle[particle_index-1][1].w*VISC*(BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].z)*(BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                    // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                    f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-BoundaryParticle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-BoundaryParticle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                    // f_adhesion -= ADHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                    f_adhesion -= adhesion * BoundaryParticle[index_j-1][1].w*xij*poly6_3d(rij, h);
                }
            }

        }
    }
    else{
        for (int j=0; j<voxel_block_size; ++j){
            // vertex index
            int index_j = Voxel2[(voxel_id-300000-1)*voxel_memory_length+32+j];// starts from 1
            if (index_j==0){ break; }// empty slot
            if (particle_index==index_j){ continue; }
            // P_j is a domain particle
            if (index_j>0){
                // vector xij
                vec3 xij = particle_pos - Particle[index_j-1][0].xyz;
                // distance rij
                float rij = length(xij);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.x += 1.0;
                    // add f_pressure and f_viscosity
                    // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                    f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (Particle[index_j-1][1].w*(Particle[index_j-1][3].w/pow(Particle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                    // f_visco  += VISC * (         P_j_mass        /         P_j_rho         ) * (        P_j_velocity       -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                    // f_viscosity += VISC * (Particle[index_j-1][1].w / Particle[index_j-1][2].w) * (Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                    //f_viscosity += Particle[particle_index-1][1].w*VISC*(Particle[index_j-1][1].w/Particle[index_j-1][2].z)*(Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                    // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                    f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (Particle[index_j-1][1].w/Particle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-Particle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                    // f_cohesion -= COHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                    f_cohesion -= cohesion * Particle[index_j-1][1].w*xij*poly6_3d(rij, h);

                    // f_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v_angluar-P_j_v_angluar)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // VISC_TRANSFER refers to mu/rho
                    // f_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][2].xyz-Particle[index_j-1][2].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));
                    // t_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v-P_j_v)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // this term needs to -2*P_i_v_angluar afrterwards
                    // t_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));

                }
            }
            // P_j is a boundary particle
            else if (index_j<0){
                // reverse index_j
                index_j = abs(index_j);
                // vector xij
                vec3 xij = particle_pos - BoundaryParticle[index_j-1][0].xyz;
                // distance rij
                float rij = length(xij);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.y += 1.0;
                    // add f_pressure and f_viscosity
                    // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                    f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (BoundaryParticle[index_j-1][1].w*(BoundaryParticle[index_j-1][3].w/pow(BoundaryParticle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                    // f_visco  += VISC * (             P_j_mass            /             P_j_rho             ) * (            P_j_velocity           -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                    // f_viscosity += VISC * (BoundaryParticle[index_j-1][1].w / BoundaryParticle[index_j-1][2].w) * (BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                    //f_viscosity += Particle[particle_index-1][1].w*VISC*(BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].z)*(BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                    // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                    f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-BoundaryParticle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-BoundaryParticle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                    // f_adhesion -= ADHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                    f_adhesion -= adhesion * BoundaryParticle[index_j-1][1].w*xij*poly6_3d(rij, h);
                }
            }

        }
    }

    // search in neighbourhood voxels
    for(int i=4; i<30; ++i){
        // its neighbourhood voxel
        int neighborhood_id;
        if(voxel_id<=300000){
            neighborhood_id = Voxel[(voxel_id-1)*voxel_memory_length+i];  // starts from 1
        }
        else{
            neighborhood_id = Voxel2[(voxel_id-300000-1)*voxel_memory_length+i];  // starts from 1
        }
        // valid neighborhood
        if(neighborhood_id!=0){
            // calculate vertices inside
            if(neighborhood_id<=300000){
                for (int j=0; j<voxel_block_size; ++j){
                    // vertex index
                    int index_j = Voxel[(neighborhood_id-1)*voxel_memory_length+32+j];// starts from 1
                    if (index_j==0){ break; }// empty slot
                    if (particle_index==index_j){ continue; }
                    // P_j is a domain particle
                    if (index_j>0){
                        // vector xij
                        vec3 xij = particle_pos - Particle[index_j-1][0].xyz;
                        // distance rij
                        float rij = length(xij);
                        // distance less than h
                        if (rij<h){
                            // counter
                            neighbourhood_counter.x += 1.0;
                            // add f_pressure and f_viscosity
                            // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                            f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (Particle[index_j-1][1].w*(Particle[index_j-1][3].w/pow(Particle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                            // f_visco  += VISC * (         P_j_mass        /         P_j_rho         ) * (        P_j_velocity       -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                            // f_viscosity += VISC * (Particle[index_j-1][1].w / Particle[index_j-1][2].w) * (Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                            //f_viscosity += Particle[particle_index-1][1].w*VISC*(Particle[index_j-1][1].w/Particle[index_j-1][2].z)*(Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                            // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                            f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (Particle[index_j-1][1].w/Particle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-Particle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                            // f_cohesion -= COHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                            f_cohesion -= cohesion * Particle[index_j-1][1].w*xij*poly6_3d(rij, h);

                            // f_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v_angluar-P_j_v_angluar)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // VISC_TRANSFER refers to mu/rho
                            // f_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][2].xyz-Particle[index_j-1][2].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));
                            // t_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v-P_j_v)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // this term needs to -2*P_i_v_angluar afrterwards
                            // t_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));

                        }
                    }
                    // P_j is a boundary particle
                    else if (index_j<0){
                        // reverse index_j
                        index_j = abs(index_j);
                        // vector xij
                        vec3 xij = particle_pos - BoundaryParticle[index_j-1][0].xyz;
                        // distance rij
                        float rij = length(xij);
                        // distance less than h
                        if (rij<h){
                            // counter
                            neighbourhood_counter.y += 1.0;
                            // add f_pressure and f_viscosity
                            // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                            f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (BoundaryParticle[index_j-1][1].w*(BoundaryParticle[index_j-1][3].w/pow(BoundaryParticle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                            // f_visco  += VISC * (             P_j_mass            /             P_j_rho             ) * (            P_j_velocity           -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                            // f_viscosity += VISC * (BoundaryParticle[index_j-1][1].w / BoundaryParticle[index_j-1][2].w) * (BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                            //f_viscosity += Particle[particle_index-1][1].w*VISC*(BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].z)*(BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                            // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                            f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-BoundaryParticle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-BoundaryParticle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                            // f_adhesion -= ADHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                            f_adhesion -= adhesion * BoundaryParticle[index_j-1][1].w*xij*poly6_3d(rij, h);

                            // f_visco = m_i*mu*lap(v_i)
                            //lap_v_i  = sum{m_i*mu * m_j/rho_j * (vi-vj) * 2||(lap_W(ij))||/||rij||}
                            // f_viscosity += m_i*VISC*(m_j/rho_j)*(vj-vi)*(2*lap_viscosity_3d(rij, h)/rij);
                            // f_viscosity += Particle[particle_index-1][1].w*VISC*(Particle[index_j-1][1].w/Particle[index_j-1][2].z)*(Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*lap_viscosity_3d(rij, h)/rij);
                        }
                    }

                }
            }
            else{
                for (int j=0; j<voxel_block_size; ++j){
                    // vertex index
                    int index_j = Voxel2[(neighborhood_id-300000-1)*voxel_memory_length+32+j];// starts from 1
                    if (index_j==0){ break; }// empty slot
                    if (particle_index==index_j){ continue; }
                    // P_j is a domain particle
                    if (index_j>0){
                        // vector xij
                        vec3 xij = particle_pos - Particle[index_j-1][0].xyz;
                        // distance rij
                        float rij = length(xij);
                        // distance less than h
                        if (rij<h){
                            // counter
                            neighbourhood_counter.x += 1.0;
                            // add f_pressure and f_viscosity
                            // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                            f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (Particle[index_j-1][1].w*(Particle[index_j-1][3].w/pow(Particle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                            // f_visco  += VISC * (         P_j_mass        /         P_j_rho         ) * (        P_j_velocity       -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                            // f_viscosity += VISC * (Particle[index_j-1][1].w / Particle[index_j-1][2].w) * (Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                            //f_viscosity += Particle[particle_index-1][1].w*VISC*(Particle[index_j-1][1].w/Particle[index_j-1][2].z)*(Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                            // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                            f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (Particle[index_j-1][1].w/Particle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-Particle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                            // f_cohesion -= COHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                            f_cohesion -= cohesion * Particle[index_j-1][1].w*xij*poly6_3d(rij, h);

                            // f_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v_angluar-P_j_v_angluar)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // VISC_TRANSFER refers to mu/rho
                            // f_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][2].xyz-Particle[index_j-1][2].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));
                            // t_transfer += MASS_i*VISC_TRANSFER * 1/rho_i * (P_i_v-P_j_v)xgrad_spiky_3d(xij.x, xij.y, xij.z, rij, h);  // this term needs to -2*P_i_v_angluar afrterwards
                            // t_transfer += Particle[particle_index-1][1].w*VISC_TRANSFER/Particle[particle_index-1][2].w * cross((Particle[particle_index-1][1].xyz-Particle[index_j-1][1].xyz), grad_spiky_3d(xij.x, xij.y, xij.z, rij, h));

                        }
                    }
                    // P_j is a boundary particle
                    else if (index_j<0){
                        // reverse index_j
                        index_j = abs(index_j);
                        // vector xij
                        vec3 xij = particle_pos - BoundaryParticle[index_j-1][0].xyz;
                        // distance rij
                        float rij = length(xij);
                        // distance less than h
                        if (rij<h){
                            // counter
                            neighbourhood_counter.y += 1.0;
                            // add f_pressure and f_viscosity
                            // f_press -= MASS_i*grad_spiky_3d(xij, rij, H)*(MASS_j*(P_j_pressure/P_j_rho**2 + P_i_pressure/P_i_rho**2))
                            f_pressure -= Particle[particle_index-1][1].w*grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (BoundaryParticle[index_j-1][1].w*(BoundaryParticle[index_j-1][3].w/pow(BoundaryParticle[index_j-1][2].w, 2) + Particle[particle_index-1][3].w/pow(Particle[particle_index-1][2].w, 2)));
                            // f_visco  += VISC * (             P_j_mass            /             P_j_rho             ) * (            P_j_velocity           -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                            // f_viscosity += VISC * (BoundaryParticle[index_j-1][1].w / BoundaryParticle[index_j-1][2].w) * (BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                            //f_viscosity += Particle[particle_index-1][1].w*VISC*(BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].z)*(BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*length(grad_spiky_3d(xij.x, xij.y, xij.z, rij, h))/rij);
                            // f_visco  += MASS_i * VISC * 2*(dimension+2) * MASS_j/P_j_rho * (P_i_v - P_j_v)*(P_i_x-P_j_x)/(rij*rij + 0.01*H**2) * grad_spiky_3d(xij, rij, H)
                            f_viscosity += Particle[particle_index-1][1].w*viscosity* 10 * (BoundaryParticle[index_j-1][1].w/BoundaryParticle[index_j-1][2].w) * dot(Particle[particle_index-1][1].xyz-BoundaryParticle[index_j-1][1].xyz, Particle[particle_index-1][0].xyz-BoundaryParticle[index_j-1][0].xyz)/(rij*rij+0.01*h2) * grad_spiky_3d(xij.x, xij.y, xij.z, rij, h);
                            // f_adhesion -= ADHESION * MASS_j*(P_i_x-P_j_x)*poly6_3d(rij, h);
                            f_adhesion -= adhesion * BoundaryParticle[index_j-1][1].w*xij*poly6_3d(rij, h);

                            // f_visco = m_i*mu*lap(v_i)
                            //lap_v_i  = sum{m_i*mu * m_j/rho_j * (vi-vj) * 2||(lap_W(ij))||/||rij||}
                            // f_viscosity += m_i*VISC*(m_j/rho_j)*(vj-vi)*(2*lap_viscosity_3d(rij, h)/rij);
                            // f_viscosity += Particle[particle_index-1][1].w*VISC*(Particle[index_j-1][1].w/Particle[index_j-1][2].z)*(Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz)*(2*lap_viscosity_3d(rij, h)/rij);
                        }
                    }

                }
            }
        }
    }

    // compute force
    //            P_i_acceleration           = (f_pressure + f_viscosity + f_external)/mass
    // t_transfer -= Particle[particle_index-1][1].w*VISC_TRANSFER*2*Particle[particle_index-1][2].xyz;
    // ParticleSubData[particle_index-1][0].xyz = t_transfer;
    Particle[particle_index-1][3].xyz = (f_pressure + f_viscosity + f_external*Particle[particle_index-1][1].w + f_cohesion + f_adhesion)/Particle[particle_index-1][1].w;
    // adapt counter to particle[i][2].xy
    //Particle[particle_index-1][2].xy = neighbourhood_counter;
}

void main() {
    ComputeParticleForce();
}
