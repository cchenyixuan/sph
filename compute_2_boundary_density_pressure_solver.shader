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



// from here
void ComputeBoundaryParticleDensityPressure(){
    // position of current particle focused
    vec3 particle_pos = BoundaryParticle[particle_index-1][0].xyz;
    // delete its density and pressure last time, optional
    BoundaryParticle[particle_index-1][2].w = 0.0;
    BoundaryParticle[particle_index-1][3].w = 0.0;
    BoundaryParticle[particle_index-1][3].xyz = vec3(0.0, 0.0, 0.0);
    // WSS direction
    vec3 wall_shear_stress = vec3(0.0, 0.0, 0.0);
    // neighbourhoods count
    vec2 neighbourhood_counter = vec2(0.0, 0.0);// .x: domain_particle; .y: boundary_particle
    // voxel_id of current particle
    int voxel_id = int(round(BoundaryParticle[particle_index-1][0].w));// starts from 1
    // find neighbourhood vertices, i.e., P_j
    // search in same voxel
    // calculate vertices inside
    if(voxel_id<=300000){
        for (int j=0; j<voxel_block_size; ++j){
            // vertex index
            int index_j = Voxel[(voxel_id-1)*voxel_memory_length+32+j];// starts from 1 or -1

            if (index_j==0){ break; }// empty slot
            // P_j is a domain particle
            else if (index_j>0){

                // distance rij
                float rij = distance(particle_pos, Particle[index_j-1][0].xyz);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.x += 1.0;
                    // add density to location (2, 2) of its mat4x4
                    //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                    BoundaryParticle[particle_index-1][2].w += Particle[index_j-1][1].w * poly6_3d(rij, h);
                    wall_shear_stress += (particle_pos - Particle[index_j-1][0].xyz) / (rij + 0.01 * h2) * poly6_3d(rij, h);
                }
            }
            // P_j is a boundary particle
            else if (index_j<0){
                // reverse index_j
                index_j = -index_j;
                // distance rij
                float rij = distance(particle_pos, BoundaryParticle[index_j-1][0].xyz);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.y += 1.0;
                    // add density to location (2, 2) of its mat4x4
                    //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                    BoundaryParticle[particle_index-1][2].w += BoundaryParticle[index_j-1][1].w * poly6_3d(rij, h);
                }
            }
        }
    }
    else{
        for (int j=0; j<voxel_block_size; ++j){
            // vertex index
            int index_j = Voxel2[(voxel_id-300000-1)*voxel_memory_length+32+j];// starts from 1 or -1
            if (index_j==0){ break; }// empty slot
            // P_j is a domain particle
            else if (index_j>0){

                // distance rij
                float rij = distance(particle_pos, Particle[index_j-1][0].xyz);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.x += 1.0;
                    // add density to location (2, 2) of its mat4x4
                    //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                    BoundaryParticle[particle_index-1][2].w += Particle[index_j-1][1].w * poly6_3d(rij, h);
                    wall_shear_stress += (particle_pos - Particle[index_j-1][0].xyz) / (rij + 0.01 * h2) * poly6_3d(rij, h);
                }
            }
            // P_j is a boundary particle
            else if (index_j<0){
                // reverse index_j
                index_j = -index_j;
                // distance rij
                float rij = distance(particle_pos, BoundaryParticle[index_j-1][0].xyz);
                // distance less than h
                if (rij<h){
                    // counter
                    neighbourhood_counter.y += 1.0;
                    // add density to location (2, 2) of its mat4x4
                    //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                    BoundaryParticle[particle_index-1][2].w += BoundaryParticle[index_j-1][1].w * poly6_3d(rij, h);
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
                for(int j=0; j<voxel_block_size; ++j){
                    // vertex index
                    int index_j = Voxel[(neighborhood_id-1)*voxel_memory_length+32+j];// starts from 1 or -1
                    if(index_j==0){break;}  // empty slot
                    // P_j is a domain particle
                    else if(index_j>0){
                        // distance rij
                        float rij = distance(particle_pos, Particle[index_j-1][0].xyz);
                        // distance less than h
                        if(rij<h){
                            // counter
                            neighbourhood_counter.x += 1.0;
                            // add density to location (2, 2) of its mat4x4
                            //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                            BoundaryParticle[particle_index-1][2].w += Particle[index_j-1][1].w * poly6_3d(rij, h);
                            wall_shear_stress += (particle_pos - Particle[index_j-1][0].xyz) / (rij + 0.01 * h2) * poly6_3d(rij, h);
                        }
                    }
                    else if(index_j<0){
                        // reverse index_j
                        index_j = -index_j;
                        // distance rij
                        float rij = distance(particle_pos, BoundaryParticle[index_j-1][0].xyz);
                        // distance less than h
                        if(rij<h){
                            // counter
                            neighbourhood_counter.y += 1.0;
                            // add density to location (2, 2) of its mat4x4
                            //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                            BoundaryParticle[particle_index-1][2].w += BoundaryParticle[index_j-1][1].w * poly6_3d(rij, h);
                        }
                    }
                }
            }
            else{
                for(int j=0; j<voxel_block_size; ++j){
                    // vertex index
                    int index_j = Voxel2[(neighborhood_id-300000-1)*voxel_memory_length+32+j];// starts from 1 or -1
                    if(index_j==0){break;}  // empty slot
                    // P_j is a domain particle
                    else if(index_j>0){
                        // distance rij
                        float rij = distance(particle_pos, Particle[index_j-1][0].xyz);
                        // distance less than h
                        if(rij<h){
                            // counter
                            neighbourhood_counter.x += 1.0;
                            // add density to location (2, 2) of its mat4x4
                            //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                            BoundaryParticle[particle_index-1][2].w += Particle[index_j-1][1].w * poly6_3d(rij, h);
                            wall_shear_stress += (particle_pos - Particle[index_j-1][0].xyz) / (rij + 0.01 * h2) * poly6_3d(rij, h);
                        }
                    }
                    else if(index_j<0){
                        // reverse index_j
                        index_j = -index_j;
                        // distance rij
                        float rij = distance(particle_pos, BoundaryParticle[index_j-1][0].xyz);
                        // distance less than h
                        if(rij<h){
                            // counter
                            neighbourhood_counter.y += 1.0;
                            // add density to location (2, 2) of its mat4x4
                            //     P_i_rho       +=         P_j_mass       * poly6_3d(rij, h)
                            BoundaryParticle[particle_index-1][2].w += BoundaryParticle[index_j-1][1].w * poly6_3d(rij, h);
                        }
                    }
                }
            }
        }
    }
    BoundaryParticle[particle_index-1][3].w = max(eos_constant * (pow(BoundaryParticle[particle_index-1][2].w/rest_dense, 7) -1), 0.0);
    BoundaryParticle[particle_index-1][3].xyz = wall_shear_stress;
}

void main(){
    ComputeBoundaryParticleDensityPressure();
}