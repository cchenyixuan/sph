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


const int voxel_memory_length = 2912;
const int voxel_block_size = 960;



// function definitions

const float PI = 3.141592653589793;
const float REST_DENS = 1000.0;
const float EOS_CONST = 276571;
const float VISC = 10.0;
const float DELTA_T = 0.00045;


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


void ComputeParticleDensityPressure(){
    // position of current particle focused
    vec3 particle_pos = Particle[particle_index-1][0].xyz;
    // delete its density and pressure last time, optional
    Particle[particle_index-1][2].zw = vec2(0.0);
    // neighbourhoods count
    vec2 neighbourhood_counter = vec2(0.0, 0.0);  // .x: domain_particle; .y: boundary_particle
    // voxel_id of current particle
    int voxel_id = int(round(Particle[particle_index-1][0].w));  // starts from 1
    // find neighbourhood vertices, i.e., P_j
    // search in same voxel
    // calculate vertices inside
    for(int j=0; j<voxel_block_size; ++j){
        // vertex index
        int index_j = Voxel[(voxel_id-1)*voxel_memory_length+32+j];  // starts from 1 or -1
        if(index_j==0){continue;}  // empty slot
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
                Particle[particle_index-1][2].z += Particle[index_j-1][1].w * poly6_3d(rij, h);
            }
        }
        // P_j is a boundary particle
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
                Particle[particle_index-1][2].z += BoundaryParticle[index_j-1][1].w * poly6_3d(rij, h);
            }
        }

    }

    // search in neighbourhood voxels
    for(int i=4; i<30; ++i){
        // its neighbourhood voxel
        int neighborhood_id = Voxel[(voxel_id-1)*voxel_memory_length+i];  // starts from 1
        // valid neighborhood
        if(neighborhood_id!=0){
            // calculate vertices inside
            for(int j=0; j<voxel_block_size; ++j){
                // vertex index
                int index_j = Voxel[(neighborhood_id-1)*voxel_memory_length+32+j];  // starts from 1 or -1
                if(index_j==0){continue;}  // empty slot
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
                        Particle[particle_index-1][2].z += Particle[index_j-1][1].w * poly6_3d(rij, h);
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
                        Particle[particle_index-1][2].z += BoundaryParticle[index_j-1][1].w * poly6_3d(rij, h);
                    }
                }

            }
        }
    }
    // compute pressure by EoS
    //   P_i_pressure    = EOS_CONST * ((P_i_rho/REST_DENS)**7 - 1)
    // EOS_CONST = rho0 * (10*v_max)**2 / gamma, where gamma = 7 in this case
    Particle[particle_index-1][2].w = EOS_CONST * (pow(Particle[particle_index-1][2].z/REST_DENS, 7) -1);
    // adapt counter to particle[i][2].xy
    Particle[particle_index-1][2].xy = neighbourhood_counter;

}

void main() {
    ComputeParticleDensityPressure();
}
