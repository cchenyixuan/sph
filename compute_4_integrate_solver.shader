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


float h2 = h * h;

void EulerMethod(){
    // calculate future position
    //   move =             P_velocity           *   dt  /2 +            P_acceleration        *   dt  /2
    vec3 move = Particle[particle_index-1][1].xyz*DELTA_T/2 + Particle[particle_index-1][3].xyz*DELTA_T/2;
    // estimate future position
    //   future_pos =             P_position            + move
    vec3 future_pos = Particle[particle_index-1][0].xyz + move;
    // which voxel this particle would go
    // current vcxel id
    int voxel_id = int(round(Particle[particle_index-1][0].w));  // starts from 1
    // current voxel center
    vec3 voxel_center = vec3(float(Voxel[voxel_id*320+1])*h, float(Voxel[voxel_id*320+2])*h, float(Voxel[voxel_id*320+3])*h);
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


}

void ComputeParticleForce(){
    // position of current particle focused
    vec3 particle_pos = Particle[particle_index-1][0].xyz;
    // voxel_id of current particle
    int voxel_id = int(round(Particle[particle_index-1][0].w));  // starts from 1
    // empty f_pressure, f_viscosity, f_external
    vec3 f_pressure = vec3(0.0);
    vec3 f_viscosity = vec3(0.0);
    vec3 f_external = vec3(0.0, -9.8, 0.0);  // gravity
    // find neighbourhood vertices, i.e., P_j
    // search in same voxel
    // calculate vertices inside
    for(int j=0; j<96; ++j){
        // vertex index
        int index_j = Voxel[(voxel_id-1)*320+32+j];  // starts from 1
        if(index_j==0){continue;}  // empty slot
        if(particle_index==index_j){continue;}
        // P_j is a domain particle
        if(index_j>0){
            // vector xij
            vec3 xij = Particle[index_j-1][0].xyz - particle_pos;
            // distance rij
            float rij = length(xij);
            // distance less than h
            if(rij<h){
                // add f_pressure and f_viscosity
                // f_press +=        grad_spiky_3d(xij, rij, H)          * (           MASS         *(      P_j_pressure      /            P_j_rho**2           +         P_i_pressure           /                 P_i_rho**2            )) *             P_i_rho
                f_pressure += grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (Particle[index_j-1][1].w*(Particle[index_j-1][2].w/pow(Particle[index_j-1][2].z, 2) + Particle[particle_index-1][2].w/pow(Particle[particle_index-1][2].z, 2))) * Particle[particle_index-1][2].z;
                // f_visco  += VISC * (         P_j_mass        /         P_j_rho         ) * (        P_j_velocity       -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                f_viscosity += VISC * (Particle[index_j-1][1].w / Particle[index_j-1][2].w) * (Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
            }
        }
        // P_j is a boundary particle
        else if(index_j<0){
            // reverse index_j
            index_j = abs(index_j);
            // vector xij
            vec3 xij = BoundaryParticle[index_j-1][0].xyz - particle_pos;
            // distance rij
            float rij = length(xij);
            // distance less than h
            if(rij<h){
                // add f_pressure and f_viscosity
                // f_press +=        grad_spiky_3d(xij, rij, H)          * (               MASS             *(          P_j_pressure          /                P_j_rho**2               +         P_i_pressure           /                 P_i_rho**2            )) *             P_i_rho
                f_pressure += grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (BoundaryParticle[index_j-1][1].w*(BoundaryParticle[index_j-1][2].w/pow(BoundaryParticle[index_j-1][2].z, 2) + Particle[particle_index-1][2].w/pow(Particle[particle_index-1][2].z, 2))) * Particle[particle_index-1][2].z;
                // f_visco  += VISC * (             P_j_mass            /             P_j_rho             ) * (            P_j_velocity           -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                f_viscosity += VISC * (BoundaryParticle[index_j-1][1].w / BoundaryParticle[index_j-1][2].w) * (BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
            }
        }

    }

    // search in neighbourhood voxels
    for(int i=4; i<32; ++i){
        // its neighbourhood voxel
        int neighborhood_id = Voxel[(voxel_id-1)*320+i];  // starts from 1
        // valid neighborhood
        if(neighborhood_id!=0){
            // calculate vertices inside
            for(int j=0; j<96; ++j){
                // vertex index
                int index_j = Voxel[(neighborhood_id-1)*320+32+j];  // starts from 1
                if(index_j==0){continue;}  // empty slot
                if(particle_index==index_j){continue;}
                // P_j is a domain particle
                if(index_j>0){
                    // vector xij
                    vec3 xij = Particle[index_j-1][0].xyz - particle_pos;
                    // distance rij
                    float rij = length(xij);
                    // distance less than h
                    if(rij<h){
                        // add f_pressure and f_viscosity
                        // f_press +=        grad_spiky_3d(xij, rij, H)          * (           MASS         *(      P_j_pressure      /            P_j_rho**2           +         P_i_pressure           /                 P_i_rho**2            )) *             P_i_rho
                        f_pressure += grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (Particle[index_j-1][1].w*(Particle[index_j-1][2].w/pow(Particle[index_j-1][2].z, 2) + Particle[particle_index-1][2].w/pow(Particle[particle_index-1][2].z, 2))) * Particle[particle_index-1][2].z;
                        // f_visco  += VISC * (         P_j_mass        /         P_j_rho         ) * (        P_j_velocity       -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                        f_viscosity += VISC * (Particle[index_j-1][1].w / Particle[index_j-1][2].w) * (Particle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                    }
                }
                // P_j is a boundary particle
                else if(index_j<0){
                    // reverse index_j
                    index_j = abs(index_j);
                    // vector xij
                    vec3 xij = BoundaryParticle[index_j-1][0].xyz - particle_pos;
                    // distance rij
                    float rij = length(xij);
                    // distance less than h
                    if(rij<h){
                        // add f_pressure and f_viscosity
                        // f_press +=        grad_spiky_3d(xij, rij, H)          * (               MASS             *(          P_j_pressure          /                P_j_rho**2               +         P_i_pressure           /                 P_i_rho**2            )) *             P_i_rho
                        f_pressure += grad_spiky_3d(xij.x, xij.y, xij.z, rij, h) * (BoundaryParticle[index_j-1][1].w*(BoundaryParticle[index_j-1][2].w/pow(BoundaryParticle[index_j-1][2].z, 2) + Particle[particle_index-1][2].w/pow(Particle[particle_index-1][2].z, 2))) * Particle[particle_index-1][2].z;
                        // f_visco  += VISC * (             P_j_mass            /             P_j_rho             ) * (            P_j_velocity           -            P_i_velocity          ) * lap_viscosity_3d(rij, h)
                        f_viscosity += VISC * (BoundaryParticle[index_j-1][1].w / BoundaryParticle[index_j-1][2].w) * (BoundaryParticle[index_j-1][1].xyz - Particle[particle_index-1][1].xyz) * lap_viscosity_3d(rij, h);
                    }
                }

            }
        }
    }

    // compute force
    //      P_i_force      = f_pressure + f_viscosity + f_external
    Particle[particle_index-1][3].xyz = (f_pressure + f_viscosity + f_external)*0.01/Particle[particle_index-1][2].z;
}

void main() {
    ComputeParticleForce();
}
