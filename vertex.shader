#version 460 core

layout(location=0) in int v_index; // vertex id
out vec4 v_color; // color output

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

uniform mat4 projection;
uniform mat4 view;


void main() {
    gl_Position = projection*view*vec4(Particle[v_index][0].xyz, 1.0); // set vertex position, w=1.0
    int voxel_id = int(Particle[v_index][0].w);
    vec3 voxel_center = Voxel[voxel_id][0].yzw;
    v_color = vec4(length(Particle[v_index][3].xyz), length(Particle[v_index][3].xyz), length(Particle[v_index][3].xyz), 0.3); // set output color by its voxel id
    //v_color = vec4(abs(sin(float(voxel_id/2))), abs(cos(float(voxel_id/3))), abs(sin(float(voxel_id/5))), 0.3);
}