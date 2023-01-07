import numpy as np

from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader

from SpaceDivision import CreateVoxels, CreateParticles


class Demo:
    def __init__(self):
        self.H = 0.1
        self.R = 0.04
        self.Domain = [[0, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [3, 3, 3]]  # 8x3
        self.voxels = CreateVoxels(domain=self.Domain, h=self.H)()
        print(self.voxels.shape)
        self.particles = CreateParticles(domain=self.Domain, h=self.H, r=self.R)()
        print(self.particles.shape)

        self.voxel_number = self.voxels.shape[0] // 80  # (n * 80, 4)
        self.particle_number = self.particles.shape[0] // 4  # (n * 4, 4)

        self.voxel_particle_numbers = np.zeros((self.voxel_number, ), dtype=np.int32)

        self.indices_buffer = np.array([i for i in range(self.particle_number)], dtype=np.int32)

        # initialize OpenGL
        # particles buffer
        self.sbo_particles = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_particles)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, self.sbo_particles)
        glNamedBufferStorage(self.sbo_particles, self.particles.nbytes, self.particles, GL_DYNAMIC_STORAGE_BIT)
        # boundary buffer
        self.sbo_boundary_particles = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_boundary_particles)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, self.sbo_boundary_particles)
        glNamedBufferStorage(self.sbo_boundary_particles, 1024, None,
                             GL_DYNAMIC_STORAGE_BIT)  # TODO: boundary is not defined and will be fixed in future
        # voxels buffer
        self.sbo_voxels = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxels)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, self.sbo_voxels)
        glNamedBufferStorage(self.sbo_voxels, self.voxels.nbytes, self.voxels, GL_DYNAMIC_STORAGE_BIT)

        # voxel_particle_numbers buffer
        self.sbo_voxel_particle_numbers = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxel_particle_numbers)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, self.sbo_voxel_particle_numbers)
        glNamedBufferStorage(self.sbo_voxel_particle_numbers, self.voxel_particle_numbers.nbytes, self.voxel_particle_numbers, GL_DYNAMIC_STORAGE_BIT)

        # vao of indices
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        self.vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, self.indices_buffer.nbytes, self.indices_buffer, GL_STATIC_DRAW)
        glEnableVertexAttribArray(0)
        glVertexAttribIPointer(0, 1, GL_INT, 4, ctypes.c_void_p(0))

        # compute shader
        # compute shader 1
        self.compute_shader_1 = compileProgram(compileShader(open("compute_1_init_domain_particles.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_1)
        self.need_init = True
        self.compute_shader_0_n_particle_loc = glGetUniformLocation(self.compute_shader_1, "n_particle")
        self.compute_shader_0_n_voxel_loc = glGetUniformLocation(self.compute_shader_1, "n_voxel")
        self.compute_shader_0_h_loc = glGetUniformLocation(self.compute_shader_1, "h")

        glUniform1i(self.compute_shader_0_n_particle_loc, int(self.particle_number))
        glUniform1i(self.compute_shader_0_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.compute_shader_0_h_loc, self.H)

        # compute shader 0
        # self.compute_shader_1 = compileProgram(compileShader(open("compute_0_init_boundary_particles.shader", "rb"), GL_COMPUTE_SHADER))
        # glUseProgram(self.compute_shader_1)
        # self.compute_shader_1_n_particle_loc = glGetUniformLocation(self.compute_shader_1, "n_particle")
        # self.compute_shader_1_n_voxel_loc = glGetUniformLocation(self.compute_shader_1, "n_voxel")
        # self.compute_shader_1_h_loc = glGetUniformLocation(self.compute_shader_1, "h")
#
        # glUniform1i(self.compute_shader_1_n_particle_loc, int(self.particle_number))
        # glUniform1i(self.compute_shader_1_n_voxel_loc, int(self.voxel_number))
        # glUniform1f(self.compute_shader_1_h_loc, self.H)

        # compute shader 2
        self.compute_shader_2 = compileProgram(compileShader(open("compute_2_density_pressure_solver.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_2)
        self.compute_shader_2_n_particle_loc = glGetUniformLocation(self.compute_shader_2, "n_particle")
        self.compute_shader_2_n_voxel_loc = glGetUniformLocation(self.compute_shader_2, "n_voxel")
        self.compute_shader_2_h_loc = glGetUniformLocation(self.compute_shader_2, "h")

        glUniform1i(self.compute_shader_2_n_particle_loc, int(self.particle_number))
        glUniform1i(self.compute_shader_2_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.compute_shader_2_h_loc, self.H)

        # compute shader 3
        self.compute_shader_3 = compileProgram(compileShader(open("compute_3_force_solver.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_3)
        self.compute_shader_3_n_particle_loc = glGetUniformLocation(self.compute_shader_3, "n_particle")
        self.compute_shader_3_n_voxel_loc = glGetUniformLocation(self.compute_shader_3, "n_voxel")
        self.compute_shader_3_h_loc = glGetUniformLocation(self.compute_shader_3, "h")

        glUniform1i(self.compute_shader_3_n_particle_loc, int(self.particle_number))
        glUniform1i(self.compute_shader_3_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.compute_shader_3_h_loc, self.H)

        # render shader
        self.render_shader = compileProgram(compileShader(open("vertex.shader", "rb"), GL_VERTEX_SHADER),
                                            compileShader(open("fragment.shader", "rb"), GL_FRAGMENT_SHADER))
        glUseProgram(self.render_shader)
        self.render_shader_n_particle_loc = glGetUniformLocation(self.render_shader, "n_particle")
        self.render_shader_n_voxel_loc = glGetUniformLocation(self.render_shader, "n_voxel")
        self.render_shader_h_loc = glGetUniformLocation(self.render_shader, "h")
        self.projection_loc = glGetUniformLocation(self.render_shader, "projection")
        self.view_loc = glGetUniformLocation(self.render_shader, "view")
#
        glUniform1i(self.render_shader_n_particle_loc, int(self.particle_number))
        glUniform1i(self.render_shader_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.render_shader_h_loc, self.H)

        # glUniformMatrix4fv(self.projection_loc, 1, GL_FALSE, self.camera.projection)
        # glUniformMatrix4fv(self.view_loc, 1, GL_FALSE, self.camera.view)

        # compute shader for voxel debug
        self.compute_shader_voxel = compileProgram(
            compileShader(open("voxel_compute.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_voxel)
        self.compute_shader_voxel_n_particle_loc = glGetUniformLocation(self.compute_shader_voxel, "n_particle")
        self.compute_shader_voxel_n_voxel_loc = glGetUniformLocation(self.compute_shader_voxel, "n_voxel")
        self.compute_shader_voxel_h_loc = glGetUniformLocation(self.compute_shader_voxel, "h")
        self.compute_shader_voxel_id_loc = glGetUniformLocation(self.compute_shader_voxel, "id")

        glUniform1i(self.compute_shader_voxel_n_particle_loc, int(self.particle_number))
        glUniform1i(self.compute_shader_voxel_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.compute_shader_voxel_h_loc, self.H)
        glUniform1i(self.compute_shader_voxel_id_loc, 0)
        # render shader for voxel
        self.render_shader_voxel = compileProgram(compileShader(open("voxel_vertex.shader", "rb"), GL_VERTEX_SHADER),
                                                  compileShader(open("voxel_geometry.shader", "rb"), GL_GEOMETRY_SHADER),
                                                  compileShader(open("voxel_fragment.shader", "rb"), GL_FRAGMENT_SHADER))
        glUseProgram(self.render_shader_voxel)

        self.render_shader_voxel_n_particle_loc = glGetUniformLocation(self.render_shader_voxel, "n_particle")
        self.render_shader_voxel_n_voxel_loc = glGetUniformLocation(self.render_shader_voxel, "n_voxel")
        self.render_shader_voxel_h_loc = glGetUniformLocation(self.compute_shader_1, "h")

        glUniform1i(self.render_shader_voxel_n_particle_loc, int(self.particle_number))
        glUniform1i(self.render_shader_voxel_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.render_shader_voxel_h_loc, self.H)

        self.voxel_projection_loc = glGetUniformLocation(self.render_shader_voxel, "projection")
        self.voxel_view_loc = glGetUniformLocation(self.render_shader_voxel, "view")

    def __call__(self, i):
        if self.need_init:
            self.need_init = False
            glUseProgram(self.compute_shader_1)
            glDispatchCompute(self.particle_number, 1, 1)
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glUseProgram(self.compute_shader_2)
        glDispatchCompute(self.particle_number, 1, 1)
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glUseProgram(self.compute_shader_3)
        glDispatchCompute(self.particle_number, 1, 1)
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        # glUseProgram(self.compute_shader_voxel)
        # glUniform1i(self.compute_shader_voxel_id_loc, i)
        # glDispatchCompute(1, 1, 1)
        # glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glBindVertexArray(self.vao)
        # glUseProgram(self.render_shader_voxel)
        # glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        # glDrawArrays(GL_POINTS, 0, self.voxel_number)

        glUseProgram(self.render_shader)
        glPointSize(15)
        glDrawArrays(GL_POINTS, 0, self.particle_number)
