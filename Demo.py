import numpy as np

from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader

from SpaceDivision import CreateVoxels, CreateParticles, CreateBoundaryParticles
import time


class Demo:
    def __init__(self):
        self.H = 0.5/2
        self.R = 0.0625/2
        self.Domain = [[0, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [3, 3, 3]]  # 8x3
        t = time.time()
        self.voxels = CreateVoxels(domain=self.Domain, h=self.H)()
        print(self.voxels.shape, "voxel {}s".format(time.time()-t))
        t = time.time()
        self.particles = CreateParticles(domain=[[0, 1, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [2, 2, 1]], h=self.H, r=self.R)()
        print(self.particles.shape, "particle {}s".format(time.time()-t))
        t = time.time()
        self.boundary_particles = CreateBoundaryParticles(domain=self.Domain, h=self.H, r=self.R)()
        print(self.boundary_particles.shape, "boundary particle {}s".format(time.time()-t))

        self.voxel_number = self.voxels.shape[0] // 80  # (n * 80, 4)
        self.particle_number = self.particles.shape[0] // 4  # (n * 4, 4)
        self.boundary_particle_number = self.boundary_particles.shape[0] // 4

        self.voxel_particle_numbers = np.zeros((self.voxel_number, ), dtype=np.int32)

        self.indices_buffer = np.array([i for i in range(max(self.particle_number, self.boundary_particle_number))], dtype=np.int32)

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
        glNamedBufferStorage(self.sbo_boundary_particles, self.boundary_particles.nbytes, self.boundary_particles,
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

        # voxel_particle_in_numbers buffer
        self.sbo_voxel_particle_in_numbers = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxel_particle_in_numbers)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, self.sbo_voxel_particle_in_numbers)
        glNamedBufferStorage(self.sbo_voxel_particle_in_numbers, self.voxel_particle_numbers.nbytes, self.voxel_particle_numbers, GL_DYNAMIC_STORAGE_BIT)

        # voxel_particle_out_numbers buffer
        self.sbo_voxel_particle_out_numbers = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxel_particle_out_numbers)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, self.sbo_voxel_particle_out_numbers)
        glNamedBufferStorage(self.sbo_voxel_particle_out_numbers, self.voxel_particle_numbers.nbytes, self.voxel_particle_numbers, GL_DYNAMIC_STORAGE_BIT)

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
        self.compute_shader_0 = compileProgram(compileShader(open("compute_0_init_boundary_particles.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_0)
        self.compute_shader_0_n_particle_loc = glGetUniformLocation(self.compute_shader_0, "n_particle")
        self.compute_shader_0_n_voxel_loc = glGetUniformLocation(self.compute_shader_0, "n_voxel")
        self.compute_shader_0_h_loc = glGetUniformLocation(self.compute_shader_0, "h")

        glUniform1i(self.compute_shader_0_n_particle_loc, int(self.particle_number))
        glUniform1i(self.compute_shader_0_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.compute_shader_0_h_loc, self.H)

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

        # compute shader 4
        self.compute_shader_4 = compileProgram(
            compileShader(open("compute_4_integrate_solver.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_4)
        self.compute_shader_4_n_particle_loc = glGetUniformLocation(self.compute_shader_4, "n_particle")
        self.compute_shader_4_n_voxel_loc = glGetUniformLocation(self.compute_shader_4, "n_voxel")
        self.compute_shader_4_h_loc = glGetUniformLocation(self.compute_shader_4, "h")

        glUniform1i(self.compute_shader_4_n_particle_loc, int(self.particle_number))
        glUniform1i(self.compute_shader_4_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.compute_shader_4_h_loc, self.H)

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

        # render shader boundary
        self.render_shader_boundary = compileProgram(compileShader(open("boundary_vertex.shader", "rb"), GL_VERTEX_SHADER),
                                                     compileShader(open("fragment.shader", "rb"), GL_FRAGMENT_SHADER))
        glUseProgram(self.render_shader_boundary)
        self.render_shader_boundary_n_particle_loc = glGetUniformLocation(self.render_shader_boundary, "n_particle")
        self.render_shader_boundary_n_voxel_loc = glGetUniformLocation(self.render_shader_boundary, "n_voxel")
        self.render_shader_boundary_h_loc = glGetUniformLocation(self.render_shader_boundary, "h")
        self.boundary_projection_loc = glGetUniformLocation(self.render_shader_boundary, "projection")
        self.boundary_view_loc = glGetUniformLocation(self.render_shader_boundary, "view")
        #
        glUniform1i(self.render_shader_boundary_n_particle_loc, int(self.particle_number))
        glUniform1i(self.render_shader_boundary_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.render_shader_boundary_h_loc, self.H)

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

            glUseProgram(self.compute_shader_0)
            glDispatchCompute(self.boundary_particle_number, 1, 1)
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

            glUseProgram(self.compute_shader_1)
            glDispatchCompute(self.particle_number, 1, 1)
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glUseProgram(self.compute_shader_2)
        glDispatchCompute(self.particle_number, 1, 1)
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glUseProgram(self.compute_shader_3)
        glDispatchCompute(self.particle_number, 1, 1)
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glUseProgram(self.compute_shader_4)
        glDispatchCompute(self.particle_number, 1, 1)
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        # glUseProgram(self.compute_shader_voxel)
        # glUniform1i(self.compute_shader_voxel_id_loc, i)
        # glDispatchCompute(1, 1, 1)
        # glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glBindVertexArray(self.vao)
        glUseProgram(self.render_shader_voxel)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        glDrawArrays(GL_POINTS, 0, self.voxel_number)

        # glUseProgram(self.render_shader_boundary)
        # glPointSize(4)
        # glDrawArrays(GL_POINTS, 0, self.boundary_particle_number)

        # glUseProgram(self.render_shader)
        # glPointSize(15)
        # glDrawArrays(GL_POINTS, 0, self.particle_number)


