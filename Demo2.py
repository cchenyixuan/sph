import numpy as np

from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader

from SpaceDivision import CreateVoxels, CreateParticles, CreateBoundaryParticles, LoadParticleObj
import time
from utils.float_to_fraction import find_fraction

from project_loader import ProjectTwoPhase, Project


class Demo:
    def __init__(self):
        self.H = 0.013
        self.R = 0.0025
        self.DELTA_T = 0.0002
        self.VISCOSITY = 0.001
        self.COHESION = 0.0001
        self.ADHESION = 0.0001
        # H = h_p/h_q and R = r_p/r_q
        self.h_p, self.h_q = find_fraction(np.format_float_positional(self.H))
        self.r_p, self.r_q = find_fraction(np.format_float_positional(self.R))
        self.t_p, self.t_q = find_fraction(np.format_float_positional(self.DELTA_T))
        self.v_p, self.v_q = find_fraction(np.format_float_positional(self.VISCOSITY))
        self.c_p, self.c_q = find_fraction(np.format_float_positional(self.COHESION))
        self.a_p, self.a_q = find_fraction(np.format_float_positional(self.ADHESION))
        """
        self.Domain = [[0, 0, 0], [2, 0.5, 0.5]]  # 8x3
        t = time.time()
        self.voxels = CreateVoxels(domain=self.Domain, h=self.H)()
        print(self.voxels.shape, "voxel {}s".format(time.time()-t))
        t = time.time()
        self.particles = CreateParticles(domain=[[0.8, 0.1, 0.1], [1.8, 0.49, 0.49]], h=self.H, r=self.R)()
        print(self.particles.shape, "particle {}s".format(time.time()-t))
        t = time.time()
        self.boundary_particles = CreateBoundaryParticles(domain=self.Domain, h=self.H, r=self.R)()
        print(self.boundary_particles.shape, "boundary particle {}s".format(time.time()-t))
        """
        # self.project = ProjectTwoPhase(0.02, r"./models/CUBE_FRAME.obj", r"./models/phase1.obj",
        #               r"./models/phase2.obj", r"./models/CUBE_COVER.obj")
        self.project = Project(self.H, self.R, r".\models\toy_model_frame.obj", r".\models\toy_model_domain.obj", r".\models\toy_model_boundary.obj", r"./models/pump_slice_recurrent.obj")
        self.voxels = self.project.voxels
        self.particles = self.project.particles
        self.boundary_particles = self.project.boundary_particles
        # self.tube = LoadParticleObj(r"./models/pumps.obj", 2.0, 0.002*1000)()
        # self.boundary_particles = np.vstack((self.boundary_particles, self.tube))

        self.voxel_number = self.voxels.shape[0] // (182*4)  # (n * (182*4), 4)
        if self.voxel_number > 300000:
            self.voxels_a = self.voxels[:182*4*300000]
            self.voxels_b = self.voxels[182*4*300000:]
        self.particle_number = self.particles.shape[0] // 4  # (n * 4, 4)
        self.boundary_particle_number = self.boundary_particles.shape[0] // 4

        self.voxel_particle_numbers = np.zeros((self.voxel_number, ), dtype=np.int32)

        self.indices_buffer = np.array([i for i in range(max(self.particle_number, self.boundary_particle_number, self.voxel_number))], dtype=np.int32)

        print(self.particle_number, self.boundary_particle_number, self.voxel_number)
        # global status buffer
        # [n_particle, n_boundary_particle, n_voxel, voxel_memory_length, voxel_block_size, h_p, h_q, r_p, r_q, max_velocity_n_times_than_r, rest_dense, eos_constant, t_p, t_q, v_p, v_q, c_p, c_q, a_p, a_q]
        self.global_status = np.array((self.particle_number, self.boundary_particle_number, self.voxel_number, 2912, 960), dtype=np.int32)
        self.global_status_float = np.array((self.H, self.R, self.DELTA_T, self.VISCOSITY, self.COHESION, self.ADHESION), dtype=np.float32)

        # particle_sub_data_buffer
        self.particles_sub_data = self.project.particles_buffer

        # initialize OpenGL
        # particles buffer
        self.sbo_particles = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_particles)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, self.sbo_particles)
        glNamedBufferStorage(self.sbo_particles, self.particles.nbytes, self.particles, GL_DYNAMIC_STORAGE_BIT)

        # particles sub data buffer
        self.sbo_particles_sub_data = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_particles_sub_data)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, self.sbo_particles_sub_data)
        glNamedBufferStorage(self.sbo_particles_sub_data, self.particles_sub_data.nbytes, self.particles_sub_data, GL_DYNAMIC_STORAGE_BIT)

        # boundary buffer
        self.sbo_boundary_particles = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_boundary_particles)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, self.sbo_boundary_particles)
        glNamedBufferStorage(self.sbo_boundary_particles, self.boundary_particles.nbytes, self.boundary_particles,
                             GL_DYNAMIC_STORAGE_BIT)
        # voxels buffer
        self.sbo_voxels = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxels)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, self.sbo_voxels)
        if self.voxel_number <= 300000:
            glNamedBufferStorage(self.sbo_voxels, self.voxels.nbytes, self.voxels, GL_DYNAMIC_STORAGE_BIT)
        else:
            glNamedBufferStorage(self.sbo_voxels, self.voxels_a.nbytes, self.voxels_a, GL_DYNAMIC_STORAGE_BIT)

        # voxels2 buffer
        self.sbo_voxels2 = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxels2)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, self.sbo_voxels2)
        if self.voxel_number <= 300000:
            pass
        else:
            glNamedBufferStorage(self.sbo_voxels2, self.voxels_b.nbytes, self.voxels_b, GL_DYNAMIC_STORAGE_BIT)

        # voxel_particle_numbers buffer
        self.sbo_voxel_particle_numbers = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxel_particle_numbers)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, self.sbo_voxel_particle_numbers)
        glNamedBufferStorage(self.sbo_voxel_particle_numbers, self.voxel_particle_numbers.nbytes, self.voxel_particle_numbers, GL_DYNAMIC_STORAGE_BIT)

        # voxel_particle_in_numbers buffer
        self.sbo_voxel_particle_in_numbers = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxel_particle_in_numbers)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, self.sbo_voxel_particle_in_numbers)
        glNamedBufferStorage(self.sbo_voxel_particle_in_numbers, self.voxel_particle_numbers.nbytes, self.voxel_particle_numbers, GL_DYNAMIC_STORAGE_BIT)

        # voxel_particle_out_numbers buffer
        self.sbo_voxel_particle_out_numbers = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_voxel_particle_out_numbers)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, self.sbo_voxel_particle_out_numbers)
        glNamedBufferStorage(self.sbo_voxel_particle_out_numbers, self.voxel_particle_numbers.nbytes, self.voxel_particle_numbers, GL_DYNAMIC_STORAGE_BIT)

        # global_status buffer
        self.sbo_global_status = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_global_status)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, self.sbo_global_status)
        glNamedBufferStorage(self.sbo_global_status, self.global_status.nbytes, self.global_status, GL_DYNAMIC_STORAGE_BIT)

        # global_status2 buffer
        self.sbo_global_status2 = glGenBuffers(1)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.sbo_global_status2)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, self.sbo_global_status2)
        glNamedBufferStorage(self.sbo_global_status2, self.global_status_float.nbytes, self.global_status_float, GL_DYNAMIC_STORAGE_BIT)

        # vao of indices
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        self.vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, self.indices_buffer.nbytes, self.indices_buffer, GL_STATIC_DRAW)
        glEnableVertexAttribArray(0)
        glVertexAttribIPointer(0, 1, GL_INT, 4, ctypes.c_void_p(0))

        # compute shader
        self.need_init = True

        # compute shader 0
        self.compute_shader_0 = compileProgram(
            compileShader(open("compute_0_init_boundary_particles.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader 1
        self.compute_shader_1 = compileProgram(
            compileShader(open("compute_1_init_domain_particles.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader 2a
        self.compute_shader_2a = compileProgram(
            compileShader(open("compute_2_boundary_density_pressure_solver.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader 2
        self.compute_shader_2 = compileProgram(
            compileShader(open("compute_2_density_pressure_solver.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader 3
        self.compute_shader_3 = compileProgram(
            compileShader(open("compute_3_force_solver.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader 4
        self.compute_shader_4 = compileProgram(
            compileShader(open("compute_4_integrate_solver.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader 5
        self.compute_shader_5 = compileProgram(
            compileShader(open("compute_5_voxel_upgrade_solver.shader", "rb"), GL_COMPUTE_SHADER))
        # compute shader a
        self.compute_shader_a = compileProgram(
            compileShader(open("./MovingBoundaryShaders/compute_a_moving_boundary.shader", "rb"), GL_COMPUTE_SHADER))
        glUseProgram(self.compute_shader_a)
        self.compute_shader_a_current_step_loc = glGetUniformLocation(self.compute_shader_a, "current_step")
        glUniform1i(self.compute_shader_a_current_step_loc, 0)

        # render shader
        self.render_shader = compileProgram(compileShader(open("vertex.shader", "rb"), GL_VERTEX_SHADER),
                                            compileShader(open("fragment.shader", "rb"), GL_FRAGMENT_SHADER))
        glUseProgram(self.render_shader)
        self.projection_loc = glGetUniformLocation(self.render_shader, "projection")
        self.view_loc = glGetUniformLocation(self.render_shader, "view")
        self.render_shader_color_type_loc = glGetUniformLocation(self.render_shader, "color_type")

        glUniform1i(self.render_shader_color_type_loc, 0)

        # render shader vector
        self.render_shader_vector = compileProgram(compileShader(open("VectorShaders/vector_vertex.shader", "rb"), GL_VERTEX_SHADER),
                                                   compileShader(open("VectorShaders/vector_geometry.shader", "rb"), GL_GEOMETRY_SHADER),
                                                   compileShader(open("VectorShaders/vector_fragment.shader", "rb"), GL_FRAGMENT_SHADER))
        glUseProgram(self.render_shader_vector)
        self.render_shader_vector_n_particle_loc = glGetUniformLocation(self.render_shader_vector, "n_particle")
        self.render_shader_vector_n_voxel_loc = glGetUniformLocation(self.render_shader_vector, "n_voxel")
        self.render_shader_vector_h_loc = glGetUniformLocation(self.render_shader_vector, "h")
        self.vector_projection_loc = glGetUniformLocation(self.render_shader_vector, "projection")
        self.vector_view_loc = glGetUniformLocation(self.render_shader_vector, "view")
        self.render_shader_vector_vector_type_loc = glGetUniformLocation(self.render_shader_vector, "vector_type")

        glUniform1i(self.render_shader_vector_n_particle_loc, int(self.particle_number))
        glUniform1i(self.render_shader_vector_n_voxel_loc, int(self.voxel_number))
        glUniform1f(self.render_shader_vector_h_loc, self.H)
        glUniform1i(self.render_shader_vector_vector_type_loc, 0)



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

    def __call__(self, i, pause=False, show_vector=False, show_voxel=False, show_boundary=False):
        if self.need_init:
            self.need_init = False

            glUseProgram(self.compute_shader_0)
            glDispatchCompute(self.boundary_particle_number, 1, 1)
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

            glUseProgram(self.compute_shader_1)
            glDispatchCompute(self.particle_number, 1, 1)
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)
        if not pause:
            glUseProgram(self.compute_shader_2a)
            glDispatchCompute(self.boundary_particle_number, 1, 1)
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

            glUseProgram(self.compute_shader_5)
            glDispatchCompute(self.voxel_number, 1, 1)
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        # glUseProgram(self.compute_shader_voxel)
        # glUniform1i(self.compute_shader_voxel_id_loc, i)
        # glDispatchCompute(1, 1, 1)
        # glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)

        glBindVertexArray(self.vao)

        if show_voxel:
            glUseProgram(self.render_shader_voxel)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glLineWidth(2)
            glDrawArrays(GL_POINTS, 0, self.voxel_number)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

        if show_boundary:
            glUseProgram(self.render_shader_boundary)
            glPointSize(4)
            glDrawArrays(GL_POINTS, 0, self.boundary_particle_number)
            #

        glUseProgram(self.render_shader)
        glPointSize(2)
        glDrawArrays(GL_POINTS, 0, self.particle_number)

        if show_vector:
            glUseProgram(self.render_shader_vector)
            glLineWidth(1)
            glDrawArrays(GL_POINTS, 0, self.particle_number)


