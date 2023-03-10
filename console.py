import time
import threading

import pyrr
import numpy as np
from OpenGL.GL import *
import glfw
import Demo
import Demo2
from camera import Camera
from PIL import Image
from Coordinates import Coord


console = False
console_buffer = """>>>"""


class DisplayPort:
    def __init__(self):
        glfw.init()
        self.window = glfw.create_window(1920, 1080, "Console", None, None)
        glfw.set_window_pos(self.window, 0, 30)
        glfw.hide_window(self.window)
        print("DisplayPort Initialized.")
        self.cursor_position = (0.0, 0.0)
        self.offset = 0
        self.left_click = False
        self.right_click = False
        self.middle_click = False
        self.pause = True
        self.show_vector = False
        self.show_boundary = False
        self.show_voxel = False
        self.record = False
        self.axis = True
        self.current_step = 0
        self.counter = 0
        self.camera = Camera()

        self.view = self.camera()
        self.view_changed = False

    def __call__(self, *args, **kwargs):
        glfw.make_context_current(self.window)

        self.coordinates = Coord()
        self.three_d_cursor = Coord(r"./3d_cursor.obj")

        self.demo = Demo2.Demo()
        glUseProgram(self.demo.render_shader_voxel)
        glUniformMatrix4fv(self.demo.voxel_projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, self.camera.view)
        glUseProgram(self.demo.render_shader)
        glUniformMatrix4fv(self.demo.projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.view_loc, 1, GL_FALSE, self.camera.view)
        glUseProgram(self.demo.render_shader_boundary)
        glUniformMatrix4fv(self.demo.boundary_projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.boundary_view_loc, 1, GL_FALSE, self.camera.view)
        glUseProgram(self.demo.render_shader_vector)
        glUniformMatrix4fv(self.demo.vector_projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.vector_view_loc, 1, GL_FALSE, self.camera.view)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_DEPTH_TEST)

        self.track_cursor()
        self.track_keyboard()

        glfw.show_window(self.window)
        i = 0
        while not glfw.window_should_close(self.window):
            glfw.poll_events()
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            # coord
            if self.axis:
                self.coordinates(projection_matrix=self.camera.projection, view_matrix=self.view,
                                 model_matrix=pyrr.matrix44.create_identity(np.float32))
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
                self.three_d_cursor(projection_matrix=self.camera.projection, view_matrix=self.view,
                                    model_matrix=pyrr.matrix44.create_from_translation(self.camera.lookat))
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

            # render codes
            self.demo(self.counter, pause=self.pause, show_vector=self.show_vector, show_boundary=self.show_boundary, show_voxel=self.show_voxel)
            if not self.pause:
                self.current_step += 1
            # glUseProgram(self.demo.compute_shader_a)
            # glUniform1i(self.demo.compute_shader_a_current_step_loc, self.current_step)
            """
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_particles)
            a0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.particles.nbytes),
                               dtype=np.float32)
            a = np.reshape(a0, (-1, 4))
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_boundary_particles)
            b0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.boundary_particles.nbytes),
                               dtype=np.float32)
            b = np.reshape(b0, (-1, 4))
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_voxels)
            c0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.voxels.nbytes),
                               dtype=np.int32)
            c = np.reshape(c0, (-1, 4))
            """
            """
            total = 0
            total_after = 0
            total_in = 0
            total_out = 0
            for i in range(self.demo.voxel_number):
                buf = c[i * 80: i * 80 + 80]
                _origin = buf[8:32].reshape((96,))
                _out = buf[32:56].reshape((96,))
                _in = buf[56:80].reshape((96,))
                n_ = 0
                n_o = 0
                n_i = 0
                for j in range(96):
                    if _origin[j] != 0:
                        n_ += 1
                    if _out[j] != 0:
                        n_o += 1
                    if _in[j] != 0:
                        n_i += 1
                print(n_, n_o, n_i)
                total += n_
                total_after += n_
                total_after += -n_o
                total_after += n_i
                total_in += n_i
                total_out += n_o
            print(total, total_after, total_in, total_out)
            """
            # count = 0
            # for item in a0:
            #     if item != 0:
            #         count +=1
            # print(count, self.demo.particle_number)
            # print(a)

            if self.view_changed:
                glUseProgram(self.demo.render_shader_voxel)
                glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, self.view)
                glUseProgram(self.demo.render_shader)
                glUniformMatrix4fv(self.demo.view_loc, 1, GL_FALSE, self.view)
                glUseProgram(self.demo.render_shader_boundary)
                glUniformMatrix4fv(self.demo.boundary_view_loc, 1, GL_FALSE, self.view)
                glUseProgram(self.demo.render_shader_vector)
                glUniformMatrix4fv(self.demo.vector_view_loc, 1, GL_FALSE, self.view)
                self.view_changed = False
            # time.sleep(0.02)
            if self.record:
                self.save_frames(f"tmp/{i}.jpg")
                # self.save_particle_data(i)
                i += 1

            glClearColor(0.0, 0.0, 0.0, 1.0)
            glfw.swap_buffers(self.window)
        glfw.terminate()

    def save_particle_data(self, i):
        import os
        os.makedirs("output", exist_ok=True)
        os.makedirs("output/group1", exist_ok=True)
        os.makedirs("output/group2", exist_ok=True)
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_particles)
        buffer = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.particles.nbytes),
                           dtype=np.float32).reshape((-1, 4))
        ptr = buffer.shape[0]//2
        self.save_as_ply(buffer[:ptr], f"output/group1/{i}.ply")
        self.save_as_ply(buffer[ptr:], f"output/group2/{i}.ply")

    @staticmethod
    def save_as_ply(particles, export_path):
        header = f"""ply\nformat ascii 1.0\nelement vertex {particles.shape[0]//4}\nproperty float x\nproperty float y\nproperty float z\nelement face 0\nproperty list uchar uint vertex_indices\nend_header\n"""
        with open(export_path, "w") as f:
            f.writelines(header)
            for i in range(particles.shape[0]//4):
                f.write("{} {} {}\n".format(*particles[i*4, :3]))
            f.close()

    @staticmethod
    def save_frames(filepath):
        x, y, width, height = glGetDoublev(GL_VIEWPORT)
        width, height = int(width), int(height)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        data = glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE)
        image = Image.frombytes("RGB", (width, height), data)
        image = image.transpose(Image.FLIP_TOP_BOTTOM)
        image.save(filepath, "JPEG")

    def track_cursor(self):
        def cursor_position_clb(*args):
            delta = np.array(args[1:], dtype=np.float32) - self.cursor_position[:]
            self.cursor_position = args[1:]
            if self.left_click:
                self.view = self.camera(pyrr.Vector3((*delta, 0.0)), "left")
                self.view_changed = True
                # glUseProgram(self.demo.render_shader_voxel)
                # glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, mat)
            elif self.middle_click:
                self.view = self.camera(pyrr.Vector3((-delta[0]*0.1, delta[1]*0.1, 0.0)), "middle")
                self.view_changed = True
                # glUseProgram(self.demo.render_shader_voxel)
                # glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, mat)

        def mouse_press_clb(window, button, action, mods):
            if button == glfw.MOUSE_BUTTON_LEFT and action == glfw.PRESS:
                self.left_click = True
                self.camera.mouse_left = True
            elif button == glfw.MOUSE_BUTTON_LEFT and action == glfw.RELEASE:
                self.left_click = False
                self.camera.mouse_left = False
            if button == glfw.MOUSE_BUTTON_RIGHT and action == glfw.PRESS:
                self.right_click = True
                self.camera.mouse_right = True
            elif button == glfw.MOUSE_BUTTON_RIGHT and action == glfw.RELEASE:
                self.right_click = False
                self.camera.mouse_right = False
            if button == glfw.MOUSE_BUTTON_MIDDLE and action == glfw.PRESS:
                self.middle_click = True
                self.camera.mouse_middle = True
            elif button == glfw.MOUSE_BUTTON_MIDDLE and action == glfw.RELEASE:
                self.middle_click = False
                self.camera.mouse_middle = False

        def scroll_clb(window, x_offset, y_offset):
            def f():
                # if sum([abs(item) for item in self.camera.position.xyz]) <= 1.01:
                #     if y_offset >= 0:
                #         return
                for i in range(20):
                    self.camera.position += self.camera.front * y_offset * 0.01
                    self.camera.position = pyrr.Vector4([*self.camera.position.xyz, 1.0])
                    self.view = self.camera(flag="wheel")
                    self.view_changed = True
                    time.sleep(0.005)
                    # if abs(sum([*self.camera.position.xyz])) <= 1.01:
                    #     if y_offset >= 0:
                    #         return

            t = threading.Thread(target=f)
            t.start()

        glfw.set_mouse_button_callback(self.window, mouse_press_clb)
        glfw.set_scroll_callback(self.window, scroll_clb)
        glfw.set_cursor_pos_callback(self.window, cursor_position_clb)

    def track_keyboard(self):
        def key_press_clb(window, key, scancode, action, mods):
            if key == glfw.KEY_SPACE and action == glfw.PRESS:
                self.pause = not self.pause
                if self.pause:
                    glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_particles)
                    a0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.particles.nbytes),
                                       dtype=np.float32)
                    a = np.reshape(a0, (-1, 4))
                    print(a[:16])
            if key == glfw.KEY_ENTER and action == glfw.PRESS:
                self.counter += 1
                self.counter %= self.demo.voxel_number
            if key == glfw.KEY_V and action == glfw.PRESS:
                self.show_vector = not self.show_vector
            if key == glfw.KEY_A and action == glfw.PRESS:
                if self.show_vector:
                    glProgramUniform1i(self.demo.render_shader_vector, self.demo.render_shader_vector_vector_type_loc, 1)
                glProgramUniform1i(self.demo.render_shader, self.demo.render_shader_color_type_loc, 1)
            if key == glfw.KEY_S and action == glfw.PRESS:
                if self.show_vector:
                    glProgramUniform1i(self.demo.render_shader_vector, self.demo.render_shader_vector_vector_type_loc, 0)
                glProgramUniform1i(self.demo.render_shader, self.demo.render_shader_color_type_loc, 0)
            if key == glfw.KEY_D and action == glfw.PRESS:
                glProgramUniform1i(self.demo.render_shader, self.demo.render_shader_color_type_loc, 3)
            if key == glfw.KEY_P and action == glfw.PRESS:
                glProgramUniform1i(self.demo.render_shader, self.demo.render_shader_color_type_loc, 2)
            if key == glfw.KEY_B and action == glfw.PRESS:
                self.show_boundary = not self.show_boundary
            if key == glfw.KEY_G and action == glfw.PRESS:
                self.show_voxel = not self.show_voxel
            if key == glfw.KEY_R and action == glfw.PRESS:
                self.record = not self.record
                if self.record:
                    print("Recording in progress.")
                else:
                    print("Recording stopped.")
            if key == glfw.KEY_C and action == glfw.PRESS:
                self.axis = not self.axis
        glfw.set_key_callback(self.window, key_press_clb)


if __name__ == "__main__":
    dp = DisplayPort()
    dp()

