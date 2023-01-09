import time
import threading

import pyrr
import numpy as np
from OpenGL.GL import *
import glfw
import Demo
from camera import Camera



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
        self.camera = Camera()

        self.view = self.camera()
        self.view_changed = False

    def __call__(self, *args, **kwargs):
        glfw.make_context_current(self.window)

        self.demo = Demo.Demo()
        glUseProgram(self.demo.render_shader_voxel)
        glUniformMatrix4fv(self.demo.voxel_projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, self.camera.view)
        glUseProgram(self.demo.render_shader)
        glUniformMatrix4fv(self.demo.projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.view_loc, 1, GL_FALSE, self.camera.view)
        glUseProgram(self.demo.render_shader_boundary)
        glUniformMatrix4fv(self.demo.boundary_projection_loc, 1, GL_FALSE, self.camera.projection)
        glUniformMatrix4fv(self.demo.boundary_view_loc, 1, GL_FALSE, self.camera.view)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_DEPTH_TEST)

        self.track_cursor()

        glfw.show_window(self.window)

        i = 0
        while not glfw.window_should_close(self.window):
            glfw.poll_events()
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

            # render codes
            self.demo(i//100)
            i += 1
            if i > self.demo.voxel_number*100:
                i = 0
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_particles)
            a0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.particles.nbytes), dtype=np.float32)
            a = np.reshape(a0, (-1, 4))
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_boundary_particles)
            b0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.boundary_particles.nbytes),
                               dtype=np.float32)
            b = np.reshape(b0, (-1, 4))
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_voxels)
            c0 = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.voxels.nbytes),
                               dtype=np.int32)
            c = np.reshape(c0, (-1, 4))
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
                self.view_changed = False
            # time.sleep(0.02)

            glClearColor(0.0, 0.0, 0.0, 1.0)
            glfw.swap_buffers(self.window)
        glfw.terminate()

    def track_cursor(self):
        def cursor_position_clb(*args):
            print(args[1:])
            delta = np.array(args[1:], dtype=np.float32) - self.cursor_position[:]
            self.cursor_position = args[1:]
            if self.left_click:
                self.view = self.camera(pyrr.Vector3((*delta, 0.0)), "left")
                self.view_changed = True
                # glUseProgram(self.demo.render_shader_voxel)
                # glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, mat)
            elif self.middle_click:
                self.view = self.camera(pyrr.Vector3((-delta[0], delta[1], 0.0)), "middle")
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
                if sum([abs(item) for item in self.camera.position.xyz]) <= 1.01:
                    if y_offset >= 0:
                        return
                for i in range(50):
                    self.camera.position += self.camera.front * y_offset * 0.02
                    self.camera.position = pyrr.Vector4([*self.camera.position.xyz, 1.0])
                    self.view = self.camera(flag="wheel")
                    self.view_changed = True
                    time.sleep(0.005)
                    if abs(sum([*self.camera.position.xyz])) <= 1.01:
                        if y_offset >= 0:
                            return

            t = threading.Thread(target=f)
            t.start()

        glfw.set_mouse_button_callback(self.window, mouse_press_clb)
        glfw.set_scroll_callback(self.window, scroll_clb)
        glfw.set_cursor_pos_callback(self.window, cursor_position_clb)


if __name__ == "__main__":
    dp = DisplayPort()
    dp()

