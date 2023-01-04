import time
import threading

import pyrr
import numpy as np
from OpenGL.GL import *
import glfw
import Demo



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
        self.click = False



    def __call__(self, *args, **kwargs):
        glfw.make_context_current(self.window)

        self.demo = Demo.Demo()
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
            # glBindBuffer(GL_SHADER_STORAGE_BUFFER, self.demo.sbo_particles)
            # a = np.frombuffer(glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, self.demo.particles.nbytes), dtype=np.float32)
            # a = np.reshape(a, (-1, 4))
            # # b = self.demo.voxels
            # print(a)

            glClearColor(0.0, 0.0, 0.0, 1.0)
            glfw.swap_buffers(self.window)
        glfw.terminate()

    def track_cursor(self):
        def cursor_position_clb(*args):
            print(args[1:])
            delta = np.array(args[1:], dtype=np.float32) - self.cursor_position[:]
            self.cursor_position = args[1:]
            if self.click:
                mat = self.demo.camera(pyrr.Vector3((*delta, 0.0)), "left")
                glUseProgram(self.demo.render_shader_voxel)
                glUniformMatrix4fv(self.demo.voxel_view_loc, 1, GL_FALSE, mat)
                glUseProgram(self.demo.render_shader)
                glUniformMatrix4fv(self.demo.view_loc, 1, GL_FALSE, mat)

        def mouse_button_clb(*args):
            if args[1] == 0 and args[2] == 1:
                self.click = True
                self.demo.camera.mouse_left = True
            else:
                self.click = False
                self.demo.camera.mouse_left = False

        glfw.set_mouse_button_callback(self.window, mouse_button_clb)

        glfw.set_cursor_pos_callback(self.window, cursor_position_clb)


if __name__ == "__main__":
    dp = DisplayPort()
    dp()

