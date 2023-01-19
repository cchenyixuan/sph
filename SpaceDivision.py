import math
import numpy as np


class CreateVoxels:
    def __init__(self, domain, h):
        self.domain = np.array(domain, dtype=np.float32)
        self.h = h

    def __call__(self, *args, **kwargs):
        return self.space_division()

    def space_division(self):
        x = math.ceil((max(self.domain[:, 0]) - min(self.domain[:, 0])) / self.h) + 1
        y = math.ceil((max(self.domain[:, 1]) - min(self.domain[:, 1])) / self.h) + 1
        z = math.ceil((max(self.domain[:, 2]) - min(self.domain[:, 2])) / self.h) + 1
        n = x * y * z

        domain_mat = np.zeros((n, 182*4, 4), dtype=np.int32)
        for i in range(n):
            index = i + 1
            # float version
            # domain_mat[i, 0, :] = [index, (i % (x * y) // y) * self.h, (i % (x * y) % y) * self.h,
            #                        i // (x * y) * self.h]
            # int version
            domain_mat[i, 0, :] = [index, (i % (x * y) % y), (i % (x * y) // y), i // (x * y)]
            back = max(0, index - (x * y))
            front = 0 if index + (x * y) > n else index + (x * y)
            pt = i % (x * y)

            if pt == 0:
                buf = np.array([0, index + 1, 0, index + y, 0, 0, 0, index + y + 1], dtype=np.int32)
            elif pt == y - 1:
                buf = np.array([index - 1, 0, 0, index + y, 0, index + y - 1, 0, 0], dtype=np.int32)
            elif pt == x * y - y:
                buf = np.array([0, index + 1, index - y, 0, 0, 0, index - y + 1, 0], dtype=np.int32)
            elif pt == x * y - 1:
                buf = np.array([index - 1, 0,index - y, 0, index - y - 1, 0,  0, 0], dtype=np.int32)
            else:
                if pt // y == 0:
                    buf = np.array([index - 1, index + 1, 0, index + y, 0, index + y - 1, 0, index + y + 1],
                                   dtype=np.int32)
                elif pt // y == x - 1:
                    buf = np.array([index - 1, index + 1, index - y, 0, index - y - 1, 0, index - y + 1, 0],
                                   dtype=np.int32)
                elif pt % y == 0:
                    buf = np.array([0, index + 1, index - y, index + y, 0, 0, index - y + 1, index + y + 1],
                                   dtype=np.int32)
                elif pt % y == x - 1:
                    buf = np.array([index - 1, 0, index - y, index + y, index - y - 1, index + y - 1, 0, 0],
                                   dtype=np.int32)
                else:
                    buf = np.array([index - 1, index + 1, index - y, index + y,
                                    index - y - 1, index + y - 1, index - y + 1, index + y + 1], dtype=np.int32)

            contents = np.zeros((2, 8), dtype=np.float32)
            if back != 0:
                contents[0, :] = [pos - x * y if pos != 0 else 0 for pos in buf]
            if front != 0:
                contents[1, :] = [pos + x * y if pos != 0 else 0 for pos in buf]
            contents = contents.T
            contents = contents.reshape(16)
            contents = np.hstack((buf[:4], back, front, buf[4:], contents, 0, 0))
            contents = contents.reshape((7, 4))
            domain_mat[i, 1:8, :] = contents

            # re-arrange
            """
            we set x, y, z have 3 status: -1, 0, 1
            Left = (-1, 0, 0)
            Right = (1, 0, 0)
            Down = (0, -1, 0)
            Up = (0, 1, 0)
            Back = (0, 0, -1)
            Front = (0, 0, 1)
            and other combinations of above, i.e.
            LeftUpBack = (-1, 1, -1) = Left + Up + Back
            
            re_arrange = [
                ["i", "x", "y", "z"],0-3
                [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0)],4-7
                [(0, 0, -1), (0, 0, 1), (-1, -1, 0), (-1, 1, 0)],8-11
                [(1, -1, 0), (1, 1, 0), (-1, 0, -1), (-1, 0, 1)],12-15
                [(1, 0, -1), (1, 0, 1), (0, -1, -1), (0, -1, 1)],16-19
                [(0, 1, -1), (0, 1, 1), (-1, -1, -1), (-1, -1, 1)],20-23
                [(-1, 1, -1), (-1, 1, 1), (1, -1, -1)", (1, -1, 1)],24-27
                [(1, 1, -1), (1, 1, 1), "0", "0"],28-31
            ]
            
            """
            re_arrange = [
                ["i", "x", "y", "z"],
                ["Left", "Right", "Down", "Up"],
                ["Back", "Front", "LeftDown", "LeftUp"],
                ["RightDown", "RightUp", "LeftBack", "LeftFront"],
                ["RightBack", "RightFront", "DownBack", "DownFront"],
                ["UpBack", "UpFront", "LeftDownBack", "LeftDownFront"],
                ["LeftUpBack", "LeftUpFront", "RightDownBack", "RightDownFront"],
                ["RightUpBack", "RightUpFront", "0", "0"],
            ]
        output_buffer = np.vstack((voxel_matrices for voxel_matrices in domain_mat))
        return output_buffer


class CreateParticles:
    def __init__(self, domain, h, r):
        self.domain = np.array(domain, dtype=np.float32)
        self.h = h
        self.r = r

    def __call__(self, *args, **kwargs):
        return self.generate_domain_particle()

    def generate_domain_particle(self):
        x_start = min(self.domain[:, 0]) + self.r
        x_end = max(self.domain[:, 0]) - self.r

        y_start = min(self.domain[:, 1]) + self.r
        y_end = max(self.domain[:, 1]) - self.r

        z_start = min(self.domain[:, 2]) + self.r
        z_end = max(self.domain[:, 2]) - self.r

        n_x = int((x_end - x_start) // (2 * self.r))
        n_y = int((y_end - y_start) // (2 * self.r))
        n_z = int((z_end - z_start) // (2 * self.r))

        start_point = np.array((x_start, y_start, z_start), dtype=np.float32)

        domain_particle_mat = np.zeros((n_x*n_y*n_z, 4, 4), dtype=np.float32)

        for i in range(n_x):
            for j in range(n_y):
                for k in range(n_z):
                    # particle id
                    particle_id = i*n_y*n_z + j*n_z + k
                    # assign a particle to its location with a random offset in scale [-r/10, r/10]
                    domain_particle_mat[particle_id][0, :3] = start_point + np.array((i*self.r*2+self.r, j*self.r*2+self.r, k*self.r*2+self.r), dtype=np.float32)# + ((np.random.ranf(3)-0.50)*2)*(0.1*self.r)
                    # assign initial mass
                    domain_particle_mat[particle_id][1, 3] = 0.27
                    # assign initial velocity

                    # assign initial density
                    # assign initial pressure
                    # assign initial color
        output_domain_particle_mat = np.vstack((particle_mat for particle_mat in domain_particle_mat))
        return output_domain_particle_mat


class CreateBoundaryParticles:
    def __init__(self, domain, h, r):
        self.domain = np.array(domain, dtype=np.float32)
        self.h = h
        self.r = r
        self.d = 2*self.r

    def __call__(self, *args, **kwargs):
        return self.generate_boundary_particle_iterative(4)

    def generate_boundary_particle_iterative(self, layers=1):
        boundary_particles = []
        min_x = np.min(self.domain[:, 0])
        min_y = np.min(self.domain[:, 1])
        min_z = np.min(self.domain[:, 2])
        max_x = np.max(self.domain[:, 0])
        max_y = np.max(self.domain[:, 1])
        max_z = np.max(self.domain[:, 2])
        for layer in range(layers):
            min_x -= self.d
            min_y -= self.d
            min_z -= self.d
            max_x += self.d
            max_y += self.d
            max_z += self.d
            quantity_x = round((max_x - min_x) / self.d + 1)
            quantity_y = round((max_y - min_y) / self.d + 1)
            quantity_z = round((max_z - min_z) / self.d + 1)
            # top bottom
            for i in [0, quantity_y - 1]:
                for j in range(quantity_z):
                    boundary_particles.append(
                        np.array([(min_x + _ * self.d, min_y + i * self.d, min_z + j * self.d) for _ in range(quantity_x)],
                                 dtype=np.float32))
            # front back
            for j in [0, quantity_z - 1]:
                for i in range(1, quantity_y - 1):
                    boundary_particles.append(
                        np.array([(min_x + _ * self.d, min_y + i * self.d, min_z + j * self.d) for _ in range(quantity_x)],
                                 dtype=np.float32))
            # left right
            for i in [0, quantity_x - 1]:
                for j in range(1, quantity_y - 1):
                    boundary_particles.append(np.array(
                        [(min_x + i * self.d, min_y + j * self.d, min_z + _ * self.d) for _ in range(1, quantity_z - 1)],
                        dtype=np.float32))

        boundary_particles = np.vstack(boundary_particles)
        buffer = np.zeros((boundary_particles.shape[0] * 4, 4), dtype=np.float32)
        for step, item in enumerate(boundary_particles):
            buffer[step * 4][:3] = item
            buffer[step * 4 + 1][3] = 10.0
            buffer[step * 4 + 2][2:] = np.array((1500, 1000), dtype=np.float32)
            buffer[step * 4 + 3][:] = np.array((1.0, 1.0, 1.0, 1.0), dtype=np.float32)

        return buffer

    def generate_boundary_particle(self):
        boundary_particles = []
        min_x = np.min(self.domain[:, 0])
        min_y = np.min(self.domain[:, 1])
        min_z = np.min(self.domain[:, 2])
        max_x = np.max(self.domain[:, 0])
        max_y = np.max(self.domain[:, 1])
        max_z = np.max(self.domain[:, 2])
        quantity_x = round((max_x - min_x) / self.d + 1)
        quantity_y = round((max_y - min_y) / self.d + 1)
        quantity_z = round((max_z - min_z) / self.d + 1)
        # top bottom
        for i in [0, quantity_y-1]:
            for j in range(quantity_z):
                boundary_particles.append(np.array([(min_x+_*self.d, min_y+i*self.d, min_z+j*self.d) for _ in range(quantity_x)], dtype=np.float32))
        # front back
        for j in [0, quantity_z-1]:
            for i in range(1, quantity_y-1):
                boundary_particles.append(np.array([(min_x+_*self.d, min_y+i*self.d, min_z+j*self.d) for _ in range(quantity_x)], dtype=np.float32))
        # left right
        for i in [0, quantity_x-1]:
            for j in range(1, quantity_y-1):
                boundary_particles.append(np.array([(min_x+i*self.d, min_y+j*self.d, min_z+_*self.d) for _ in range(1, quantity_z-1)], dtype=np.float32))

        boundary_particles = np.vstack(boundary_particles)
        buffer = np.zeros((boundary_particles.shape[0]*4, 4), dtype=np.float32)
        for step, item in enumerate(boundary_particles):
            buffer[step*4][:3] = item
            buffer[step*4+1][3] = 1000.0
            buffer[step*4+2][2:] = np.array((500, 1000), dtype=np.float32)
            buffer[step * 4 + 3][:] = np.array((1.0, 1.0, 1.0, 1.0), dtype=np.float32)

        return buffer


if __name__ == "__main__":
    H = 0.25
    R = 0.2
    Domain = [[0,0,0], [1,1,1]]  # 8x3
    voxels = CreateVoxels(domain=Domain, h=H)()
    print(voxels)
    print(voxels.shape)
    particles = CreateParticles(domain=Domain, h=H, r=R)()
    print(particles)
    print(particles.shape)
