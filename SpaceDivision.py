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

        domain_mat = np.zeros((n, 80, 4), dtype=np.float32)
        for i in range(n):
            index = i + 1
            domain_mat[i, 0, :] = [index, (i % (x * y) // y) * self.h, (i % (x * y) % y) * self.h,
                                   i // (x * y) * self.h]
            front = max(0, index - (x * y))
            back = 0 if index + (x * y) > n else index + (x * y)
            pt = i % (x * y)

            if pt == 0:
                buf = np.array([index, index + 1, index + y, index + y + 1], dtype=np.float32)
            elif pt == y - 1:
                buf = np.array([index, index - 1, index + y, index + y - 1], dtype=np.float32)
            elif pt == x * y - y:
                buf = np.array([index, index + 1, index - y, index - y + 1], dtype=np.float32)
            elif pt == x * y - 1:
                buf = np.array([index, index - 1, index - y, index - y - 1], dtype=np.float32)
            else:
                if pt // y == 0:
                    buf = np.array([index, index - 1, index + 1, index + y, index + y - 1, index + y + 1],
                                   dtype=np.float32)
                elif pt // y == x - 1:
                    buf = np.array([index, index - 1, index + 1, index - y, index - y - 1, index - y + 1],
                                   dtype=np.float32)
                elif pt % y == 0:
                    buf = np.array([index, index + 1, index - y, index - y + 1, index + y, index + y + 1],
                                   dtype=np.float32)
                elif pt % y == x - 1:
                    buf = np.array([index, index - 1, index - y, index - y - 1, index + y, index + y - 1],
                                   dtype=np.float32)
                else:
                    buf = np.array([index, index - 1, index + 1,
                                    index - y, index - y - 1, index - y + 1,
                                    index + y, index + y - 1, index + y + 1], dtype=np.float32)

            contents = buf
            if front != 0:
                contents = np.hstack((contents, buf - x * y))
            if back != 0:
                contents = np.hstack((contents, buf + x * y))
            contents = np.pad(contents[1:], [(0, 29 - len(contents))])
            contents = contents.reshape((7, 4))

            domain_mat[i, 1:8, :] = contents
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
                    domain_particle_mat[particle_id][0, :3] = start_point + np.array((i*self.r*2+self.r, j*self.r*2+self.r, k*self.r*2+self.r), dtype=np.float32) + ((np.random.ranf(3)-0.50)*2)*(0.1*self.r)
                    # assign initial mass
                    domain_particle_mat[particle_id][1, 3] = 1.0
                    # assign initial velocity

                    # assign initial density
                    # assign initial pressure
                    # assign initial color
        output_domain_particle_mat = np.vstack((particle_mat for particle_mat in domain_particle_mat))
        return output_domain_particle_mat

    def generate_boundary_particle(self):
        ...


if __name__ == "__main__":
    H = 1.0
    R = 0.2
    Domain = [[-2, -1, -3], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [3, 1, 3]]  # 8x3
    voxels = CreateVoxels(domain=Domain, h=H)()
    print(voxels)
    print(voxels.shape)
    particles = CreateParticles(domain=Domain, h=H, r=R)()
    print(particles)
    print(particles.shape)
