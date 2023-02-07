import numpy as np
from SpaceDivision import CreateVoxels


class Project:
    def __init__(self, h, frame, domain, boundary, moving_parts):
        self.h = h
        self.r = h/4
        self.rho = 1000.0
        self.particle_mass = (2*self.r)**3*self.rho * 1.1  # scaled by 1.5
        self.frame = self.get_geometry(self.load_file(frame))
        self.voxels = CreateVoxels(self.frame, h)()
        self.particles = self.load_domain(self.load_file(domain))
        self.boundary_particles = self.load_boundary(self.load_file(boundary))
        self.pump_particles = self.load_pump_boundary(self.load_file(moving_parts))

    @staticmethod
    def load_file(file):
        import re
        find_vertex = re.compile(r"v (\+?-?[\d.]+) (\+?-?[\d.]+) (\+?-?[\d.]+)\n", re.S)
        data = []
        with open(file, "r") as f:
            for row in f:
                ans = re.findall(find_vertex, row)
                if ans:
                    ans = [float(ans[0][i]) for i in range(3)]
                    data.append([ans[0], ans[1], -ans[2]])
            f.close()
        return np.array(data, dtype=np.float32)

    @staticmethod
    def get_geometry(particles):
        return np.array([
            [np.min(particles[:, 0]), np.min(particles[:, 1]), np.min(particles[:, 2])],
            [np.max(particles[:, 0]), np.max(particles[:, 1]), np.max(particles[:, 2])]
        ], dtype=np.float32)

    def load_domain(self, particles):
        output = np.zeros((particles.shape[0]*4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step*4][:3] = vertex
            output[step*4+1][3] = self.particle_mass
            output[step*4+2][3] = self.rho
            output[step*4+3][3] = 0.0  # initial pressure
        return output

    def load_boundary(self, particles):
        output = np.zeros((particles.shape[0] * 4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step * 4][:3] = vertex
            output[step * 4 + 1][3] = self.particle_mass*4  # boundary has 4 times mass as usual particles
            output[step * 4 + 2][3] = self.rho
            output[step * 4 + 3][3] = 1000.0  # initial pressure
        return output

    def load_moving_boundary(self, particles):
        output = np.zeros((particles.shape[0] * 4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step * 4][:3] = vertex
            output[step * 4 + 1][3] = self.particle_mass * 4  # boundary has 4 times mass as usual particles
            output[step * 4 + 2][3] = self.rho
            output[step * 4 + 3][3] = 100.0  # initial pressure

            output[step * 4 + 2][:3] = vertex  # original position
        return output

    def load_pump_boundary(self, particles):
        output = np.zeros((particles.shape[0] * 4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step * 4][:3] = vertex
            output[step * 4 + 1][3] = 0.0  # boundary has 4 times mass as usual particles
            output[step * 4 + 2][3] = self.rho
            output[step * 4 + 3][3] = 100.0  # initial pressure

            output[step * 4 + 2][0] = 1.0  # type
            output[step * 4 + 2][1] = round(abs(vertex[0]-1.06)*100)  # group_id

        return output


class ProjectTwoPhase:
    def __init__(self, h, frame, domain1, domain2, boundary):
        self.h = h
        self.r = h/4
        self.rho = 1000.0
        self.particle_mass = (2*self.r)**3*self.rho * 1.1  # scaled by 1.5
        self.frame = self.get_geometry(self.load_file(frame))
        self.voxels = CreateVoxels(self.frame, h)()
        self.particles_1 = self.load_domain(self.load_file(domain1))
        self.particles_2 = self.load_domain(self.load_file(domain2))
        self.particles_1_buffer = self.create_particle_sub_buffer(self.particles_1, 1)
        self.particles_2_buffer = self.create_particle_sub_buffer(self.particles_2, 2)
        self.boundary_particles = self.load_boundary(self.load_file(boundary))

    @staticmethod
    def create_particle_sub_buffer(particles, group_id):
        buffer = np.zeros_like(particles)
        for i in range(buffer.shape[0]//4):
            buffer[i*4+3][-1] = group_id
        return buffer

    @staticmethod
    def load_file(file):
        import re
        find_vertex = re.compile(r"v (\+?-?[\d.]+) (\+?-?[\d.]+) (\+?-?[\d.]+)\n", re.S)
        data = []
        with open(file, "r") as f:
            for row in f:
                ans = re.findall(find_vertex, row)
                if ans:
                    ans = [float(ans[0][i]) for i in range(3)]
                    data.append([ans[0], ans[1], -ans[2]])
            f.close()
        return np.array(data, dtype=np.float32)

    @staticmethod
    def get_geometry(particles):
        return np.array([
            [np.min(particles[:, 0]), np.min(particles[:, 1]), np.min(particles[:, 2])],
            [np.max(particles[:, 0]), np.max(particles[:, 1]), np.max(particles[:, 2])]
        ], dtype=np.float32)

    def load_domain(self, particles):
        output = np.zeros((particles.shape[0]*4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step*4][:3] = vertex
            output[step*4+1][3] = self.particle_mass
            output[step*4+2][3] = self.rho
            output[step*4+3][3] = 0.0  # initial pressure
        return output

    def load_boundary(self, particles):
        output = np.zeros((particles.shape[0] * 4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step * 4][:3] = vertex
            output[step * 4 + 1][3] = self.particle_mass*4  # boundary has 4 times mass as usual particles
            output[step * 4 + 2][3] = self.rho
            output[step * 4 + 3][3] = 1000.0  # initial pressure
        return output

    def load_moving_boundary(self, particles):
        output = np.zeros((particles.shape[0] * 4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step * 4][:3] = vertex
            output[step * 4 + 1][3] = self.particle_mass * 4  # boundary has 4 times mass as usual particles
            output[step * 4 + 2][3] = self.rho
            output[step * 4 + 3][3] = 100.0  # initial pressure

            output[step * 4 + 2][:3] = vertex  # original position
        return output

    def load_pump_boundary(self, particles):
        output = np.zeros((particles.shape[0] * 4, 4), dtype=np.float32)
        for step, vertex in enumerate(particles):
            output[step * 4][:3] = vertex
            output[step * 4 + 1][3] = 0.0  # boundary has 4 times mass as usual particles
            output[step * 4 + 2][3] = self.rho
            output[step * 4 + 3][3] = 100.0  # initial pressure

            output[step * 4 + 2][0] = 1.0  # type
            output[step * 4 + 2][1] = round(abs(vertex[0]-1.06)*100)  # group_id

        return output


if __name__ == "__main__":
    # project = Project(0.02, r"./models/frame_recurrent.obj", r"./models/domain_recurrent.obj", r"./models/pump_cover_recurrent.obj", r"./models/pump_slice_recurrent.obj")
    project = ProjectTwoPhase(0.02, r"./models/CUBE_FRAME.obj", r"./models/phase1.obj",
                      r"./models/phase2.obj", r"./models/CUBE_COVER.obj")