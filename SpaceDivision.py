import math
import numpy as np

H = 1.0
Domain = [[0, 0, 0], [0, 3, 0], [0, 0, 3], [0, 3, 3], [4, 0, 0], [4, 0, 3], [4, 3, 0], [4, 3, 3]]

Domain = np.array(Domain)
x = math.ceil((max(Domain[:, 0]) - min(Domain[:, 0])) / H) + 1
y = math.ceil((max(Domain[:, 1]) - min(Domain[:, 1])) / H) + 1
z = math.ceil((max(Domain[:, 2]) - min(Domain[:, 2])) / H) + 1
n = x * y * z

DomainMat = np.zeros((n, 80, 4))
for i in range(n):
    index = i + 1
    DomainMat[i, 0, :] = [index, (i % (x * y) // y) * H, (i % (x * y) % y) * H, i // (x * y) * H]
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
            buf = np.array([index, index - 1, index + 1, index + y, index + y - 1, index + y + 1], dtype=np.float32)
        elif pt // y == x - 1:
            buf = np.array([index, index - 1, index + 1, index - y, index - y - 1, index - y + 1], dtype=np.float32)
        elif pt % y == 0:
            buf = np.array([index, index + 1, index - y, index - y + 1, index + y, index + y + 1], dtype=np.float32)
        elif pt % y == x - 1:
            buf = np.array([index, index - 1, index - y, index - y - 1, index + y, index + y - 1], dtype=np.float32)
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

    DomainMat[i, 1:8, :] = contents
    # print(DomainMat[i, 0:8, :])
