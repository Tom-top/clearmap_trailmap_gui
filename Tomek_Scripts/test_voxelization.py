import numpy as np

diam = 10
xdiam = diam
ydiam = diam
zdiam = diam

points = np.array([(50, 50, 50), (50, 50, 50), (30, 30, 5),
                   (10, 10, 10), (20, 10, 50), (30, 30, 10), (20, 20, 5)])
xsize, ysize, zsize = 100, 100, 100
voximg = np.zeros([xsize, ysize, zsize], dtype=np.int32)

iCentroid = 0
nCentroid = points.shape[0]
nSphereIndices = int(xdiam * ydiam * zdiam)

# precompute indices centered at 0,0,0
xs = np.zeros([nSphereIndices], dtype=np.int)
ys = np.zeros([nSphereIndices], dtype=np.int)
zs = np.zeros([nSphereIndices], dtype=np.int)
buffer = np.zeros([xdiam, ydiam, zdiam], dtype=np.int)
ns = 0

xdiam2 = (xdiam - 1) * (xdiam - 1) / 4
ydiam2 = (ydiam - 1) * (ydiam - 1) / 4
zdiam2 = (zdiam - 1) * (zdiam - 1) / 4

plt.figure()
x_plot_data = np.arange(int(-(xdiam / 2 + 1)), int(xdiam / 2 + 1) + 1, 1)
y_plot_data = (x_plot_data - 1) * (x_plot_data - 1) / 4
plt.plot(x_plot_data, y_plot_data)
plt.plot((-(xdiam / 2), (xdiam / 2)), (1,1), color='red' )

for x in range(int(-(xdiam / 2 + 1)), int(xdiam / 2 + 1)):
    for y in range(int(-(ydiam / 2 + 1)), int(ydiam / 2 + 1)):
        for z in range(int(-(zdiam / 2 + 1)), int(zdiam / 2 + 1)):
            if x * x / xdiam2 + y * y / ydiam2 + z * z / zdiam2 < 1:
                # buffer[ns, ns, ns] = x
                xs[ns] = x
                ys[ns] = y
                zs[ns] = z
                ns += 1

# for x in range(int(-xdiam / 2 + 1), int(xdiam / 2 + 1)):
#     for y in range(int(-ydiam / 2 + 1), int(ydiam / 2 + 1)):
#         for z in range(int(-zdiam / 2 + 1), int(zdiam / 2 + 1)):
#             if x * x / xdiam2 + y * y / ydiam2 + z * z / zdiam2 < 1:
#                 # buffer[ns, ns, ns] = x
#                 xs[ns] = x
#                 ys[ns] = y
#                 zs[ns] = z
#                 ns += 1

import matplotlib.pyplot as plt
plt.figure()
ax0 = plt.subplot(311)
ax0.plot(xs)
ax1 = plt.subplot(312)
ax1.plot(ys)
ax2 = plt.subplot(313)
ax2.plot(zs)
plt.show()

iss = 0
cx0 = 0.0
cy0 = 0.0
cz0 = 0.0

cxf = 0.0
cyf = 0.0
czf = 0.0

cx = 0
cy = 0
cz = 0

for iCentroid in range(nCentroid):
    if ((iCentroid % 25000) == 0):
        print
        "\nProcessed %d/%d\n" % (iCentroid, nCentroid);

    cx0 = points[iCentroid, 0]
    cy0 = points[iCentroid, 1]
    cz0 = points[iCentroid, 2]

    for iss in range(ns):
        cxf = cx0 + xs[iss]
        cyf = cy0 + ys[iss]
        czf = cz0 + zs[iss]

        if cxf >= 0 and cxf < xsize:
            if cyf >= 0 and cyf < ysize:
                if czf >= 0 and czf < zsize:
                    cx = int(cxf)
                    cy = int(cyf)
                    cz = int(czf)

                    voximg[cx, cy, cz] = voximg[cx, cy, cz] + 1

pos_data = np.argwhere(voximg > 0)
# values = [voximg[pos_data[n]] for n in pos_data]
plt.figure()
ax0 = plt.axes(projection='3d')
ax0.set_box_aspect([1,1,1])
ax0.scatter(pos_data[:,0],
            pos_data[:,1],
            pos_data[:,2],
            c=voximg[pos_data[:,0],
                     pos_data[:,1],
                     pos_data[:,2]],
            cmap='viridis',
            s=1,
            linewidth=0.5,
            alpha=0.5)
ax0.set_xlim3d(0, xsize)
ax0.set_ylim3d(0, ysize)
ax0.set_zlim3d(0, ysize)
plt.show()