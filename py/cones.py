#------------------------------------------------------------------------------------------
# Gaussian-Product Subdivision Surfaces - Python Demo
# Copyright (c) 2020 Reinhold Preiner
#
# This code is licensed under the MIT license. 
# https://opensource.org/licenses/mit-license.html
#------------------------------------------------------------------------------------------


import numpy as np
import igl
import os
root_folder = ".."


def loop_gps(v, f, c, n):
    # transform to 9D dual-space vertices
    qq1 = np.zeros((len(v), 3))
    qq2 = np.zeros((len(v), 3))
    qlin = np.zeros((len(v), 3))
    for i, cov in enumerate(c):
        ic = np.linalg.inv(cov)
        icf = ic.flatten()
        qq1[i] = [icf[0],icf[1],icf[2]]
        qq2[i] = [icf[4],icf[5],icf[8]]
        qlin[i] = ic @ v[i]

    # perform Gaussian-product subdivision
    # note: igl.loop only handles 3D subdivs, so we split the 9D meshes into three 3D ones
    for _ in range(n):
        qq1, dmy = igl.loop(qq1, f)
        qq2, dmy = igl.loop(qq2, f)
        qlin, f  = igl.loop(qlin, f)
    
    # transform back to 3D
    v = np.zeros((len(qlin), 3))
    for i, ql in enumerate(qlin):
        icov = [qq1[i],
                [qq1[i][1], qq2[i][0], qq2[i][1]],
                [qq1[i][2], qq2[i][1], qq2[i][2]]]
        v[i] = np.linalg.inv(icov) @ ql
        
    return v, f


#------------------------------------------------------------------------------------------

# create a cone
v = np.array([
    [0, 1.5, 0],
    [1.41421, -1.5, -1.41421],[0, -1.5, -2], [-1.41421, -1.5, -1.41421],[-2, -1.5, 0],
    [-1.41421, -1.5, 1.41421], [0, -1.5, 2], [1.41421, -1.5, 1.41421],[2, -1.5, 0]
])
f = np.array([[0, 1, 2],[0, 2, 3], [0, 3, 4], [0, 4, 5],[0, 5, 6],[0, 6, 7],[0, 7, 8],[0, 8, 1]])



# initialize with isotropic covariances
c = [np.identity(3)*0.01 for i in range(len(v))]
vs, fs = loop_gps(v, f, c, 5)

# make apex covariance sharper
c[0] *= 0.001
vs2, fs2 = loop_gps(v, f, c, 5)

# make apex covariance vertically anisotropic (reduce variance in x direction)
c[0] = [[0.001,0,0],
        [0,0.01,0],
        [0,0,0.01]]
vs3, fs3 = loop_gps(v, f, c, 5)

# make apex covariance horizontally anisotropic (reduce variance in y direction)
c[0] = [[0.01,0,0],
        [0,0.001,0],
        [0,0,0.01]]
vs4, fs4 = loop_gps(v, f, c, 5)

# make apex covariance horizontally anisotropic and super flat
c[0] = [[1,0,0],
        [0,0.00001,0],
        [0,0,1]]
vs5, fs5 = loop_gps(v, f, c, 5)



# save meshes to .off files
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone.off"), v, f)
print("Wrote control mesh to", os.path.join(root_folder, "data", "cone.off"))

igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-pointy-gps.off"), vs2, fs2)
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-concave-gps.off"), vs3, fs3)
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-plateau-gps.off"), vs4, fs4)
igl.write_triangle_mesh(os.path.join(root_folder, "data", "cone-cylinder-gps.off"), vs5, fs5)
print("Wrote GPS meshes to", os.path.join(root_folder, "data", "cone-<variant>-gps.off"))
