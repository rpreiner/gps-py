#------------------------------------------------------------------------------------------
# Gaussian-Product Subdivision Surfaces - Python Demo
# Copyright (c) 2020 Reinhold Preiner
#
# This code is licensed under the MIT license. 
# https://opensource.org/licenses/mit-license.html
#------------------------------------------------------------------------------------------

import igl
import scipy as sp
import numpy as np
import os
import copy
root_folder = ".."


# define mesh struct
class Mesh: 
    v = []
    f = []   
    
#------------------------------------------------------------------------------------------


# load input triangle mesh
mesh = Mesh()
mesh.v, mesh.f = igl.read_triangle_mesh(os.path.join(root_folder, "data", "tweety.off"))


# perform ordinary Loop subdivison
lmesh = copy.deepcopy(mesh);
for _ in range(4):
    lmesh.v, lmesh.f = igl.loop(lmesh.v, lmesh.f)


#------------------------------------------------------------------------------------------

# Gaussian Inference via product of face Gaussians (Eq. 16 and 17)

cm = Mesh()
cm.v = mesh.v.copy()
cm.f = mesh.f.copy()

# create list of empty (inverse) covariances
icov = np.zeros((len(mesh.v),3,3))


# 1. infer face covs and inverse vertex covs
for f in cm.f:
    # compute covariance of face vertices
    cov = np.zeros((3,3))
    sum = np.zeros(3)
    for j in [1,2]:
        v = cm.v[f[j]] - cm.v[f[0]]
        sum += v
        cov += np.outer(v,v)
    cov = cov/3 - np.outer(sum,sum)/9
    # bias covariance by some fraction of its dominant eigenvalue
    bias = np.linalg.eigvalsh(cov)[2] * 0.05
    cov += np.identity(3) * bias
    # inverse cov at vertices is given by the sum of inverse of surrounding face covs
    for fv in f:
        icov[fv] += np.linalg.inv(cov)
    

# 2. transform to 9D dual-space vertices
qq1 = np.zeros((len(cm.v), 3))
qq2 = np.zeros((len(cm.v), 3))
qlin = np.zeros((len(cm.v), 3))
for i, ic in enumerate(icov):
    icf = ic.flatten()
    qq1[i] = [icf[0],icf[1],icf[2]]
    qq2[i] = [icf[4],icf[5],icf[8]]
    qlin[i] = ic @ cm.v[i]


# 3. perform Gaussian-product subdivision
#    note: igl.loop only handles 3D subdivs -> split the 9D mesh into three 3D ones
for _ in range(4):
    qq1, f = igl.loop(qq1, cm.f)
    qq2, f = igl.loop(qq2, cm.f)
    qlin, cm.f = igl.loop(qlin, cm.f)

    
# 4. transform back to 3D
cm.v = np.zeros((len(qlin),3))
for i, ql in enumerate(qlin):
    icov = [qq1[i],
            [qq1[i][1], qq2[i][0], qq2[i][1]],
            [qq1[i][2], qq2[i][1], qq2[i][2]]]
    cm.v[i] = np.linalg.inv(icov) @ ql


#------------------------------------------------------------------------------------------

# save meshes to .off file

path = os.path.join(root_folder, "data", "tweety-loop.off")
igl.write_triangle_mesh(path, lmesh.v, lmesh.f)
print("Wrote Loop mesh to", path)

path = os.path.join(root_folder, "data", "tweety-gps.off")
igl.write_triangle_mesh(path, cm.v, cm.f)
print("Wrote GPS mesh to", path)

