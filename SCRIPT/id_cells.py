#!/usr/bin/env python-mpacts

import os
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial.transform import Rotation as R

"""
build database Lineager savefiles and use positional data to name cells
@author: micsla
"""
data_folder = Path(os.getcwd())  / "embryosim" / "lineage"  #wth
#data_folder = os.path.abspath(os.getcwd() + "/lineage")
# data_folder = "C:/Users/mcmic/Documents/Kitematic/mpactsc/home/docker/shared/embryosim/lineage"


def identify(coords):
    """gets cell nucleus locations and hopefully returns their correct names"""
    sets = fetch(coords.shape[0])
    options = []
    for si in sets:
        er, al = aligned_error(coords, np.array(si[['x', 'y', 'z']]))
        snames = np.array(si[['cell']])
        if len(options) == 0 or er < options[0]:
            options[:] = [er, [snames[al[i]][0] for i in range(len(al))]]
    return options


def fetch(nc):
    """get all sets with nc cells"""
    path_list = Path(data_folder).glob('*aligned.gr')
    datasets = []
    for path in path_list:  # We have several datasets present
        nuclei = pd.read_csv(str(path), sep="\t")
        grouped = nuclei[['cell', 'time']].groupby('time').count()
        times = grouped.index[grouped.cell == nc]  # Time points that have nc cells
        for time in times:
            nucl = nuclei[nuclei['time'] == time]  # Slice one time
            if nucl.shape[0] == nc:
                nucl = nucl.drop(columns=["radius", "self", "parent", "child1", "child2"])
                datasets.append(nucl.copy())  # Assume we will get a good dataset
                nucl["y"] = -nucl["y"]
                datasets.append(nucl.copy())  # Assume we will get a flipped dataset
    return datasets


def aligned_error(c1_, c2_):
    """calculates the mimimal error between these clouds after aligning and rotating"""
    # mtx1, mtx2, disparity = procrustes(c1_, c2_)
    # return error(mtx1, mtx2)
    c1 = c1_.copy() - np.mean(c1_, axis=0)  # Center
    c2 = c2_.copy() - np.mean(c2_, axis=0)
    size1 = np.sum([np.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2) for p in c1])  # Scale
    size2 = np.sum([np.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2) for p in c2])
    c2 = (size1 / size2) * c2
    val1, vec1 = np.linalg.eig(np.cov(c1.T))  # PCA
    val2, vec2 = np.linalg.eig(np.cov(c2.T))
    dir1 = vec1[np.argmax(val1)]  # First vector
    dir2 = vec2[np.argmax(val2)]

    best_score_align = []
    c2al = rotate_align(c2, dir1, dir2)  # Assume front
    best_score_align = rotate_loop(best_score_align, c1, c2al, dir1)
    c2al = rotate_align(c2, dir1, -dir2)  # Assume back
    best_score_align = rotate_loop(best_score_align, c1, c2al, dir1)
    return best_score_align


def rotate_align(cloud, dir1, dir2):
    """rotate cloud from dir2 to dir1"""
    a, b = (dir2 / np.linalg.norm(dir2)).reshape(3), (dir1 / np.linalg.norm(dir1)).reshape(3)
    v = np.cross(a, b)  # Line to rotate around
    c = np.dot(a, b)  # Angle
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rota = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    cloud = (rota @ cloud.T).T  # c2 now aligned to c1
    return cloud


def rotate_loop(best_score_align, c1, c2, normalr):
    """check all rotations of c2 around normal for align error"""
    for ang in np.linspace(0, 2 * np.pi, 150):
        rot = R.from_rotvec(ang * np.array(normalr))
        c2r = rot.apply(c2)
        er, al = error(c1, c2r)
        if len(best_score_align) == 0 or er < best_score_align[0]:
            best_score_align[:] = [er, al]
    return best_score_align


def error(c1, c2):
    """calculates the error between these clouds"""
    err = 0
    used = [False] * len(c2)
    mapping = list(range(len(c1)))
    options_ = [10 * np.max(c1)] * len(c2)
    for pi in range(len(c1)):
        options = options_.copy()
        for pj in range(len(c2)):  # To find the closest point
            if not used[pj]:
                options[pj] = (c1[pi][0] - c2[pj][0]) ** 2 + (c1[pi][1] - c2[pj][1]) ** 2 + (c1[pi][2] - c2[pj][2]) ** 2
        si = np.argmin(options)
        err += options[int(si)]
        used[int(si)] = True
        mapping[pi] = si  # Mapping has indices of cloud 1 and as value their reference index in set 2
    return err, mapping


# print(identify(np.array([[0, 1, 0], [1, 0, 0], [0, 1.125, 1], [2, 1, 1], [-1, 1, 0], [-1, 0, 1], [-1.8812, 1, 1]])))
