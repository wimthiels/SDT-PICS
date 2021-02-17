#!/usr/bin/env python-mpacts

import numpy as np
import math
import sys

"""
uses nodes to create a convex hull using gift wrapping algorithm
@author: micsla
"""


def giftwrap(cloud, ignore_half=False, quantile_gut=0.5, divide_remove=1, seedlen=50, debug=False):
    """filter cloud input and build a complete hull with progressive gift wrapping, increase params to be faster
    ignore_half is about twice as fast, increased quantile ]0-1] and divide_remove [1-100[ will reduce cloud"""
    print('\nconvex hull: cloud original size', len(cloud))
    cloud = remove_inner(cloud, quantile_gut)
    cloud = cloud[0::int(divide_remove)]  # Simplify by skipping points
    if len(cloud) < seedlen:
        sys.exit("Cloud too small or reduced too much")
    if debug:
        print(cloud)
    cloud0 = cloud  # backup of entire cloud
    print('convex hull: cloud filtered size', len(cloud))
    xs, ys, zs = list(zip(*cloud))
    mxs, mys, mzs = np.mean(xs), np.mean(ys), np.mean(zs)  # Find center of cloud
    xpriority = sorted(range(len(xs)), key=lambda k: xs[k])  # We want the lowest X points for triangle seed
    lowxoptions = xpriority[1:seedlen]
    x1 = x2 = x3 = xpriority[0]
    for i in lowxoptions:
        for j in lowxoptions:  # Look for third valid triangle point
            if x1 == x3 and i != j and isborder(cloud[x1], cloud[i], cloud[j], cloud):
                x2 = i
                x3 = j
    if x1 == x2 or x1 == x3 or x2 == x3:
        print('Seed', x1, x2, x3)
        sys.exit("Failed to make first triangle")
    nodes = np.array([cloud[x1], cloud[x2], cloud[x3]])  # All nodes of the mesh
    vectornode1 = np.array([nodes[0][0] - mxs, nodes[0][1] - mys, nodes[0][2] - mzs])  # Vector to the first node
    normaltri1 = normal(*nodes)  # Normal of the first triangle plane
    if np.sign(np.dot(vectornode1, normaltri1)) < 0:  # Make sure normal is facing outwards
        nodes = np.array([cloud[x2], cloud[x1], cloud[x3]])
        normaltri1 = normal(*nodes)  # Normal of the first triangle plane
    triangles = np.array([[0, 1, 2]])  # Triangles that refer to node indices
    connections = np.array([[np.array([0, 1]), 0, normaltri1], [np.array([0, 2]), 0, normaltri1],
                            [np.array([1, 2]), 0, normaltri1]])  # Open connections: sorted node indices, first
    # triangle index and its normal vector. Need a second triangle to be closed and deleted.
    grave = np.array([-1, -1])  # To blacklist connections completely
    while len(connections) > 0:
        nodes, triangles, connections, cloud, grave = \
            solve_side(nodes, triangles, connections, cloud, grave, [mxs, mys, mzs], ignore_half, debug)
        last_tri = triangles[-1]
        if len(triangles) % 50 == 0:
            print("Number of triangles: {}. Latest triangle: {}.".format(len(triangles), last_tri))
        if not isborder(nodes[last_tri[0]], nodes[last_tri[1]], nodes[last_tri[2]], cloud0):
            # plot(nodes, triangles, cloud)
            print('f')  # sys.exit("Internal triangle, convexhull messed up again")
    # plot(nodes, triangles, cloud)
    return nodes, triangles, np.array([mxs, mys, mzs])  # This is the requided information for a mesh output, (wth: added centroid needed for scaling)


def solve_side(nodes, triangles, connections, cloud, grave, means, ignore_half, debug):
    """makes and adds one new triangle from the provided mesh and open side indices"""
    con = connections[0]  # We will process only the first connection
    # con consists of 0) side node indices, 1) connected triangle index, 2) triangle normal vector
    index0 = arrayindex(triangles[con[1]], con[0][0])[0]
    index1 = arrayindex(triangles[con[1]], con[0][1])[0]
    indexp = 3 - index1 - index0  # Index 0, 1, 2 always add up to 3
    tplane = np.array([nodes[con[0][1]], nodes[con[0][1]], nodes[con[0][1]]])
    tplane[index1] = nodes[con[0][0]]  # By swapping these elements, the normal is guaranteed pointing out again
    pa = nodes[triangles[con[1]][0]]
    pb = nodes[triangles[con[1]][1]]
    pc = nodes[triangles[con[1]][2]]
    sameside = sides(tplane[index0], tplane[index1], nodes[triangles[con[1]][indexp]], means, cloud)

    angles = np.full(len(cloud), 5.)  # Initialized on something higher than pi, also i like 5
    return_dict = {}
    requested = [i for i in range(len(cloud)) if not ignore_half or sameside[i]]  # ignore_half lets us filter half
    for i in requested:
        collect(cloud, i, tplane, indexp, con, pa, pb, pc, return_dict)
    for key in return_dict.keys():
        angles[key] = return_dict[key]

    j = np.argmin(angles)  # Lowest angle is wanted
    nodes, index = add_node(nodes, cloud[j])  # Add the new node if needed
    while len(arrayindex(grave, np.sort([index, con[0][0]]))) > 0 \
            or len(arrayindex(grave, np.sort([index, con[0][1]]))) > 0 \
            or np.isnan(angles[j]):  # If invalid choose the next lowest angle
        if angles[j] == 5.:  # This should not happen
            sys.exit("no valid angle was found")
        angles[j] = 5.
        j = np.argmin(angles)  # Lowest angle is wanted
        if debug:
            print(angles)
        nodes, index = add_node(nodes, cloud[j])  # Add the new node if needed
    if con[0][0] == con[0][1] or con[0][0] == index or index == con[0][1]:
        sys.exit("impossible triangle")
    chosen_triangle = np.array([con[0][1], con[0][1], con[0][1]])
    chosen_triangle[index1] = con[0][0]
    chosen_triangle[indexp] = index
    if debug:
        print('new tri', chosen_triangle)
    triangles = np.vstack((triangles, chosen_triangle))  # Add the triangle
    tplane[indexp] = nodes[index]
    trinormal = normal(*tplane)

    for ni in triangles[-1]:  # Check nodes of the last triangle
        if node_complete(ni, triangles):  # If the node is surrounded right, it can be deleted from cloud db
            cloud = np.array([cloud[i] for i in range(len(cloud)) if not np.array_equal(cloud[i], nodes[ni])])
            if debug:
                print('point deleted from cloud')
    grave = np.vstack((grave, con[0]))  # Blacklist connection
    newconnections = np.array([connections[c] for c in range(len(connections)) if c > 0])  # Clear the processed one

    for s in range(2):
        newside = np.sort([con[0][s], index])  # Potential new open side of the newest triangle
        aindex = arrayindex(newconnections[:, 0], newside)
        if len(aindex) == 0:  # If it was not in connections; get it in connections
            newconnections = np.vstack((newconnections, np.array([[newside, len(triangles) - 1, trinormal]])))
            if debug:
                print('border con added', newconnections[-1])
        else:  # But if it was there waiting for a match, then get it out
            if debug:
                print('border con deleted', newconnections[aindex[0]])
            grave = np.vstack((grave, newconnections[aindex[0]][0]))  # Blacklist connection
            newconnections = np.array([newconnections[c] for c in range(len(newconnections)) if c != aindex[0]])

    if debug:
        print("Number of nodes: {}, number of triangles: {}".format(len(nodes), len(triangles)))
        print(newconnections)
        print(grave)
    return nodes, triangles, newconnections, cloud, grave


def arrayindex(arr, sub):
    """returns list with indices where the subarray is present"""
    return [x for x in range(len(arr)) if np.array_equal(arr[x], sub)]


def sides(pA, pB, pC, means, cloud):
    """returns which points of the cloud are on the opposite side of (pA,pB,means) compared to pC"""
    cloudsize = np.shape(cloud)[0]
    cloudm = np.c_[cloud, np.ones(cloudsize)]  # Add column of 1
    divplane = plane(pA, pB, means)
    csides = np.round(np.dot(cloudm, divplane), 12)
    counterside = np.dot(np.append(pC, 1), divplane)
    return np.sign(csides) != np.sign(counterside)


def isborder(pa, pb, pc, cloud):
    """returns true if all points of the cloud are on one side, or inside the plane defined by the 3 given points"""
    cloudsize = np.shape(cloud)[0]
    cloudm = np.c_[cloud, np.ones(cloudsize)]  # Add column of 1
    csides = np.round(np.dot(cloudm, plane(pa, pb, pc)), 12)  # sometimes float error would turn opposite sign
    sign = np.sign(csides)
    gsign = np.sign(np.sum(sign))  # Global sign
    if any(np.sign(s) == -gsign for s in csides):  # If any on the other side
        return False
    return True


def plane(pa, pb, pc):
    """returns the parameters for a plane through these points : ax + by + cz + d = 0"""
    v1 = np.subtract(pc, pa)
    v2 = np.subtract(pb, pa)
    cp = np.cross(v1, v2)  # normal
    a, b, c = cp
    d = -np.dot(cp, pc)
    return np.array([a, b, c, d])


def get_angle(cloud, i, tplane, indexp, con, pa, pb, pc):
    """gets some data and calculates the angle of the plane with point i from cloud"""
    p = cloud[i]
    if np.array_equal(p, pa) or np.array_equal(p, pb) or np.array_equal(p, pc):  # Invalid filter
        return np.array([i, 5.])
    else:
        tplane[indexp] = p
        newnormal = normal(*tplane)
        angle = math.acos(np.dot(con[2], newnormal) / (np.linalg.norm(con[2]) * np.linalg.norm(newnormal)))
        return np.array([i, angle])


def collect(cloud, i, tplane, indexp, con, pa, pb, pc, return_dict):
    """multiprocessing collector of angle stuff"""
    arangle = get_angle(cloud, i, tplane, indexp, con, pa, pb, pc)
    index, angle = arangle
    return_dict[int(index)] = angle


def add_node(nodes, node):
    """adds a new node to the array if not already present and also gives the index"""
    for i in range(len(nodes)):
        if np.array_equal(nodes[i], node):
            return nodes, i
    return np.vstack((nodes, node)), len(nodes)


def node_complete(index, triangles):
    """returns wheter a node is complete by checking if all surrounding triangles connect"""
    trinodes = np.array([index])  # Collection initialized by the point of interest
    for t in triangles:
        if index in t:
            trinodes = np.append(trinodes, t)
    trinodes = [i for i in trinodes if i != index]  # Filter out all occurences of point of interest
    unique, counts = np.unique(trinodes, return_counts=True)
    for c in counts:  # Every connected point needs to be present in 2 triangles for the circle to be closed
        if c % 2 == 1:
            return False
    return True


def normal(a, b, c):
    """returns normal vector of a plane through these points"""
    return np.cross(b - a, c - a)


def remove_inner(cloud1, quantile=0.5):
    """removes a quantile (higher = del more) of the points that are on the inside of the cloud through PCA scores"""
    cloud = cloud1.copy()
    means = np.mean(cloud1, axis=0)
    cloud -= means  # Normalize
    cov_matrix = np.cov(cloud.T)
    values, vectors = np.linalg.eig(cov_matrix)  # Principal component analysis
    dist = np.array([p.T @ vectors @ np.diag(values) @ vectors.T @ p for p in cloud])  # Ellips distance score
    meandist = np.quantile(dist, quantile)
    return np.array([cloud1[i] for i in range(len(cloud)) if dist[i] >= meandist])


# from esimutils.render import plot
# import time
# import pickle
# start_time = time.time()
# np.random.seed(13)
# nodd = np.array([np.random.uniform(-1e-1, 2e-1, 3) for i in range(1000)])
# nodd = pickle.load(open("nodes.p", "rb")) * 1e6
# print(giftwrap(nodd, True, 0.3, 4, 50, False))
# print("--- %s seconds ---" % (time.time() - start_time))
