#!/usr/bin/python
""" 
CREDITS:
Code inspired by the https://github.com/alsignoriello/periodic_voronoi/blob/master/tile.py 
code of GitHub user "alsignoriello".

Modifications made:
    - Assumes points have 0<x<x_len ; 0<y<y_len instead of 0<x<1, 0<y<1.
    - Since our input points have been generated forming a regular triangular 
    network, we add to each box the variation over the initial positions, which 
    are provided by "mod" matrix.
    
"""
import numpy as np


def tile_points(points, N, Nx, Ny, mod):
    #x_len = np.max(points[:,0])
    #y_len = np.max(points[:,1])
    x_len = (points[1,0]-points[0,0])*Nx
    y_len = (points[2*Nx,1]-points[0,1])*Ny
    
	# 9 tiles surrounding individual coordinates
    point_tile = np.zeros((9 * N, 2)) 

	# original coordinates
    point_tile[:N] = points + mod

	# upper left 
    point_tile[N:2*N, 0] = points[:,0] - x_len + mod[:,0]
    point_tile[N:2*N, 1] = points[:,1] + y_len + mod[:,1]

	# directly above
    point_tile[2*N:3*N, 0] = points[:,0]         + mod[:,0]
    point_tile[2*N:3*N, 1] = points[:,1] + y_len + mod[:,1]

	# upper right
    point_tile[3*N:4*N, 0] = points[:,0] + x_len + mod[:,0]
    point_tile[3*N:4*N, 1] = points[:,1] + y_len + mod[:,1]

	# right
    point_tile[4*N:5*N, 0] = points[:,0] + x_len + mod[:,0]
    point_tile[4*N:5*N, 1] = points[:,1]         + mod[:,1]

	# lower right
    point_tile[5*N:6*N, 0] = points[:,0] + x_len + mod[:,0]
    point_tile[5*N:6*N, 1] = points[:,1] - y_len + mod[:,1]

	# under
    point_tile[6*N:7*N, 0] = points[:,0]         + mod[:,0]
    point_tile[6*N:7*N, 1] = points[:,1] - y_len + mod[:,1]

	# lower left
    point_tile[7*N:8*N,0] = points[:,0] - x_len + mod[:,0]
    point_tile[7*N:8*N,1] = points[:,1] - y_len + mod[:,1]

	# left 
    point_tile[8*N:,0] = points[:,0] - x_len + mod[:,0]
    point_tile[8*N:,1] = points[:,1]         + mod[:,1]

    return point_tile, x_len, y_len



# Only for regular network (with no irregularisation)
def tile_points_reg(points, N, Nx, Ny):
    #x_len = np.max(points[:,0])
    #y_len = np.max(points[:,1])
    x_len = (points[1,0]-points[0,0])*Nx
    y_len = (points[Nx*(Ny-1)*2,1]-points[0,1])
    
	# 9 tiles surrounding individual coordinates
    point_tile = np.zeros((9 * N, 2)) 

	# original coordinates
    point_tile[:N] = points

	# upper left 
    point_tile[N:2*N, 0] = points[:,0] - x_len 
    point_tile[N:2*N, 1] = points[:,1] + y_len 

	# directly above
    point_tile[2*N:3*N, 0] = points[:,0]         
    point_tile[2*N:3*N, 1] = points[:,1] + y_len 

	# upper right
    point_tile[3*N:4*N, 0] = points[:,0] + x_len 
    point_tile[3*N:4*N, 1] = points[:,1] + y_len 

	# right
    point_tile[4*N:5*N, 0] = points[:,0] + x_len 
    point_tile[4*N:5*N, 1] = points[:,1]         

	# lower right
    point_tile[5*N:6*N, 0] = points[:,0] + x_len 
    point_tile[5*N:6*N, 1] = points[:,1] - y_len 

	# under
    point_tile[6*N:7*N, 0] = points[:,0]         
    point_tile[6*N:7*N, 1] = points[:,1] - y_len 

	# lower left
    point_tile[7*N:8*N,0] = points[:,0] - x_len
    point_tile[7*N:8*N,1] = points[:,1] - y_len

	# left 
    point_tile[8*N:,0] = points[:,0] - x_len 
    point_tile[8*N:,1] = points[:,1]         

    return point_tile