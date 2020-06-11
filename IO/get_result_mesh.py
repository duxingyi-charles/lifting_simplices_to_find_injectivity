#!/usr/bin/env python3

import sys
import meshio
import numpy as np  # version >= 1.16.0

def read_TLC_input(filename):
	"""read rest mesh, init mesh and handles from TLC input file"""
	with open(filename) as file:
		# rest vertices
		line = file.readline().split()
		nV = int(line[0])
		restV = np.loadtxt(file, dtype = np.float64, max_rows = nV)

		# initial vertices
		line = file.readline().split()
		nV = int(line[0])
		initV = np.loadtxt(file, dtype = np.float64, max_rows = nV)

		# simplex cells
		line = file.readline().split()
		nF = int(line[0])
		cells = np.loadtxt(file, dtype = int, max_rows = nF)

		# handles
		line = file.readline()
		nHdls = int(line)
		handles = np.loadtxt(file, dtype = int, max_rows = nHdls)
	return {"restV":restV, "initV":initV, "cells":cells, "handles":handles}


def read_TLC_result(filename):
	"""read TLC program results into a dictionary"""
	result_dict = {}

	with open(filename) as file:
		while True:
			head = file.readline()
			if head.strip() == "":
				break
			head = head.split()
			data_name = head[0]
			array_shape = [int(i) for i in head[1:]]

			data = np.loadtxt(file, dtype=np.float64, max_rows=1)
			if len(array_shape) == 1:
				result_dict[data_name] = data
			else:
				result_dict[data_name] = np.reshape(data, array_shape)

	return result_dict



def main():
	if len(sys.argv) != 4:
		print(f'Usage: {sys.argv[0]} [inputFile] [resultFile] [outFile]')
		exit()

	inputFile, resultFile, outFile = sys.argv[1:]

	# read TLC data
	input_dict  = read_TLC_input(inputFile)
	result_dict = read_TLC_result(resultFile)

	# compute TLC result mesh
	pts = result_dict["resV"]
	nV, nDim = pts.shape
	if nDim < 3:  # pad to 3D points
		pts = np.pad(pts,((0,0),(0,3-nDim)))

	raw_cells = input_dict["cells"]
	nF, simplex_size = raw_cells.shape
	
	if simplex_size == 3: # tri
		cells = [ ("triangle", raw_cells)]
	else: # tet
		cells = [ ("tetra", raw_cells) ]
	
	resMesh = meshio.Mesh(pts, cells)

	# write result mesh to file
	meshio.write(outFile, resMesh) 



if __name__ == '__main__':
	main()