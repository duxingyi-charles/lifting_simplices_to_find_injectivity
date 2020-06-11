#!/usr/bin/env python3

import sys
import meshio
import numpy as np

def read_handles(filename):
	"""read handles.txt into a numpy array"""
	with open(filename) as file:
		hdls = [int(line) for line in file.readlines() if line.strip()]
	return np.array(hdls)

def write_TLC_input(filename, restVert, initVert, cells, handles):
	"""write data (numpy array) as TLC program input format"""
	with open(filename, 'w') as file:
		# rest vertices
		nV, nDim = restVert.shape
		file.write(f'{nV} {nDim}\n')
		np.savetxt(file, restVert)
		# init vertices
		nV, nDim = initVert.shape
		file.write(f'{nV} {nDim}\n')
		np.savetxt(file, initVert)
		# simplex cells
		nF, simplex_size = cells.shape
		file.write(f'{nF} {simplex_size}\n')
		np.savetxt(file, cells, "%d")
		# handles
		nhdls = handles.size
		file.write(f'{nhdls}\n')
		np.savetxt(file, handles, "%d")
	

def main():
	if len(sys.argv) != 5:
		print(f'Usage: {sys.argv[0]} [restFile] [initFile] [handleFile] [outFile]')
		exit()

	restFile, initFile, handleFile, outFile = sys.argv[1:]

	# read data
	restMesh = meshio.read(restFile)
	initMesh = meshio.read(initFile)
	hdls     = read_handles(handleFile)

	# process data
	restV = restMesh.points
	
	initV = initMesh.points

	cells = restMesh.cells[0][1]

	assert np.array_equal(initMesh.cells[0][1], cells)
	
	# write data
	write_TLC_input(outFile, restV, initV, cells, hdls)



if __name__ == "__main__":
	main()

