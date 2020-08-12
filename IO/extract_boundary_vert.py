#!/usr/bin/env python3

import sys
import meshio
import numpy as np

def extract_boundary_vert_2D(cells):
	"""extract boundary vertex indices from triangle cells"""
	edge_count = {}

	# count edges
	for cell in cells:
		cell.sort()
		i,j,k = cell
		edges = [(i,j),(j,k),(i,k)]
		for e in edges:
			if e in edge_count:
				edge_count[e] += 1
			else:
				edge_count[e] = 1

	# collect boundary vertices
	bndry_vert = set()
	for e,cnt in edge_count.items():
		if cnt == 1:
			for v in e:
				bndry_vert.add(v)

	return np.array(list(bndry_vert))

def extract_boundary_vert_3D(cells):
	"""extract boundary vertex indices from tet cells"""
	face_count = {}

	# count faces
	for cell in cells:
		cell.sort()
		i,j,k,l = cell
		faces = [(i,j,k),(i,j,l),(i,k,l),(j,k,l)]
		for f in faces:
			if f in face_count:
				face_count[f] += 1
			else:
				face_count[f] = 1

	# collect boundary vertices
	bndry_vert = set()
	for f,cnt in face_count.items():
		if cnt == 1:
			for v in f:
				bndry_vert.add(v)
			
	return np.array(list(bndry_vert))

def write_handles(filename, handles):
	"""write handle data to file"""
	with open(filename, 'w') as file:
		np.savetxt(file, handles, "%d")


def main():
	if len(sys.argv) != 3:
		print(f'Usage: {sys.argv[0]} [inputMeshFile] [outputHandleFile]')
		exit()

	meshFile, handleFile = sys.argv[1:]

	# read in mesh
	mesh = meshio.read(meshFile)

	verts = mesh.points
	cells = mesh.cells[0][1]

	# extract boundary vertex indices
	nF, simplex_size = cells.shape

	if simplex_size == 3:
		handles = extract_boundary_vert_2D(cells)
	elif simplex_size == 4:
		handles = extract_boundary_vert_3D(cells)
	else:
		handles = np.array([])

	# write data
	write_handles(handleFile, handles)


if __name__ == '__main__':
	main()

