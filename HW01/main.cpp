//
// Created by 徐溶延 on 2020/7/20.
//
#include <vector>
#include <fstream>

#include <SurfaceMesh/SurfaceMesh.h>
#include <SurfaceMesh/IO.h>
#include <dbg.h>

#include "../utils/SurfaceMeshUtils.h"
#include "../vtk.h"

using Vertex = Surface_Mesh::SurfaceMesh::Vertex;

int main(int argc, char *argv[]) {
	Surface_Mesh::SurfaceMesh mesh;
	mesh.read(argv[1]);
	dbg(argv[1]);
	dbg(mesh.n_vertices());
	dbg(mesh.n_faces());
	dbg(mesh.n_edges());

	std::vector<Vertex> path;

	auto shortest_distance = xry_mesh::getShortestPath<double>(mesh,
															   Vertex(std::stoi(argv[2])),
															   Vertex(std::stoi(argv[3])),
															   path);

	std::ofstream ofs("shortest_path.vtk");

	std::vector<double> nodes;
	for (const auto v : path) {
		const auto &p = mesh.position(v);
		nodes.push_back(p[0]);
		nodes.push_back(p[1]);
		nodes.push_back(p[2]);
	}

	std::vector<int> line;
	for (size_t i = 0; i < path.size() - 1; i++) {
		line.push_back(i);
		line.push_back(i + 1);
	}

	line2vtk(ofs, nodes.data(), nodes.size() / 3, line.data(), line.size() / 2);
	dbg(shortest_distance);
	dbg(path);
	return 0;
}

