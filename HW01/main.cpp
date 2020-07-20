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
using Edge = Surface_Mesh::SurfaceMesh::Edge;

int main(int argc, char *argv[]) {
	Surface_Mesh::SurfaceMesh mesh;
	mesh.read(argv[1]);
	dbg(argv[1]);
	dbg(mesh.n_vertices());
	dbg(mesh.n_faces());
	dbg(mesh.n_edges());

	std::vector<Vertex> path;
	std::vector<Edge> edges;

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

	Surface_Mesh::SurfaceMesh plane;
	plane.read(argv[4]);
	double mst_val = xry_mesh::getMST(plane, edges);
	nodes.clear();
	line.clear();
	for (const auto &p : plane.points()) {
		nodes.push_back(p[0]);
		nodes.push_back(p[1]);
		nodes.push_back(p[2]);
	}
	for (const auto &e : edges) {
		line.push_back(plane.vertex(e, 0).idx());
		line.push_back(plane.vertex(e, 1).idx());
	}
	std::ofstream ofs2("mst.vtk");
	line2vtk(ofs2, nodes.data(), nodes.size() / 3, line.data(), line.size() / 2);
	dbg(edges);
	dbg(mst_val);
	return 0;
}

