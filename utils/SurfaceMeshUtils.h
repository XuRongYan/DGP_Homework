//
// Created by 徐溶延 on 2020/5/30.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H
#define LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H

#include <vector>

#include <Eigen/Dense>
#include <SurfaceMesh/SurfaceMesh.h>

#include <dbg.h>

template<typename Scalar>
void dijkstra(Surface_Mesh::SurfaceMesh &mesh,
			  const Surface_Mesh::SurfaceMesh::Vertex v_begin) {
	auto visit = mesh.get_vertex_property<bool>("v:visited");
	auto distances = mesh.get_vertex_property<Scalar>("v:distances");
	auto pre = mesh.get_vertex_property<Surface_Mesh::SurfaceMesh::Vertex>("v:pre");

	assert(visit);
	assert(distances);
	assert(pre);

	for (const auto &v : mesh.vertices()) {
		distances[v] = std::numeric_limits<Scalar>::infinity();
	}

	distances[v_begin] = 0;

	for (size_t i = 0; i < mesh.n_vertices(); i++) {
		Surface_Mesh::SurfaceMesh::Vertex vu(-1);
		assert(!mesh.is_valid(vu));
		Scalar MIN = std::numeric_limits<Scalar>::infinity();
		for (const auto &v : mesh.vertices()) {
			if (!visit[v] && distances[v] < MIN) {
				vu = v;
				MIN = distances[v];
			}
		}
		if (!mesh.is_valid(vu)) return;
		visit[vu] = true;

		for (const auto &v : mesh.vertices(vu)) {
			if (!visit[v] && distances[v] > distances[vu] + mesh.edge_length(mesh.find_edge(vu, v))) {
				distances[v] = distances[vu] + mesh.edge_length(mesh.find_edge(vu, v));
				pre[v] = vu;
			}
		}
	}
}

namespace xry_mesh {
	/**
	 * update SurfaceMesh points' positions
	 * @tparam Scalar
	 * @param mesh
	 * @param V
	 */
	template<typename Scalar>
	void rebuild3DMesh(Surface_Mesh::SurfaceMesh &mesh, const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &V) {
		size_t idx = 0;
		for (auto &p : mesh.points()) {
			p = V.col(idx++);
		}
	}

	template<typename Scalar>
	void rebuild2DMesh(Surface_Mesh::SurfaceMesh &mesh, const Eigen::Matrix<Scalar, 2, Eigen::Dynamic> &V) {
		size_t idx = 0;
		for (auto &p : mesh.points()) {
			Eigen::Matrix<Scalar, 2, 1> p2d = V.col(idx++);
			p = Surface_Mesh::Point(p2d(0, 0), p2d(1, 0), 0);
		}
	}

	/**
	 * compute a triangle's area
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<typename Scalar>
	Scalar computeArea(const Eigen::Matrix<Scalar, 3, 1> &p1,
					   const Eigen::Matrix<Scalar, 3, 1> &p2,
					   const Eigen::Matrix<Scalar, 3, 1> &p3) {
		auto A = 0.5 * (p2 - p1).cross(p3 - p1);
		return A.norm();
	}


	/**
	 * compute a triangle's area
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<typename Scalar>
	Scalar computeAreaDet(const Eigen::Matrix<Scalar, 3, 1> &p1,
						  const Eigen::Matrix<Scalar, 3, 1> &p2,
						  const Eigen::Matrix<Scalar, 3, 1> &p3) {
		Eigen::Matrix<Scalar, 3, 3> D;
		D.col(0) << 1, 1, 1;
		D.col(1) = p2 - p1;
		D.col(2) = p3 - p1;
		return 0.5 * D.determinant();
	}

	/**
	 * compute a triangle's area
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<typename Scalar>
	Scalar computeArea2D(const Eigen::Matrix<Scalar, 2, 1> &p1,
						 const Eigen::Matrix<Scalar, 2, 1> &p2,
						 const Eigen::Matrix<Scalar, 2, 1> &p3) {
		Eigen::Matrix<Scalar, 2, 2> D;
		D.col(0) = p2 - p1;
		D.col(1) = p3 - p1;
		return 0.5 * D.determinant();
	}

	/**
	 * compute a triangle's area
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<typename Scalar>
	Scalar computeArea(const Surface_Mesh::SurfaceMesh &mesh,
					   const Surface_Mesh::SurfaceMesh::Face &f) {
		std::vector<Eigen::Matrix<Scalar, 3, 1>> points;
		for (auto v : mesh.vertices(f)) {
			points.emplace_back(mesh.position(v));
		}
		assert(points.size() == 3);
		return computeArea<Scalar>(points[0], points[1], points[2]);
	}


	/**
	 * compute a triangle's area
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<typename Scalar>
	Scalar computeAreaDet(const Surface_Mesh::SurfaceMesh &mesh,
						  const Surface_Mesh::SurfaceMesh::Face &f) {
		std::vector<Eigen::Matrix<Scalar, 3, 1>> points;
		for (auto v : mesh.vertices(f)) {
			points.emplace_back(mesh.position(v));
		}
		assert(points.size() == 3);
		return computeAreaDet<Scalar>(points[0], points[1], points[2]);
	}

	/**
	 * compute areas for each triangle of mesh
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<class Scalar>
	void computeAreas(Surface_Mesh::SurfaceMesh &mesh) {
		auto areas = mesh.add_face_property<Scalar>("f:areas");
		for (auto f : mesh.faces()) {
			areas[f] = computeArea<Scalar>(mesh, f);
		}
	}


	/**
	 * compute areas for each triangle of mesh
	 * @tparam Scalar double or float
	 * @param p1
	 * @param p2
	 * @param p3
	 * @return
	 */
	template<class Scalar>
	void computeAreasDet(Surface_Mesh::SurfaceMesh &mesh) {
		auto areas = mesh.add_face_property<Scalar>("f:areas");
		for (auto f : mesh.faces()) {
			areas[f] = computeAreaDet<Scalar>(mesh, f);
		}
	}

	/**
	 * 将SurfaceMesh中的face_property转换为std::vector形式
	 * @tparam T
	 * @param mesh
	 * @param name
	 * @return
	 */
	template<class T>
	std::vector<T> faceProperty2StdVector(Surface_Mesh::SurfaceMesh &mesh, const std::string &name) {
		std::vector<T> res;
		res.reserve(mesh.n_faces());
		auto property = mesh.get_face_property<T>(name);
		for (const auto &f : mesh.faces()) {
			res.push_back(property[f]);
		}
		return res;
	}

	/**
	 * get vertices geometry position matrix
	 * @tparam Scalar
	 * @param mesh
	 * @return
	 */
	template<typename Scalar>
	Eigen::Matrix<Scalar, 3, Eigen::Dynamic> getGeometryMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
		std::vector<Scalar> points;
		for (const auto &p : mesh.points()) {
			points.push_back(p.x());
			points.push_back(p.y());
			points.push_back(p.z());
		}
		return Eigen::Map<Eigen::Matrix<Scalar, 3, Eigen::Dynamic>>(points.data(), 3, points.size() / 3);
	}

	/**
	 * get shortest distance for given pair of vertices
	 * @tparam Scalar
	 * @param mesh
	 * @param v_begin
	 * @param v_end
	 * @return
	 */
	template<typename Scalar>
	Scalar getShortestDistance(Surface_Mesh::SurfaceMesh &mesh,
							   const Surface_Mesh::SurfaceMesh::Vertex v_begin,
							   const Surface_Mesh::SurfaceMesh::Vertex v_end) {
		auto visit = mesh.add_vertex_property<bool>("v:visited");
		auto distances = mesh.add_vertex_property<Scalar>("v:distances");
		auto pre = mesh.add_vertex_property<Surface_Mesh::SurfaceMesh::Vertex>("v:pre");

		dijkstra<Scalar>(mesh, v_begin);

		Scalar res = distances[v_end];
		mesh.remove_vertex_property(pre);
		mesh.remove_vertex_property(distances);
		mesh.remove_vertex_property(visit);
		return res;
	}

	/**
	 * get shortest path for given pair of vertices
	 * @tparam Scalar
	 * @param mesh
	 * @param v_begin
	 * @param v_end
	 * @param edges
	 * @return
	 */
	template<typename Scalar>
	Scalar getShortestPath(Surface_Mesh::SurfaceMesh &mesh,
						   const Surface_Mesh::SurfaceMesh::Vertex v_begin,
						   const Surface_Mesh::SurfaceMesh::Vertex v_end,
						   std::vector<Surface_Mesh::SurfaceMesh::Vertex> &path) {
		auto visit = mesh.add_vertex_property<bool>("v:visited");
		auto distances = mesh.add_vertex_property<Scalar>("v:distances");
		auto pre = mesh.add_vertex_property<Surface_Mesh::SurfaceMesh::Vertex>("v:pre");

		path.clear();

		dijkstra<Scalar>(mesh, v_begin);

		auto v = v_end;
		while (pre[v] != v_begin) {
			path.emplace_back(v);
			v = pre[v];
		}
		path.emplace_back(v_begin);

		std::reverse(path.begin(), path.end());

		Scalar res = distances[v_end];
		mesh.remove_vertex_property(pre);
		mesh.remove_vertex_property(distances);
		mesh.remove_vertex_property(visit);
		return res;
	}

	/**
	 * get face id matrix
	 * @param mesh
	 * @return
	 */
	Eigen::Matrix3Xi getTopoMatrix(const Surface_Mesh::SurfaceMesh &mesh);

	void ensurePositiveAreas(Surface_Mesh::SurfaceMesh &mesh);

} // namespace xry_mesh




#endif //LOCALLYINJECTIVEMAPPINGS_SURFACEMESHUTILS_H
