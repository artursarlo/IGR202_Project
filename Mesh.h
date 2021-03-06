// --------------------------------------------------------------------------
// Copyright(C) 2009-2015
// Tamy Boubekeur
//
// All rights reserved.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (http://www.gnu.org/licenses/gpl.txt)
// for more details.
// --------------------------------------------------------------------------

#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include "Vec3.h"
#include <algorithm>

class Edge {
 public:
  inline Edge () {
    edge_vertex[0] = edge_vertex[1] = 0;
    edge_mid_point = 0;
    used = true;
  }
  inline Edge (const Edge & e) {
    edge_vertex[0] = e.edge_vertex[0];
    edge_vertex[1] = e.edge_vertex[1];
    t = e.t;
    edge_mid_point = e.edge_mid_point;
    used = e.used;
  }
  inline Edge (unsigned int va, unsigned int vb, unsigned int vc =0, bool u =true) {
    edge_vertex[0] = va;
    edge_vertex[1] = vb;
    edge_mid_point = vc;
    used = u;
  }
  inline Edge & operator= (const Edge & e) {
    edge_vertex[0] = e.edge_vertex[0];
    edge_vertex[1] = e.edge_vertex[1];
    t = e.t;
    edge_mid_point = e.edge_mid_point;
    used = e.used;
    return (*this);
  }
  inline bool operator == (const Edge & e) const {
    return(((edge_vertex[0] == e.edge_vertex[0]) && (edge_vertex[1] == e.edge_vertex[1])) ||
           ((edge_vertex[1] == e.edge_vertex[0]) && (edge_vertex[0] == e.edge_vertex[1]))
           );
  }
  unsigned int edge_vertex[2];
  unsigned int edge_mid_point;
  std::vector<unsigned int> t;
  bool used;

};

/// A simple vertex class storing position and normal
class Vertex {
 public:
  inline Vertex () {
    used = true;
  }
  inline Vertex (const Vec3f & p_i, const Vec3f & n_i, bool u_i = true){
    p = p_i;
    n= n_i;
    used = u_i;
  }
  inline virtual ~Vertex () {
    used = true;
  }
  Vec3f p;
  Vec3f n;
  Vec3f g;
  float area;
  std::vector<unsigned int> Neighbor;
  std::vector<unsigned int> edges;
  bool used;

  inline bool operator == (const Vertex & v) const {
    return(p == v.p && n == v.n);
  }

  inline void add_neighbor(unsigned int v) {
    if (std::find(Neighbor.begin(), Neighbor.end(), v) == Neighbor.end()) {
      Neighbor.push_back(v);
    }
  }
  inline void add_edge(unsigned int e) {
    if (std::find(edges.begin(), edges.end(), e) == edges.end()) {
      edges.push_back(e);
    }
  }
};

/// A Triangle class expressed as a triplet of indices (over an external vertex list)
class Triangle {
 public:
  inline Triangle () {
    v[0] = v[1] = v[2] = 0;
    e[0] = e[1] = e[2] = 0;
    used = true;
  }
  inline Triangle (const Triangle & t) {
    v[0] = t.v[0];
    v[1] = t.v[1];
    v[2] = t.v[2];
    e[0] = t.e[0];
    e[1] = t.e[1];
    e[2] = t.e[2];
    used = t.used;
  }
  inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2,
                   unsigned int e0=0, unsigned int e1=0, unsigned int e2=0,
                   bool u= true) {
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
    used = u;
  }
  inline virtual ~Triangle () {}
  inline Triangle & operator= (const Triangle & t) {
    v[0] = t.v[0];
    v[1] = t.v[1];
    v[2] = t.v[2];
    e[0] = t.e[0];
    e[1] = t.e[1];
    e[2] = t.e[2];
    used = t.used;
    return (*this);
  }

  unsigned int v[3];
  unsigned int e[3];
  bool used;
};

/// A Mesh class, storing a list of vertices and a list of triangles indexed over it.
class Mesh {
 public:
	std::vector<Vertex> V;
	std::vector<Triangle> T;
  std::vector<Edge> E;

  inline Mesh& operator= (const Mesh & M) {
    T = M.T;
    V = M.V;
    E = M.E;
    return (*this);
  };
  inline unsigned int add_vertex(Vertex v) {
    std::vector<Vertex>::iterator i = std::find(V.begin(), V.end(), v);
    if (i == V.end()) {
      V.push_back(v);
      return V.size() -1;
    }
    else{
      return i - V.begin();
    }
  }

  /// Loads the mesh from a <file>.off
	void loadOFF (const std::string & filename);

  /// Compute smooth per-vertex normals
  void recomputeNormals ();

  /// scale to the unit cube and center at original
  void centerAndScaleToUnit ();

  void calculate_Voronoi_areas();
  void recomputeNeighbors ();

  bool is_obtuse_triangle(Vec3f p1, Vec3f p2, Vec3f p3);
  void recomputeEdges ();

  float area_triangle(Vec3f p1, Vec3f p2, Vec3f p3);

  void do_tangential_smoothing();

  float zero_step();

  void first_step(float l);

  void second_step(float l);

  void third_step ();
};
