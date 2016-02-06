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

class Edge {
 public:
  inline Edge () {
    edge_vertex[0] = edge_vertex[1] = 0;
    edge_mid_point = 0;
  }
  inline Edge (const Edge & e) {
    edge_vertex[0] = e.edge_vertex[0];
    edge_vertex[1] = e.edge_vertex[1];
    edge_mid_point = e.edge_mid_point;
  }
  inline Edge (unsigned int va, unsigned int vb, unsigned int vc =0) {
    edge_vertex[0] = va;
    edge_vertex[1] = vb;
    edge_mid_point = vc;
  }
  inline Edge & operator= (const Edge & e) {
    edge_vertex[0] = e.edge_vertex[0];
    edge_vertex[1] = e.edge_vertex[1];
    edge_mid_point = e.edge_mid_point;
    return (*this);
  }
  inline bool operator == (const Edge & e) const {
    return(((edge_vertex[0] == e.edge_vertex[0]) && (edge_vertex[1] == e.edge_vertex[1])) ||
           ((edge_vertex[1] == e.edge_vertex[0]) && (edge_vertex[0] == e.edge_vertex[1]))
           );
  }
  unsigned int edge_vertex[2];
  unsigned int edge_mid_point;
};

/// A simple vertex class storing position and normal
class Vertex {
 public:
  inline Vertex () {}
  inline Vertex (const Vec3f & p, const Vec3f & n) : p (p), n (n) {}
  inline virtual ~Vertex () {}
  Vec3f p;
  Vec3f n;
  std::vector<unsigned int> Neighbor;

  inline bool operator == (const Vertex & v) const {
    return(p == v.p && n == v.n);
  }

  inline void add_neighbor(unsigned int v) {
    if (std::find(Neighbor.begin(), Neighbor.end(), v) == Neighbor.end()) {
      Neighbor.push_back(v);
    }
  }
};

/// A Triangle class expressed as a triplet of indices (over an external vertex list)
class Triangle {
 public:
  inline Triangle () {
    v[0] = v[1] = v[2] = 0;
    e[0] = e[1] = e[2] = Edge();
    edges_calculated = false;
  }
  inline Triangle (const Triangle & t) {
    v[0] = t.v[0];
    v[1] = t.v[1];
    v[2] = t.v[2];
    e[0] = t.e[0];
    e[1] = t.e[1];
    e[2] = t.e[2];
    edges_calculated = t.edges_calculated;
  }
  inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2, Edge e0 = Edge(), Edge e1 = Edge(), Edge e2 = Edge(), bool c= false) {
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
    edges_calculated = c;
  }
  inline virtual ~Triangle () {}
  inline Triangle & operator= (const Triangle & t) {
    v[0] = t.v[0];
    v[1] = t.v[1];
    v[2] = t.v[2];
    e[0] = t.e[0];
    e[1] = t.e[1];
    e[2] = t.e[2];
    edges_calculated = t.edges_calculated;
    return (*this);
  }

  unsigned int v[3];
  Edge e[3];
  bool edges_calculated;
};

/// A Mesh class, storing a list of vertices and a list of triangles indexed over it.
class Mesh {
 public:
	std::vector<Vertex> V;
	std::vector<Triangle> T;

  inline Mesh& operator= (const Mesh & M) {
    T = M.T;
    V = M.V;
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

  void recomputeNeighbors ();

  void recomputeEdges ();

  float zero_step();

  void first_step(float l);

  void second_step(float l);
};
