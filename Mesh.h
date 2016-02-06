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
#include "Vec3.h"
#include <algorithm>

/// A simple vertex class storing position and normal
class Vertex {
 public:
  inline Vertex () {}
  inline Vertex (const Vec3f & p, const Vec3f & n) : p (p), n (n) {}
  inline virtual ~Vertex () {}
  Vec3f p;
  Vec3f n;
  Vec3f g;
  float area;
  std::vector<unsigned int> Neighbor;

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
  }
  inline Triangle (const Triangle & t) {
    v[0] = t.v[0];
    v[1] = t.v[1];
    v[2] = t.v[2];
  }
  inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2) {
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
  }
  inline virtual ~Triangle () {}
  inline Triangle & operator= (const Triangle & t) {
    v[0] = t.v[0];
    v[1] = t.v[1];
    v[2] = t.v[2];
    return (*this);
  }
  unsigned int v[3];
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

  /// Loads the mesh from a <file>.off
	void loadOFF (const std::string & filename);

  /// Compute smooth per-vertex normals
  void recomputeNormals ();

  /// scale to the unit cube and center at original
  void centerAndScaleToUnit ();

  void calculate_Voronoi_areas();

  bool is_obtuse_triangle(Vec3f p1, Vec3f p2, Vec3f p3);

  float area_triangle(Vec3f p1, Vec3f p2, Vec3f p3);

  void do_tangential_smoothing();

  float zero_step();

  void first_step(float l);

  void second_step(float l);
};
