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

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

typedef struct{
  unsigned int edge_vertexes[2];
  unsigned int vertex_edge_mid_point;
} ereaseable_edge;

void Mesh::loadOFF (const std::string & filename) {
	ifstream in (filename.c_str ());
  if (!in)
    exit (1);
	string offString;
  unsigned int sizeV, sizeT, tmp;
  in >> offString >> sizeV >> sizeT >> tmp;
  V.resize (sizeV);
  T.resize (sizeT);
  for (unsigned int i = 0; i < sizeV; i++)
    in >> V[i].p;
  int s;
  for (unsigned int i = 0; i < sizeT; i++) {
    in >> s;
    for (unsigned int j = 0; j < 3; j++)
      in >> T[i].v[j];
  }
  in.close ();
  centerAndScaleToUnit ();
  recomputeNormals ();
}

void Mesh::recomputeNormals () {
  for (unsigned int i = 0; i < V.size (); i++)
    V[i].n = Vec3f (0.0, 0.0, 0.0);
  for (unsigned int i = 0; i < T.size (); i++) {
    Vec3f e01 = V[T[i].v[1]].p -  V[T[i].v[0]].p;
    Vec3f e02 = V[T[i].v[2]].p -  V[T[i].v[0]].p;
    Vec3f n = cross (e01, e02);
    n.normalize ();
    for (unsigned int j = 0; j < 3; j++)
      V[T[i].v[j]].n += n;
  }
  for (unsigned int i = 0; i < V.size (); i++)
    V[i].n.normalize ();
}

void Mesh::centerAndScaleToUnit () {
  Vec3f c;
  for  (unsigned int i = 0; i < V.size (); i++)
    c += V[i].p;
  c /= V.size ();
  float maxD = dist (V[0].p, c);
  for (unsigned int i = 0; i < V.size (); i++){
    float m = dist (V[i].p, c);
    if (m > maxD)
      maxD = m;
  }
  for  (unsigned int i = 0; i < V.size (); i++)
    V[i].p = (V[i].p - c) / maxD;
}

// Return the average edge lenght in the mesh
// Tere is an error associated with this calculus
// The bord edges will have a less importante weight in the calculus
// No problem, since the error is small and a more complicated implementation
// Will use much more disc space and time
float  Mesh::zero_step (){

  float l_average = 0.0f;

  for (unsigned int i = 0; i < T.size (); i++){
    l_average += dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    l_average += dist (V[T[i].v[0]].p, V[T[i].v[2]].p);
    l_average += dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
  }

  return l_average/ (T.size() *3.0f);
}

// First step of the algorithm
void Mesh::first_step (float l) {

  unsigned int T_original_size = T.size();
  std::vector<int> triangles_2b_erased;
  std::vector<int> long_edges;

  for (unsigned int i = 0; i < T_original_size; i++){
    float dist0 = dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    float dist1 = dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
    float dist2 = dist (V[T[i].v[2]].p, V[T[i].v[0]].p);

    if (dist0 > l){
      long_edges.push_back(0);
    }
    if (dist1 > l){
      long_edges.push_back(1);
    }
    if (dist2 > l){
      long_edges.push_back(2);
    }

    int v0, v1, v2, v3, v4, v5;

    switch (long_edges.size()){
    case 1:
      triangles_2b_erased.push_back(i);
      T.resize(T.size() +2);
      V.resize(V.size() +1);

      V[V.size() -1] = Vertex((V[T[i].v[long_edges[0]]].p +V[T[i].v[(long_edges[0] +1)%3]].p) *0.5f,
                              (V[T[i].v[long_edges[0]]].n +V[T[i].v[(long_edges[0] +1)%3]].n) *0.5f
                              );

      v0 = V.size() -1;
      v1 = T[i].v[(long_edges[0] +1)%3];
      v2 = T[i].v[(long_edges[0] +2)%3];
      v3 = T[i].v[long_edges[0]];

      T[T.size() -2] = Triangle(v0, v1, v2);
      T[T.size() -1] = Triangle(v0, v2, v3);
      break;
    case 2:
      triangles_2b_erased.push_back(i);
      T.resize(T.size() +3);
      V.resize(V.size() +2);

      V[V.size() -2] = Vertex((V[T[i].v[long_edges[0]]].p +V[T[i].v[(long_edges[0] +1)%3]].p) *0.5f,
                              (V[T[i].v[long_edges[0]]].n +V[T[i].v[(long_edges[0] +1)%3]].n) *0.5f
                              );
      V[V.size() -1] = Vertex((V[T[i].v[long_edges[1]]].p +V[T[i].v[(long_edges[1] +1)%3]].p) *0.5f,
                              (V[T[i].v[long_edges[1]]].n +V[T[i].v[(long_edges[1] +1)%3]].n) *0.5f
                              );

      if((long_edges[0] == 0) && (long_edges[1] == 2)){
        v0 = V.size() -1;
        v1 = T[i].v[long_edges[0]];
        v2 = V.size() -2;
        v3 = T[i].v[(long_edges[0] +1)%3];
        v4 = T[i].v[(long_edges[0] +2)%3];
      }
      else{
        v0 = V.size() -2;
        v1 = T[i].v[(long_edges[0] +1)%3];
        v2 = V.size() -1;
        v3 = T[i].v[(long_edges[0] +2)%3];
        v4 = T[i].v[long_edges[0]];
      }
      T[T.size() -3] = Triangle(v0, v1, v2);
      T[T.size() -2] = Triangle(v0, v2, v3);
      T[T.size() -1] = Triangle(v0, v3, v4);
      break;
    case 3:
      triangles_2b_erased.push_back(i);
      T.resize(T.size() +4);
      V.resize(V.size() +3);

      V[V.size() -3] = Vertex((V[T[i].v[long_edges[0]]].p +V[T[i].v[(long_edges[0] +1)%3]].p) *0.5f,
                              (V[T[i].v[long_edges[0]]].n +V[T[i].v[(long_edges[0] +1)%3]].n) *0.5f
                              );
      V[V.size() -2] = Vertex((V[T[i].v[long_edges[1]]].p +V[T[i].v[(long_edges[1] +1)%3]].p) *0.5f,
                              (V[T[i].v[long_edges[1]]].n +V[T[i].v[(long_edges[1] +1)%3]].n) *0.5f
                              );
      V[V.size() -1] = Vertex((V[T[i].v[long_edges[2]]].p +V[T[i].v[(long_edges[2] +1)%3]].p) *0.5f,
                              (V[T[i].v[long_edges[2]]].n +V[T[i].v[(long_edges[2] +1)%3]].n) *0.5f
                              );

      v0 = V.size() -3;
      v1 = T[i].v[(long_edges[0] +1)%3];
      v2 = V.size() -2;
      v3 = V.size() -1;
      v4 = T[i].v[long_edges[0]];
      v5 = T[i].v[(long_edges[0] +2)%3];

      T[T.size() -4] = Triangle(v0, v1, v2);
      T[T.size() -3] = Triangle(v0, v2, v3);
      T[T.size() -2] = Triangle(v0, v3, v4);
      T[T.size() -1] = Triangle(v2, v5, v3);
      break;
    default:
      break;
    }

    long_edges.resize(0);
  }

  for (unsigned int i = 0; i < triangles_2b_erased.size(); i++){
    T.erase(T.begin() +triangles_2b_erased[i] -i);
  }
}

// Second step of the algorithm
void Mesh::second_step (float l) {

  std::vector<int> triangles_2b_erased;
  std::vector<ereaseable_edge> edges_black_list;
  ereaseable_edge e;

  unsigned int *va;
  unsigned int *vb;
  bool black_edge_already_marked;
  bool black_edge_touching_other;

  for (unsigned int i = 0; i < T.size(); i++){
    float dist0 = dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    float dist1 = dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
    float dist2 = dist (V[T[i].v[2]].p, V[T[i].v[0]].p);

    if (dist0 < l){
      e.edge_vertexes[0] = T[i].v[0];
      e.edge_vertexes[1] = T[i].v[1];

      black_edge_already_marked = false;
      black_edge_touching_other = false;

      for (unsigned int j = 0; j < edges_black_list.size(); j++){
        va = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[0]);
        vb = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[1]);

        if ((va != T[i].v +3) && (vb != T[i].v +3)){
          black_edge_already_marked = true;
        }
        else if((va != T[i].v +3) || (vb != T[i].v +3)){
          black_edge_touching_other = true;
          e.vertex_edge_mid_point = edges_black_list[j].vertex_edge_mid_point;
        }
      }

      if(!black_edge_already_marked){
        if(!black_edge_touching_other){
          V.resize(V.size() +1);

          V[V.size() -1] = Vertex((V[T[i].v[0]].p + V[T[i].v[1]].p) *0.5f,
                                  (V[T[i].v[0]].n +V[T[i].v[1]].n) *0.5f
                                  );

          e.vertex_edge_mid_point = V.size() -1;
        }
      }
      edges_black_list.push_back(e);
    }
    else if (dist1 < l){
      e.edge_vertexes[0] = T[i].v[1];
      e.edge_vertexes[1] = T[i].v[2];

      black_edge_already_marked = false;
      black_edge_touching_other = false;

      for (unsigned int j = 0; j < edges_black_list.size(); j++){
        va = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[0]);
        vb = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[1]);

        if ((va != T[i].v +3) && (vb != T[i].v +3)){
          black_edge_already_marked = true;
        }
        else if((va != T[i].v +3) || (vb != T[i].v +3)){
          black_edge_touching_other = true;
          e.vertex_edge_mid_point = edges_black_list[j].vertex_edge_mid_point;
        }
      }

      if(!black_edge_already_marked){
        if(!black_edge_touching_other){
          V.resize(V.size() +1);

          V[V.size() -1] = Vertex((V[T[i].v[1]].p + V[T[i].v[2]].p) *0.5f,
                                  (V[T[i].v[1]].n +V[T[i].v[2]].n) *0.5f
                                  );

          e.vertex_edge_mid_point = V.size() -1;
        }
      }
      edges_black_list.push_back(e);
    }
    else if (dist2 < l){
      e.edge_vertexes[0] = T[i].v[0];
      e.edge_vertexes[1] = T[i].v[2];

      black_edge_already_marked = false;
      black_edge_touching_other = false;

      for (unsigned int j = 0; j < edges_black_list.size(); j++){
        va = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[0]);
        vb = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[1]);

        if ((va != T[i].v +3) && (vb != T[i].v +3)){
          black_edge_already_marked = true;
        }
        else if((va != T[i].v +3) || (vb != T[i].v +3)){
          black_edge_touching_other = true;
          e.vertex_edge_mid_point = edges_black_list[j].vertex_edge_mid_point;
        }
      }

      if(!black_edge_already_marked){
        if(!black_edge_touching_other){
          V.resize(V.size() +1);

          V[V.size() -1] = Vertex((V[T[i].v[0]].p + V[T[i].v[2]].p) *0.5f,
                                  (V[T[i].v[0]].n +V[T[i].v[2]].n) *0.5f
                                  );

          e.vertex_edge_mid_point = V.size() -1;
        }
      }
      edges_black_list.push_back(e);
    }
  }

  for (unsigned int i = 0; i < T.size(); i++){
    for (unsigned int j = 0; j < edges_black_list.size(); j++){
      va = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[0]);
      vb = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[1]);

      if ((va != T[i].v +3) && (vb != T[i].v +3)){
        triangles_2b_erased.push_back(i);
      }
      else if (va != T[i].v +3){
        T[i].v[va -T[i].v] = edges_black_list[j].vertex_edge_mid_point;
      }
      else if (vb != T[i].v +3){
        T[i].v[vb -T[i].v] = edges_black_list[j].vertex_edge_mid_point;
      }
    }
  }

  for (unsigned int i = 0; i < triangles_2b_erased.size(); i++){
    T.erase(T.begin() +triangles_2b_erased[i] -i);
  }
}
