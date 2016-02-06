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
#include <vector>
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

    V[T[i].v[0]].add_neighbor(T[i].v[1]);
    V[T[i].v[0]].add_neighbor(T[i].v[2]);
    V[T[i].v[1]].add_neighbor(T[i].v[0]);
    V[T[i].v[1]].add_neighbor(T[i].v[2]);
    V[T[i].v[2]].add_neighbor(T[i].v[0]);
    V[T[i].v[2]].add_neighbor(T[i].v[1]);
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

void Mesh::recomputeNeighbors () {
  for (unsigned int i = 0; i < V.size (); i++)
    V[i].Neighbor.resize(0);
  for (unsigned int i = 0; i < T.size (); i++) {
    V[T[i].v[0]].add_neighbor(T[i].v[1]);
    V[T[i].v[0]].add_neighbor(T[i].v[2]);
    V[T[i].v[1]].add_neighbor(T[i].v[0]);
    V[T[i].v[1]].add_neighbor(T[i].v[2]);
    V[T[i].v[2]].add_neighbor(T[i].v[0]);
    V[T[i].v[2]].add_neighbor(T[i].v[1]);
  }
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

  std::vector<ereaseable_edge> edges_2b_erased;
  ereaseable_edge e;

  bool apply_edge_collapse, edge_2b_killed;
  std::vector<unsigned int>::iterator ia, ib, ic, id;
  unsigned int *va;
  unsigned int *vb;
  unsigned int *vc;
  unsigned int *vd;

  for (unsigned int i = 0; i < T.size(); i++){
    // std::cerr << "Analyzing triangle "  << i << " of :" << T.size() << std::endl;

    float dist0 = dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    float dist1 = dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
    float dist2 = dist (V[T[i].v[2]].p, V[T[i].v[0]].p);

    apply_edge_collapse = false;

    if (dist0 < l){
      e.edge_vertexes[0] = T[i].v[0];
      e.edge_vertexes[1] = T[i].v[1];

      apply_edge_collapse = true;
    }
    else if (dist1 < l){
      e.edge_vertexes[0] = T[i].v[1];
      e.edge_vertexes[1] = T[i].v[2];

      apply_edge_collapse = true;
    }
    else if (dist2 < l){
      e.edge_vertexes[0] = T[i].v[0];
      e.edge_vertexes[1] = T[i].v[2];

      apply_edge_collapse = true;
    }

    if(apply_edge_collapse){
      edge_2b_killed = true;

      for (unsigned int k = 0; k < edges_2b_erased.size(); k++){
        ia = std::find(V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.begin(),
                       V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.end(),
                       e.edge_vertexes[0]
                       );
        ib = std::find(V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.begin(),
                       V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.end(),
                       e.edge_vertexes[1]
                       );
        ic = std::find(V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.begin(),
                       V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.end(),
                       e.edge_vertexes[0]
                       );
        id = std::find(V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.begin(),
                       V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.end(),
                       e.edge_vertexes[1]
                       );
        va = std::find(std::begin(edges_2b_erased[k].edge_vertexes), std::end(edges_2b_erased[k].edge_vertexes), e.edge_vertexes[0]);
        vb = std::find(std::begin(edges_2b_erased[k].edge_vertexes), std::end(edges_2b_erased[k].edge_vertexes), e.edge_vertexes[1]);

        if ((ia != V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.end()) ||
            (ib != V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.end()) ||
            (ic != V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.end()) ||
            (id != V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.end())){
          edge_2b_killed = false;
        }

        if ((va != std::end(edges_2b_erased[k].edge_vertexes)) ||
            (vb != std::end(edges_2b_erased[k].edge_vertexes))){
          edge_2b_killed = false;
        }
      }

      if(edge_2b_killed){


        V.resize(V.size() +1);

        V[V.size() -1] = Vertex((V[e.edge_vertexes[0]].p + V[e.edge_vertexes[1]].p) *0.5f,
                                (V[e.edge_vertexes[0]].n +V[e.edge_vertexes[1]].n) *0.5f
                                );

        e.vertex_edge_mid_point = V.size() -1;

        edges_2b_erased.push_back(e);

      }
    }
  }

  for (unsigned int j = 0; j < T.size(); j++){
    for (unsigned int k = 0; k < edges_2b_erased.size(); k++){
      va = std::find(std::begin(T[j].v), std::end(T[j].v), edges_2b_erased[k].edge_vertexes[0]);
      vb = std::find(std::begin(T[j].v), std::end(T[j].v), edges_2b_erased[k].edge_vertexes[1]);
      if ((va != std::end(T[j].v)) && (vb != std::end(T[j].v))){
        triangles_2b_erased.push_back(j);
        // std::cerr << "Triangles 2b erased size: " << triangles_2b_erased.size() << std::endl;
      }
      if (va != std::end(T[j].v)){
        // T[j].v[va -std::begin(T[j].v)] = edges_2b_erased[k].vertex_edge_mid_point;
        *va = edges_2b_erased[k].vertex_edge_mid_point;
      }
      if (vb != std::end(T[j].v)){
        // T[j].v[vb -std::begin(T[j].v)] = edges_2b_erased[k].vertex_edge_mid_point;
        *vb = edges_2b_erased[k].vertex_edge_mid_point;
      }
    }
  }

  for (unsigned int k = 0; k < edges_2b_erased.size(); k++){
    std::cerr << "Erased Edge:  "  << edges_2b_erased[k].edge_vertexes[0] << "," << edges_2b_erased[k].edge_vertexes[1] << std::endl;
  }

  std::sort (triangles_2b_erased.begin(), triangles_2b_erased.end());
  for (unsigned int i = 0; i < triangles_2b_erased.size(); i++){
    T.erase(T.begin() +triangles_2b_erased[i] -i);
  }
}

//   for (unsigned int i = 0; i < T.size(); i++){
//     for (unsigned int j = 0; j < edges_black_list.size(); j++){
//       va = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[0]);
//       vb = std::find(T[i].v, T[i].v +3, edges_black_list[j].edge_vertexes[1]);

//       if ((va != T[i].v +3) && (vb != T[i].v +3)){
//         triangles_2b_erased.push_back(i);
//       }
//       else if (va != T[i].v +3){
//         T[i].v[va -T[i].v] = edges_black_list[j].vertex_edge_mid_point;
//       }
//       else if (vb != T[i].v +3){
//         T[i].v[vb -T[i].v] = edges_black_list[j].vertex_edge_mid_point;
//       }
//     }
//   }

//   for (unsigned int i = 0; i < triangles_2b_erased.size(); i++){
//     T.erase(T.begin() +triangles_2b_erased[i] -i);
//   }
// }
