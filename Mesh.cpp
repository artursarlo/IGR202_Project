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

#define REMESH_FACTOR_L_PERCENTAGE 0.9f

using namespace std;

typedef struct{
  unsigned int edge_vertexes[2];
  unsigned int edge_mid_point;
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
  recomputeEdges ();
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

void Mesh::recomputeEdges () {
  std::vector<Edge> computed_edges;
  std::vector<Edge>::iterator It;

  for (unsigned int i = 0; i < T.size (); i++) {
    if(!T[i].edges_calculated){
      T[i].e[0] = Edge (T[i].v[0], T[i].v[1]);
      T[i].e[1] = Edge (T[i].v[1], T[i].v[2]);
      T[i].e[2] = Edge (T[i].v[2], T[i].v[0]);

      for (unsigned int k =0; k<3; k++){
        It = std::find(computed_edges.begin(), computed_edges.end(), T[i].e[k]);
        if (It == computed_edges.end()) {
          V.push_back(Vertex((V[T[i].v[k]].p +V[T[i].v[(k+1)%3]].p) *0.5f,
                             (V[T[i].v[k]].n +V[T[i].v[(k+1)%3]].n) *0.5f
                             ));
          T[i].e[k].edge_mid_point =  V.size() -1;
          computed_edges.push_back(T[i].e[k]);
        }
        else{
          T[i].e[k].edge_mid_point =
            computed_edges[It -computed_edges.begin()].edge_mid_point;
        }
      }
      T[i].edges_calculated = true;
    }
  }
}


/**
 * Calculates the average edge length in the mesh.
 *
 * @return Average edge length in the mesh
 */
float Mesh::zero_step (){

  Edge e[3];
  std::vector<Edge> mesh_edges;
  std::vector<Edge>::iterator It;
  float l_average = 0.0f;
  unsigned int edge_count = 0;

  for (unsigned int i = 0; i < T.size (); i++) {
    e[0] = Edge (T[i].v[0], T[i].v[1]);
    e[1]= Edge (T[i].v[1], T[i].v[2]);
    e[2] = Edge (T[i].v[2], T[i].v[0]);

      for (unsigned int k =0; k<3; k++){
        It = std::find(mesh_edges.begin(), mesh_edges.end(), e[k]);
        if (It == mesh_edges.end()) {
          mesh_edges.push_back(e[k]);
          edge_count += 1;
          l_average += dist (V[e[k].edge_vertex[0]].p,
                             V[e[k].edge_vertex[1]].p);
        }
      }
  }

  return l_average /edge_count;
}

/**
 * Applies an edge split in the mesh for all edges bigger than (4/3)*average
 * edge length.
 *
 * @param l Average edge length 
 */
void Mesh::first_step (float l) {

  float l_roof = l *REMESH_FACTOR_L_PERCENTAGE *4.0f /3.0f;
  unsigned int T_original_size = T.size();
  std::vector<int> triangles_2b_erased;
  std::vector<int> long_edges;

  for (unsigned int i = 0; i < T_original_size; i++){
    float dist0 = dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    float dist1 = dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
    float dist2 = dist (V[T[i].v[2]].p, V[T[i].v[0]].p);

    if (dist0 > l_roof){
      long_edges.push_back(0);
    }
    if (dist1 > l_roof){
      long_edges.push_back(1);
    }
    if (dist2 > l_roof){
      long_edges.push_back(2);
    }

    int v0, v1, v2, v3, v4, v5;

    switch (long_edges.size()){
    case 1:
      triangles_2b_erased.push_back(i);
      T.resize(T.size() +2);

      v0 = T[i].e[long_edges[0]].edge_mid_point;
      v1 = T[i].v[(long_edges[0] +1)%3];
      v2 = T[i].v[(long_edges[0] +2)%3];
      v3 = T[i].v[long_edges[0]];

      T[T.size() -2] = Triangle(v0, v1, v2);
      T[T.size() -1] = Triangle(v0, v2, v3);
      break;
    case 2:
      triangles_2b_erased.push_back(i);
      T.resize(T.size() +3);

      if((long_edges[0] == 0) && (long_edges[1] == 2)){
        v0 = T[i].e[long_edges[1]].edge_mid_point;
        v1 = T[i].v[long_edges[0]];
        v2 = T[i].e[long_edges[0]].edge_mid_point;
        v3 = T[i].v[(long_edges[0] +1)%3];
        v4 = T[i].v[(long_edges[0] +2)%3];
      }
      else{
        v0 = T[i].e[long_edges[0]].edge_mid_point;;
        v1 = T[i].v[(long_edges[0] +1)%3];
        v2 = T[i].e[long_edges[1]].edge_mid_point;
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

      v0 = T[i].e[long_edges[0]].edge_mid_point;
      v1 = T[i].v[(long_edges[0] +1)%3];
      v2 = T[i].e[long_edges[1]].edge_mid_point;
      v3 = T[i].e[long_edges[2]].edge_mid_point;
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

/**
 * Applies an edge collapse in the mesh for all edges smaller than (4/5)*average
 * edge length.
 *
 * @param l Average edge length 
 */
void Mesh::second_step (float l) {

  float l_floor = l *REMESH_FACTOR_L_PERCENTAGE *4.0f /5.0f;

  std::vector<int> triangles_2b_erased;

  std::vector<ereaseable_edge> edges_2b_erased;
  ereaseable_edge e;

  bool apply_edge_collapse, edge_2b_killed;
  std::vector<unsigned int>::iterator ia, ib, ic, id;
  unsigned int *va;
  unsigned int *vb;

  for (unsigned int i = 0; i < T.size(); i++){

    float dist0 = dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    float dist1 = dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
    float dist2 = dist (V[T[i].v[2]].p, V[T[i].v[0]].p);

    apply_edge_collapse = false;

    if (dist0 < l_floor){
      e.edge_vertexes[0] = T[i].v[0];
      e.edge_vertexes[1] = T[i].v[1];

      apply_edge_collapse = true;
    }
    else if (dist1 < l_floor){
      e.edge_vertexes[0] = T[i].v[1];
      e.edge_vertexes[1] = T[i].v[2];

      apply_edge_collapse = true;
    }
    else if (dist2 < l_floor){
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

        // if ((ia != V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.end()) ||
        //     (ib != V[edges_2b_erased[k].edge_vertexes[0]].Neighbor.end()) ||
        //     (ic != V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.end()) ||
        //     (id != V[edges_2b_erased[k].edge_vertexes[1]].Neighbor.end())){
        //   edge_2b_killed = false;
        // }

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

        e.edge_mid_point = V.size() -1;

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
      }
      if (va != std::end(T[j].v)){
        *va = edges_2b_erased[k].edge_mid_point;
      }
      if (vb != std::end(T[j].v)){
        *vb = edges_2b_erased[k].edge_mid_point;
      }
    }
  }

  std::sort (triangles_2b_erased.begin(), triangles_2b_erased.end());
  for (unsigned int i = 0; i < triangles_2b_erased.size(); i++){
    T.erase(T.begin() +triangles_2b_erased[i] -i);
  }
}
