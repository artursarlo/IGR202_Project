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

#if REMESH_VERBOSE
#include <ctime>
#endif

#define REMESH_FACTOR_L_PERCENTAGE 0.9f

using namespace std;

typedef struct{
  unsigned int edge_vertexes[2];
  unsigned int edge_mid_point;
} ereaseable_edge;

void Mesh::loadOFF (const std::string & filename) {
  
#if REMESH_VERBOSE
  std::clock_t begin_time, end_time;
  double elapsed_secs;
	std::cerr << "loadOFF Begin..." << std::endl;
  begin_time = clock();
  #endif

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

  #if REMESH_VERBOSE
  end_time = clock();
  elapsed_secs = double(end_time - begin_time) /CLOCKS_PER_SEC;
  std::cerr << "loadOFF End... Elapsed time (seconds): "
  << elapsed_secs << std::endl;
  #endif
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

  #if REMESH_VERBOSE
  std::clock_t begin_time, end_time;
  double elapsed_secs;

  std::cerr << "recomputeEdges Begin..." << std::endl;
  begin_time = clock();
  #endif

  std::vector<Edge>::iterator It;
  Edge new_edges[3];

  E.resize(0);
  for (unsigned int i = 0; i < V.size (); i++) {
    V[i].edges.resize(0);
  }

  for (unsigned int i = 0; i < T.size (); i++) {
    new_edges[0] = Edge (T[i].v[0], T[i].v[1]);
    new_edges[1] = Edge (T[i].v[1], T[i].v[2]);
    new_edges[2] = Edge (T[i].v[2], T[i].v[0]);

    for (unsigned int k =0; k<3; k++){
      It = std::find(E.begin(), E.end(), new_edges[k]);
      if (It == E.end()) {
        V.push_back(Vertex((V[T[i].v[k]].p +V[T[i].v[(k+1)%3]].p) *0.5f,
                           (V[T[i].v[k]].n +V[T[i].v[(k+1)%3]].n) *0.5f
                           ));
        new_edges[k].edge_mid_point =  V.size() -1;

        E.push_back(new_edges[k]);
        T[i].e[k] = E.size() -1;
      }
      else{
        T[i].e[k] = It -E.begin();
      }
    }

    V[T[i].v[0]].add_edge(T[i].e[0]);
    V[T[i].v[0]].add_edge(T[i].e[2]);
    V[T[i].v[1]].add_edge(T[i].e[0]);
    V[T[i].v[1]].add_edge(T[i].e[1]);
    V[T[i].v[2]].add_edge(T[i].e[1]);
    V[T[i].v[2]].add_edge(T[i].e[2]);

  }

  #if REMESH_VERBOSE
  end_time = clock();
  elapsed_secs = double(end_time - begin_time) /CLOCKS_PER_SEC;
  std::cerr << "recomputeEdges End... Elapsed time (seconds): "
            << elapsed_secs << std::endl;
  #endif

}


/**
 * Calculates the average edge length in the mesh.
 *
 * @return Average edge length in the mesh
 */
float Mesh::zero_step (){

  #if REMESH_VERBOSE
  std::clock_t begin_time, end_time;
  double elapsed_secs;

  std::cerr << "zero_step Begin..." << std::endl;
  begin_time = clock();
  #endif

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

  #if REMESH_VERBOSE
  end_time = clock();
  elapsed_secs = double(end_time - begin_time) /CLOCKS_PER_SEC;
  std::cerr << "zero_step End... Elapsed time (seconds): "
            << elapsed_secs << std::endl;
  #endif

  return l_average /edge_count;
}

/**
 * Applies an edge split in the mesh for all edges bigger than (4/3)*average
 * edge length.
 *
 * @param l Average edge length 
 */
void Mesh::first_step (float l) {

  #if REMESH_VERBOSE
  std::clock_t begin_time, end_time;
  double elapsed_secs;

  std::cerr << "fisrt_step Begin..." << std::endl;
  begin_time = clock();
  #endif

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

      v0 = E[T[i].e[long_edges[0]]].edge_mid_point;
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
        v0 = E[T[i].e[long_edges[1]]].edge_mid_point;
        v1 = T[i].v[long_edges[0]];
        v2 = E[T[i].e[long_edges[0]]].edge_mid_point;
        v3 = T[i].v[(long_edges[0] +1)%3];
        v4 = T[i].v[(long_edges[0] +2)%3];
      }
      else{
        v0 = E[T[i].e[long_edges[0]]].edge_mid_point;;
        v1 = T[i].v[(long_edges[0] +1)%3];
        v2 = E[T[i].e[long_edges[1]]].edge_mid_point;
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

      v0 = E[T[i].e[long_edges[0]]].edge_mid_point;
      v1 = T[i].v[(long_edges[0] +1)%3];
      v2 = E[T[i].e[long_edges[1]]].edge_mid_point;
      v3 = E[T[i].e[long_edges[2]]].edge_mid_point;
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

  #if REMESH_VERBOSE
  end_time = clock();
  elapsed_secs = double(end_time - begin_time) /CLOCKS_PER_SEC;
  std::cerr << "first_step End... Elapsed time (seconds): "
            << elapsed_secs << std::endl;
  #endif
}

/**
 * Applies an edge collapse in the mesh for all edges smaller than (4/5)*average
 * edge length.
 *
 * @param l Average edge length 
 */
void Mesh::second_step (float l) {

  #if REMESH_VERBOSE
  std::clock_t begin_time, end_time;
  double elapsed_secs;

  std::cerr << "second_step Begin..." << std::endl;
  begin_time = clock();
  #endif

  float l_floor = l *REMESH_FACTOR_L_PERCENTAGE *4.0f /5.0f;

  std::vector<int> triangles_2b_erased;
  std::vector<Edge> edges_2b_erased;
  std::vector<Edge>::iterator It;
  Edge *e;
  bool apply_edge_collapse, edge_2b_killed;

  unsigned int *va;
  unsigned int *vb;
  unsigned int active_edges;

  for (unsigned int i = 0; i < T.size(); i++){

    float dist0 = dist (V[T[i].v[0]].p, V[T[i].v[1]].p);
    float dist1 = dist (V[T[i].v[1]].p, V[T[i].v[2]].p);
    float dist2 = dist (V[T[i].v[2]].p, V[T[i].v[0]].p);

    apply_edge_collapse = false;

    if (dist0 < l_floor){
      e = &E[T[i].e[0]];
      apply_edge_collapse = true;
    }
    else if (dist1 < l_floor){
      e = &E[T[i].e[1]];
      apply_edge_collapse = true;
    }
    else if (dist2 < l_floor){
      e = &E[T[i].e[2]];
      apply_edge_collapse = true;
    }

    if(apply_edge_collapse){
      if(e->used){
        if(V[e->edge_vertex[0]].used && V[e->edge_vertex[1]].used){
          e->used = false;
          V[e->edge_vertex[0]].used = false;
          V[e->edge_vertex[1]].used = false;
        }
      }    
    }
  }

  for (unsigned int j = 0; j < T.size(); j++){
    active_edges = 0;

    if(E[T[j].e[0]].used)
      active_edges += 1;
    if(E[T[j].e[1]].used)
      active_edges += 1;
    if(E[T[j].e[2]].used)
      active_edges += 1;

    if (active_edges < 3){
      triangles_2b_erased.push_back(j);
    }
    else{
      for (unsigned int k = 0; k < 3; k++){
        if (!V[T[j].v[k]].used){
          for (unsigned int w = 0; w < V[T[j].v[k]].edges.size(); w++){
            if (!E[V[T[j].v[k]].edges[w]].used){
              T[j].v[k] = E[V[T[j].v[k]].edges[w]].edge_mid_point;
            }
          }
        }
      }
    }
  }

  std::sort (triangles_2b_erased.begin(), triangles_2b_erased.end());
  for (unsigned int i = 0; i < triangles_2b_erased.size(); i++){
    T.erase(T.begin() +triangles_2b_erased[i] -i);
  }

  #if REMESH_VERBOSE
  end_time = clock();
  elapsed_secs = double(end_time - begin_time) /CLOCKS_PER_SEC;
  std::cerr << "second_step End... Elapsed time (seconds): "
            << elapsed_secs << std::endl;
  #endif
}
