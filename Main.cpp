// ----------------------------------------------
// Informatique Graphique 3D & R�alit� Virtuelle.
// Travaux Pratiques
// Traitement G�om�trique
// Copyright (C) 2015 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <GL/glut.h>

#include "Vec3.h"
#include "Camera.h"
#include "Mesh.h"

using namespace std;

static const unsigned int DEFAULT_SCREENWIDTH = 1024;
static const unsigned int DEFAULT_SCREENHEIGHT = 768;
static const string DEFAULT_MESH_FILE ("MeshCollection/max_50k.off");

static string appTitle ("Informatique Graphique & Realite Virtuelle - Travaux Pratiques - Traitement G�om�trique");
static GLint window;
static unsigned int FPS = 0;
static bool fullScreen = false;

static Camera camera;
static Mesh mesh, o_mesh, one_mesh, two_mesh, three_mesh, four_mesh;

static float l_average;
static int remesh_counter = 0;

void printUsage () {
	std::cerr << std::endl
            << appTitle << std::endl
            << "Author: Tamy Boubekeur" << std::endl << std::endl
            << "Usage: ./main [<file.off>]" << std::endl
            << "Commands:" << std::endl
            << "------------------" << std::endl
            << " ?: Print help" << std::endl
            << " w: Toggle wireframe mode" << std::endl
            << " +: Increment last iteration step" << std::endl
            << " -: Decrement last iteration step" << std::endl
            << " <drag>+<left button>: rotate model" << std::endl
            << " <drag>+<right button>: move model" << std::endl
            << " <drag>+<middle button>: zoom" << std::endl
            << " q, <esc>: Quit" << std::endl << std::endl;

}
void initLighting () {
  GLfloat light_position1[4] = {52, 16, 50, 0};
  GLfloat light_position2[4] = {26, 48, 50, 0};
  GLfloat light_position3[4] = {-16, 52, 50, 0};
  GLfloat direction1[3] = {-52,-16,-50};
  GLfloat direction2[3] = {-26,-48,-50};
  GLfloat direction3[3] = {16,-52,-50};
  GLfloat color1[4] = {1,0, 0, 1};
  GLfloat color2[4] = {0, 1, 0, 1};
  GLfloat color3[4] = {0, 0, 1, 1};
  GLfloat ambient[4] = {0.1f, 0.1f, 0.1f, 0.5f};
  glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
  glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
  glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
  glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
  glLightfv (GL_LIGHT2, GL_POSITION, light_position2);
  glLightfv (GL_LIGHT2, GL_SPOT_DIRECTION, direction2);
  glLightfv (GL_LIGHT2, GL_DIFFUSE, color2);
  glLightfv (GL_LIGHT2, GL_SPECULAR, color2);
  glLightfv (GL_LIGHT3, GL_POSITION, light_position3);
  glLightfv (GL_LIGHT3, GL_SPOT_DIRECTION, direction3);
  glLightfv (GL_LIGHT3, GL_DIFFUSE, color3);
  glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
  glEnable (GL_LIGHT1);
  glEnable (GL_LIGHT2);
  glEnable (GL_LIGHT3);
  glDisable (GL_LIGHT0);
  glEnable (GL_LIGHTING);
}

void initMaterial () {
  GLfloat material_color[4] = {1,1,1,1.0f};
  GLfloat material_specular[4] = {0.85,0.85,0.85,1.0};
  GLfloat material_ambient[4] = {0.0,0.0,0.0,0.0};
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
  glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, 128.0f);
}

void init (const char * modelFilename) {
  glCullFace (GL_BACK);     // Specifies the faces to cull (here the ones pointing away from the camera)
  glEnable (GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
  glDepthFunc (GL_LESS); // Specify the depth test for the z-buffer
  glEnable (GL_DEPTH_TEST); // Enable the z-buffer in the rasterization
	glLineWidth (2.0); // Set the width of edges in GL_LINE polygon mode
  glClearColor (0.0f, 0.0f, 0.0f, 1.0f); // Background color
  glClearColor (0.0f, 0.0f, 0.0f, 1.0f);

	initLighting ();
	initMaterial ();
	glDisable (GL_COLOR_MATERIAL);

	mesh.loadOFF (modelFilename);

  o_mesh = mesh;
  one_mesh = mesh;
  two_mesh = mesh;
  three_mesh = mesh;
  four_mesh = mesh;

  l_average = one_mesh.zero_step();

  for (unsigned int k = 0; k < 10; k++) {
    one_mesh = four_mesh;

    one_mesh.first_step(l_average);
    one_mesh.recomputeEdges();

    two_mesh = one_mesh;
    two_mesh.second_step(l_average);
    two_mesh.recomputeEdges();
    two_mesh.recomputeNeighbors();

    three_mesh = two_mesh;
    three_mesh.third_step();
    three_mesh.recomputeEdges();
    three_mesh.recomputeNeighbors();
    three_mesh.recomputeNormals();


    four_mesh = three_mesh;
    four_mesh.calculate_Voronoi_areas();
		// The 4th step in the algorithm, but using area-based tangential smoothing
	  four_mesh.do_tangential_smoothing();
  }

  one_mesh.recomputeNormals();
  two_mesh.recomputeNormals();
  three_mesh.recomputeNormals();
  four_mesh.recomputeNormals();

  camera.resize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
}

void drawScene () {
  glBegin (GL_TRIANGLES);
  for (unsigned int i = 0; i < mesh.T.size (); i++)
    for (unsigned int j = 0; j < 3; j++) {
      const Vertex & v = mesh.V[mesh.T[i].v[j]];
      glNormal3f (v.n[0], v.n[1], v.n[2]); // Specifies current normal vertex
      glVertex3f (v.p[0], v.p[1], v.p[2]); // Emit a vertex (one triangle is emitted each time 3 vertices are emitted)
    }
  glEnd ();
}

void reshape(int w, int h) {
  camera.resize (w, h);
}

void display () {
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  camera.apply ();
  drawScene ();
  glFlush ();
  glutSwapBuffers ();
}

void key (unsigned char keyPressed, int x, int y) {
  switch (keyPressed) {
  case 'f':
    if (fullScreen) {
      glutReshapeWindow (camera.getScreenWidth (), camera.getScreenHeight ());
      fullScreen = false;
    } else {
      glutFullScreen ();
      fullScreen = true;
    }
    break;
  case '+':
    remesh_counter += 1;
    if (remesh_counter == 5)
      remesh_counter = 4;
    switch (remesh_counter){
    case 0:
      mesh = o_mesh;
      break;
    case 1:
      mesh = one_mesh;
      break;
    case 2:
      mesh = two_mesh;
      break;
    case 3:
      mesh = three_mesh;
      break;
    case 4:
      mesh = four_mesh;
      break;
    default:
      break;
    }
    break;
  case '-':
    remesh_counter -= 1;
    if (remesh_counter == -1)
      remesh_counter = 0;
    switch (remesh_counter){
    case 0:
      mesh = o_mesh;
      break;
    case 1:
      mesh = one_mesh;
      break;
    case 2:
      mesh = two_mesh;
      break;
    case 3:
      mesh = three_mesh;
      break;
    case 4:
      mesh = four_mesh;
      break;
    default:
      break;
    }
    break;
  case 'q':
  case 27:
    exit (0);
    break;
  case 'w':
    GLint mode[2];
		glGetIntegerv (GL_POLYGON_MODE, mode);
		glPolygonMode (GL_FRONT_AND_BACK, mode[1] ==  GL_FILL ? GL_LINE : GL_FILL);
    break;
    break;
  default:
    printUsage ();
    break;
  }
}

void mouse (int button, int state, int x, int y) {
  camera.handleMouseClickEvent (button, state, x, y);
}

void motion (int x, int y) {
  camera.handleMouseMoveEvent (x, y);
}

void idle () {
  static float lastTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
  static unsigned int counter = 0;
  counter++;
  float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
  if (currentTime - lastTime >= 1000.0f) {
    FPS = counter;
    counter = 0;
    static char winTitle [128];
    unsigned int numOfTriangles = mesh.T.size ();
    sprintf (winTitle, "Number Of Triangles: %d - FPS: %d", numOfTriangles, FPS);
    glutSetWindowTitle (winTitle);
    lastTime = currentTime;
  }
  glutPostRedisplay ();
}

int main (int argc, char ** argv) {
  if (argc > 2) {
    printUsage ();
    exit (1);
  }
  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
  glutInitWindowSize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
  window = glutCreateWindow (appTitle.c_str ());
  init (argc == 2 ? argv[1] : DEFAULT_MESH_FILE.c_str ());
  glutIdleFunc (idle);
  glutReshapeFunc (reshape);
  glutDisplayFunc (display);
  glutKeyboardFunc (key);
  glutMotionFunc (motion);
  glutMouseFunc (mouse);
  printUsage ();
  glutMainLoop ();
  return 0;
}
