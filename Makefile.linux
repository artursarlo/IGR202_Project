CIBLE = main
SRCS =  Main.cpp Camera.cpp Mesh.cpp
LIBS =  -lglut -lGLU -lGL -lm 

CC = g++
CPP = g++

FLAGS = -Wall -O2 

CFLAGS = $(FLAGS)
CXXFLAGS = $(FLAGS)

OBJS = $(SRCS:.cpp=.o)   

$(CIBLE): $(OBJS)
	g++ $(LDFLAGS) -o $(CIBLE) $(OBJS) $(LIBS)
clean:
	rm -f  *~  $(CIBLE) $(OBJS)

Camera.o: Camera.cpp Camera.h Vec3.h
Mesh.o: Mesh.cpp Mesh.h Vec3.h
Main.o: Main.cpp Vec3.h Camera.h Mesh.h



