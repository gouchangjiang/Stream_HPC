path_only = `dirname $(realpath $(lastword $(MAKEFILE_LIST)))`

all:heuristics

LIB_PATH = ${path_only}/lib
INC_PATH = ${path_only}/include
SRC_PATH = ${path_only}/src
OBJ_PATH = ${path_only}/

CPP=g++
PEDANTIC_PARANOID_FREAK =       -O0 -Wshadow -Wcast-align \
				-Waggregate-return -Wmissing-prototypes -Wmissing-declarations \
				-Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations \
				-Wmissing-noreturn -Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions

REASONABLY_CAREFUL_DUDE =	-Wall
NO_PRAYER_FOR_THE_WICKED =	-w -O2 
WARNINGS = $(REASONABLY_CAREFUL_DUDE)
CFLAGS = $(WARNINGS) -m64 -g -DNOASSERT -std=c++14
#CFLAGS = $(PEDANTIC_PARANOID_FREAK) -std=c++11 -g
INCLUDES = -I${INC_PATH}
LIBS = ${LIB_PATH}/heuristics.a 

# all source files have associated object files
node.o: src/node.cpp include/node.hpp
	$(CPP) $(INCLUDES) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

graph.o: src/graph.cpp include/graph.hpp
	$(CPP) $(INCLUDES) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

utils.o: src/utils.cpp include/utils.hpp
	$(CPP) $(INCLUDES) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

heuristics.o: src/algorithms.cpp include/algorithms.hpp include/common.hpp graph.o utils.o node.o
	$(CPP) $(INCLUDES) $(CFLAGS) -o ${OBJ_PATH}/$@ -c $< -fPIC

heuristics:heuristics.o 
	rm -f ${LIB_PATH}/$@.a ${LIB_PATH}/$@.so
	ar rcs ${LIB_PATH}/$@.a $(OBJ_PATH)/heuristics.o $(OBJ_PATH)/graph.o $(OBJ_PATH)/utils.o $(OBJ_PATH)/node.o

clean:
	rm -f ${OBJ_PATH}/*.o *~ 

.PHONY: clean all

.SUFFIXES:
.SECONDARY:
