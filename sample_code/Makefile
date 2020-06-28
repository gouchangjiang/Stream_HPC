path_only = `dirname $(realpath $(lastword $(MAKEFILE_LIST)))`

all: call-heuristics

LIB_PATH = ${path_only}/../lib
INC_PATH = ${path_only}/../include
BIN_PATH = ${path_only}/

CPP = g++
PEDANTIC_PARANOID_FREAK =       -O0 -Wshadow -Wcast-align \
				-Waggregate-return -Wmissing-prototypes -Wmissing-declarations \
				-Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations \
				-Wmissing-noreturn -Wredundant-decls -Wnested-externs \
				-Wpointer-arith -Wwrite-strings -finline-functions
REASONABLY_CAREFUL_DUDE =	-Wall
NO_PRAYER_FOR_THE_WICKED =	-w -O2 
WARNINGS = $(REASONABLY_CAREFUL_DUDE)
CFLAGS = $(WARNINGS) -m64 -g -DNOASSERT -std=c++11
INCLUDES = -I${INC_PATH}
DEFS = 
LDADD =
LIBS = ${LIB_PATH}/heuristics.a 


call-heuristics: main.cpp 
	$(CPP) $(INCLUDES) $(DEFS) $(CFLAGS) $< $(LIBS) $(LDADD) -o ${BIN_PATH}/$@

clean:
	rm -f ${BIN_PATH}/call-heuristics  *.o *~ 

.PHONY: clean all

.SUFFIXES:
.SECONDARY:
