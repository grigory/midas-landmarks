# (c) Microsoft Corporation. All rights reserved.

# 
# source files 
#

SOURCES = midas.cpp \
          mt_random.cpp \
	  midas_timer.cpp

#
# parameters for various compilers
#

GCC_NAME    = g++
GCC_FLAGS   = -Wall -O4 
GCC_LIBS    = -lm -L/usr/lib32/
GCC_DEFINES = 
GCC_OBJECTS = $(SOURCES:.cpp=.o)
GCC_REMOVE  = 'rm' *.o

# CHANGE THESE LINES TO USE YOUR FAVORITE COMPILER
#
CCC      = $(GCC_NAME)
FLAGS    = $(GCC_FLAGS)
LIBS     = $(GCC_LIBS)
DEFINES  = $(GCC_DEFINES)
OBJECTS  = $(GCC_OBJECTS)
REMOVE   = $(GCC_REMOVE)
INCLUDES = -I.

.SUFFIXES: .cpp

midas.exe: $(OBJECTS)
	$(CCC) $(FLAGS) $(DEFINES) $(INCLUDES) $(OBJECTS) $(LIBS) -o midas.exe

all: clean midas.exe

clean: 
	$(REMOVE)	

.cpp.o:
	$(CCC) $(FLAGS) $(DEFINES) $(INCLUDES) -c $<
             
# OBJECT FILES DEPENDENCIES
#
midas.o: midas.cpp graph.h vpair.h typmidas.h mt_random.h evaluator.h \
 dijkstra.h binheap.h myfinder.h landfind.h midas_timer.h

midas_timer.o: midas_timer.cpp midas_timer.h

mt_random.o: mt_random.cpp mt_random.h

