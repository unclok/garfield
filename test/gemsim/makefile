OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	 `root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) -v
#-fopenmp 

#CXX = /usr/bin/g++48
# Debug flags
CFLAGS += -g

LDFLAGS = -lgfortran `root-config --glibs` -lGeom -lm
LDFLAGS += -L$(LIBDIR) -lGarfield # -fopenmp
LDFLAGS += -g
gasfile: gasfile.C 
	$(CXX) $(CFLAGS) gasfile.C
	$(CXX) -o gasfile gasfile.o $(LDFLAGS) 
	rm gasfile.o
