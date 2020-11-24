# Makefile for ARSPI

# Paths/Environmental Variables
SRCPATH	  = ./src/
BINPATH	  = ./bin/
DATPATH	  = ./dat/
INCPATH   = ./inc/

# MAC PATHS
GRBMACPATH = /Library/gurobi902/mac64
INCPATHMAC = $(GRBMACPATH)/include
CPPLIBMAC = -L$(GRBMACPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm

# *****************************************************
# Variables to control Makefile operation

CXXCCR = g++
CXXFLAGSCCR = -m64 -std=c++14 -g -Wall -Wextra -pedantic

CXXMAC = clang++
CXXFLAGSMAC = -Wall -g -std=c++11

# ****************************************************

# FOR MAC
all:
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)rspi.cpp -o $(BINPATH)main $(SRCPATH)main.cpp -I$(INCPATHMAC) $(CPPLIBMAC)

clean:
	rm -rf *.o $(BINPATH)*.dSYM $(BINPATH)main