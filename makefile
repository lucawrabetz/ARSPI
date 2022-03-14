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
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)M3.cpp -o $(BINPATH)main $(SRCPATH)main.cpp -I$(INCPATHMAC) $(CPPLIBMAC) 

sp:
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)M3.cpp -o $(BINPATH)sp $(SRCPATH)sp.cpp -I$(INCPATHMAC) $(CPPLIBMAC) 

enum:
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)M3.cpp -o $(BINPATH)enum $(SRCPATH)enum.cpp -I$(INCPATHMAC) $(CPPLIBMAC) 

clean:
	rm -rf *.o $(BINPATH)*.dSYM $(BINPATH)main
