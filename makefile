# Makefile for ARSPI

# Paths/Environmental Variables
# Project Folder Directories
SRCPATH	  = ./src/
BINPATH	  = ./bin/
DATPATH	  = ./dat/
INCPATH   = ./inc/

# LUCA MACBOOKPRO PATHS (MBP)
GRBMBPPATH = /Library/gurobi902/mac64
GRBINCMBP = $(GRBMBPPATH)/include
CPPLIBMBP = -L$(GRBMBPPATH)/lib/ -lgurobi_c++ -lgurobi90 $(CPPSTDLIB) -lpthread -lm

# BENEDUM SMOD04 PATHS (SMD)
GRBSMDPATH = /home/luw28/gurobi950/linux64
GRBINCSMD = $(GRBSMDPATH)/include
CPPLIBSMD = -L$(GRBSMDPATH)/lib/ -lgurobi_g++5.2 -lgurobi95 $(CPPSTDLIB) -lpthread -lm

# *****************************************************
# Variables to control Makefile operation

CXXMAC = clang++
CXXFLAGSMAC = -Wall -g -std=c++11

CXXSMD = g++
CXXFLAGSSMD = -m64 -std=c++14 -g -Wall -Wextra -pedantic

# ****************************************************

# FOR MBP 
main:
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)M3.cpp -o $(BINPATH)main $(SRCPATH)main.cpp -I$(GRBINCMBP) $(CPPLIBMBP) 

sp:
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)M3.cpp -o $(BINPATH)sp $(SRCPATH)sp.cpp -I$(GRBINCMBP) $(CPPLIBMBP) 

enum:
	$(CXXMAC) $(CXXFLAGSMAC) $(SRCPATH)M3.cpp -o $(BINPATH)enum $(SRCPATH)enum.cpp -I$(GRBINCMBP) $(CPPLIBMBP) 

# FOR SMD 
# main:
# 	$(CXXSMD) $(CXXFLAGSSMD) $(SRCPATH)M3.cpp -o $(BINPATH)main $(SRCPATH)main.cpp -I$(GRBINCSMD) $(CPPLIBSMD) 
# 
# sp:
# 	$(CXXSMD) $(CXXFLAGSSMD) $(SRCPATH)M3.cpp -o $(BINPATH)sp $(SRCPATH)sp.cpp -I$(GRBINCSMD) $(CPPLIBSMD) 
# 
# enum:
# 	$(CXXSMD) $(CXXFLAGSSMD) $(SRCPATH)M3.cpp -o $(BINPATH)enum $(SRCPATH)enum.cpp -I$(GRBINCSMD) $(CPPLIBSMD) 
# 
# clean:
# 	rm -rf *.lp *.o $(BINPATH)*.dSYM $(BINPATH)main
