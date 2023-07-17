SRC_DIR = src
CXX=g++-13
CXXFLAGS= -Wall -std=c++11 -O2 -fopenmp
SRC_FILES = $(wildcard $(SRC_DIR)/*.cc) main.cc
OBJ_NAME = run
INCLUDE_PATHS = -I include 
LINKER_FLAGS = -fopenmp


all:
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) $(SRC_FILES) -o $(OBJ_NAME)

clean:
	rm -f $(OBJ)
