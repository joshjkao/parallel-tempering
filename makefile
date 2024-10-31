SRC_DIR=src
CXX=g++
CXXFLAGS=-Wall -std=c++11 -O3 -fopenmp
SRC_FILES=$(wildcard $(SRC_DIR)/*.cpp) main.cpp
OBJ_NAME=run
INCLUDE_PATHS=-I include


all:
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) $(SRC_FILES) -o $(OBJ_NAME)

debug:
	$(CXX) -Wall -g -fopenmp $(INCLUDE_PATHS) $(SRC_FILES) -o $(OBJ_NAME)

clean:
	rm -f $(OBJ)
