CXX = g++
CXXFLAGS = -Wall -std=c++11
SRC = test_SurfMesh3D SurfMesh3D Matrix iterative_solvers Plateau

SRC := $(SRC:=.cpp)
HDR := $(SRC:.cpp=.hpp)
OBJ := $(SRC:.cpp=.o)

all: test_SurfMesh3D

test_SurfMesh3D: $(OBJ)
	$(CXX) $(CXXFLAGS) -o obj/$@ $(addprefix obj/, $(OBJ))

print:
	@echo $(SRC)
	@echo $(HDR)


$(OBJ) :%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -o obj/$@ -c $<

.PHONY: all clean
