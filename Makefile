CXX = g++
CXXFLAGS = -Wall -std=c++11
SRC = Matrix iterative_solvers Plateau

SRC := $(SRC:=.cpp)
HDR := $(SRC:.cpp=.hpp)
OBJ := $(SRC:.cpp=.o)

all: $(OBJ)

print:
	@echo $(SRC)
	@echo $(HDR)


$(OBJ) :%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -o obj/$@ -c $<

.PHONY: all clean
