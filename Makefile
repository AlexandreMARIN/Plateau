CXX = g++
CXXFLAGS = -Wall -std=c++11
EXEC = solvePlateau
MOD = SurfMesh3D Matrix iterative_solvers Plateau

SRC := $(MOD:=.cpp) $(EXEC:=.cpp)
HDR := $(MOD:=.hpp) R23.hpp
OBJ := $(MOD:=.o) $(EXEC:=.o)



all: init $(EXEC)
	@echo "\nThe executable file is \"solvePlateau\"\n"

init:
	-@mkdir obj
	-@mkdir output
	@echo "\nTwo directories have been created (obj/ and output/).\n"
	@echo "\nPut results in /output.\n"

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(addprefix obj/, $(OBJ))

print:
	@echo "\nSources and headers:\n"
	@echo $(SRC)
	@echo $(HDR)
	@echo "Contents of the directories:\n"
	@ls src/
	@ls hdr/


$(MOD:=.o) :%.o : src/%.cpp hdr/%.hpp
	$(CXX) $(CXXFLAGS) -o obj/$@ -c $(filter %.cpp, $<)

$(EXEC:=.o) :%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -o obj/$@ -c $<

clean:
	-@rm $(EXEC)
	-@rm obj/*.o
	-@rm *~

.PHONY: all clean init print
