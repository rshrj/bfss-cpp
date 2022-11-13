CXX = g++
CXXFLAGS = -Wall

matrix : src/matrix.cpp src/matrix_type.hpp
	$(CXX) $(CXXFLAGS) src/matrix.cpp -o bin/matrix
dump :
	bin/matrix | tee data/output.dat

run :
	bin/matrix

clean :
	rm bin/*
