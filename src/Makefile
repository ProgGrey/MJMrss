CXX := g++
#CXX := clang++
# -fsanitize=address
#CXXFLAGS := -Og -march=native -std=c++14 -g -fsanitize=address -flto -Wall -Wextra -pedantic
CXXFLAGS := -O3 -ffp-contract=fast -march=native -flto -std=c++14 -Wall -Wextra -pedantic

main:main.o state.o dist.o matrices.o
	${CXX} $(CXXFLAGS) main.o state.o dist.o matrices.o -o main

main.o:main.cpp main.hpp state.hpp dist.hpp
	${CXX} $(CXXFLAGS) main.cpp -c -o main.o

state.o:state.cpp state.hpp
	${CXX} $(CXXFLAGS) state.cpp -c -o state.o

dist.o:dist.cpp dist.hpp
	${CXX} $(CXXFLAGS) dist.cpp -c -o dist.o

matrices.o:matrices.cpp matrices.hpp
	${CXX} $(CXXFLAGS) matrices.cpp -c -o matrices.o

# --bug-hunting
check:
	cppcheck --enable=all --language=c++ --max-ctu-depth=100 --std=c++14 .

clean:
	rm *.o main