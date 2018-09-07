compile:
	g++-7 -g -std=c++17 -Wall -Wno-unused-variable -Wno-unused-function -Wno-reorder src/poly_grav.cpp src/poly_run.cpp -o main

compile2:
	gcc -g -Wall -Wno-unused-variable -Wno-unused-function -Wno-reorder src/volInt.c -o bin/volInt


run:
	./main

