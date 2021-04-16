Compile=g++ -g -std=c++17
dtmc_lc_3d: main.o dtmc_lc_3d.o metropolis.o measurement.o experiment.o
	${Compile} -o dtmc_lc_3d main.o dtmc_lc_3d.o metropolis.o  measurement.o experiment.o
main.o: main.cpp dtmc_lc_3d.h
	${Compile} -o main.o -c main.cpp
dtmc_lc_3d.o: dtmc_lc_3d.cpp dtmc_lc_3d.h
	${Compile} -o dtmc_lc_3d.o -c dtmc_lc_3d.cpp
metropolis.o: metropolis.cpp dtmc_lc_3d.h
	${Compile} -o metropolis.o -c metropolis.cpp
measurement.o: measurement.cpp dtmc_lc_3d.h
	${Compile} -o measurement.o -c measurement.cpp
experiment.o: experiment.cpp dtmc_lc_3d.h
	${Compile} -o experiment.o -c experiment.cpp
clean:
	rm *.o
