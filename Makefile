all: main.cpp creatingModelValues.o creatingModelValues.h sofa
	g++ main.cpp $(regression_src) creatingModelValues.o $(sofa_lib) -o main

creatingModelValues.o: creatingModelValues.cpp Ephemeris.h creatingModelValues.h sofa
	g++ creatingModelValues.cpp -l $(sofa_lib) -c

processingObservingData: processingObservingData.cpp Ephemeris.h sofa
	g++ processingObservingData.cpp -l $(sofa_lib) -o processingObservingData
clean:
	rm -f *.o main

regression_src = Regression/Cholesky.cpp Regression/GaussNewton.cpp Regression/Matrix.cpp
sofa_lib = sofa/c/src/libsofa_c.a
