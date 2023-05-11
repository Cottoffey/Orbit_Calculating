all: main.cpp creatingModelValues.o creatingModelValues.h 
	g++ main.cpp $(regression_src) creatingModelValues.o -l:libsofa_c.a -o main

creatingModelValues.o: creatingModelValues.cpp Ephemeris.h creatingModelValues.h
	g++ creatingModelValues.cpp -l:libsofa_c.a -c 

processingObservingData: processingObservingData.cpp Ephemeris.h
	g++ processingObservingData.cpp -l:libsofa_c.a -o processingObservingData
	
clean:
	rm *.o main
regression_src = Regression/Cholesky.cpp Regression/GaussNewton.cpp Regression/Matrix.cpp
