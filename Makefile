all: main.cpp creatingModelValues.o creatingModelValues.h libsofa.a
	g++ -std=c++17 main.cpp $(regression_src) creatingModelValues.o $(sofa_lib) -o main

creatingModelValues.o: creatingModelValues.cpp Ephemeris.h creatingModelValues.h libsofa.a
	g++ -std=c++17 creatingModelValues.cpp -l  $(sofa_lib) -c

<<<<<<< HEAD
processingObservingData: processingObservingData.cpp Ephemeris.h libsofa.a
=======
processingObservingData: processingObservingData.cpp Ephemeris.h
>>>>>>> 7e14a1c (Added high percision ephemerises.)
	g++ -std=c++17 processingObservingData.cpp  $(sofa_lib) -o processingObservingData

libsofa.a:
	$(MAKE) -C ./sofa/c/src

clean:
	rm -f *.o main creatingModelValues processingObservingData
	rm -f ./sofa/c/src/libsofa.a
	$(MAKE) clean -C ./sofa/c/src

regression_src = Regression/Cholesky.cpp Regression/GaussNewton.cpp Regression/Matrix.cpp
sofa_lib = sofa/c/src/libsofa_c.a
