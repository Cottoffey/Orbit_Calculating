all: main.cpp creatingModelValues.o creatingModelValues.h 
	g++ main.cpp creatingModelValues.o -l:libsofa_c.a -o main

creatingModelValues.o: creatingModelValues.cpp Ephemeris.h creatingModelValues.h
	g++ creatingModelValues.cpp -l:libsofa_c.a -c 

processingObservingData: processingObservingData.cpp Ephemeris.h
	g++ processingObservingData.cpp -l:libsofa_c.a -o processingObservingData
