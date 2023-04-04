creatingModelValues: creatingModelValues.cpp Ephemeris.h
	g++ creatingModelValues.cpp -l:libsofa_c.a -o creatingModelValues

processingObservingData: processingObservingData.cpp Ephemeris.h
	g++ processingObservingData.cpp -l:libsofa_c.a -o processingObservingData
