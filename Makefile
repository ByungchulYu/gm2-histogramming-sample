all: Byu2Histograms.o runHistogramming 

Byu2Histograms.o: Byu2Histograms.cc Makefile
	g++ -c -Wall -Wextra Byu2Histograms.cc $(shell root-config --cflags) -ffast-math -O2

runHistogramming: runHistogramming.o
	g++ -o runHistogramming Byu2Histograms.o runHistogramming.o $(shell root-config --libs) -lMinuit

runHistogramming.o: runHistogramming.cc Makefile
	g++ -c -Wall -Wextra runHistogramming.cc $(shell root-config --cflags) -ffast-math -O2

clean:
	rm *.o runHistogramming 