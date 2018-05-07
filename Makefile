CC      = g++
CFLAGS  = -std=c++11 -Wall 
LDFLAGS = -fopenmp

all: cross_decomposition 


cross_decomposition: cross_decomposition.o 
	$(CC) -o $@ $^ $(LDFLAGS)

cross_decomposition.o: cross_decomposition.cpp
	$(CC) -c $(CFLAGS) $<  $(LDFLAGS)

.PHONY: clean cleanest

clean:
	rm *.o

cleanest: clean
	rm cross_decomposition	
