CFLAGS = -O
CC = g++

poisson2d: main.o poisson2d.o
	$(CC) $(CFLAGS) -o poisson main.o poisson2d.o

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

poisson2d.o: poisson2d.cpp
	$(CC) $(CFLAGS) -c poisson2d.cpp 

clean: 
	rm *.o
