default:main

CC = mpiicc
CFLAGS = -O3 -xAVX -openmp
CLIBS  = -openmp

main:main.o
	$(CC) $(CLIBS) main.o -o main

main.o:main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

clean:
	rm -rf main *.o