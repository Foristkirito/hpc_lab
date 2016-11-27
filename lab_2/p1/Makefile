default:stencil

CC = mpiicc
CFLAGS = -O3 -xAVX -openmp
CLIBS  = -openmp

stencil:main.o
	$(CC) $(CLIBS) main.o wxl.a -o stencil

main.o:main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

clean:
	rm -rf stencil *.o out_*

