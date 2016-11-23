CC=gcc
VPATH=../../engine
test_bign: main.o
	$(CC) -o $@ $<
main.o: main.c proj_engine.h
	$(CC) -c main.c
clean:
	rm -f *.o