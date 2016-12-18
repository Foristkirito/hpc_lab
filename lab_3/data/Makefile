CC=gcc
main: main.o
	$(CC) -o $@ $<
main.o: main.c
	$(CC) -c main.c
read_a: read_a.o
	$(CC) -o $@ $<
read_a.o: read_a.c
	$(CC) -c read_a.c
read_b: read_b.o
	$(CC) -o $@ $<
read_b.o: read_b.c
	$(CC) -c read_b.c
clean:
	rm -f *.o