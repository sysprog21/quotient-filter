CFLAGS = -Wall -O2 -std=gnu99 -g

all: bench

quotient-filter.o: quotient-filter.c quotient-filter.h
	$(CC) $(CFLAGS) -c -o $@ $<

bench: quotient-filter.o bench.c
	$(CC) $(CFLAGS) -o $@ $^

clean:
	$(RM) quotient-filter.o bench
