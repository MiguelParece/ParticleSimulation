CC=gcc
CFLAGS=-O2
LDFLAGS=-lm

all: parsim

parsim: parsim.c
	$(CC) $(CFLAGS) -o parsim parsim.c $(LDFLAGS)

clean:
	rm -f parsim