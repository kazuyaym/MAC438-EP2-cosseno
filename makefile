CC     = gcc
CFLAGS = -O2 -Wall -pedantic
LIBS   = -lgmp -pthread -lm

all:
	$(CC) $(CFLAGS) ep2.c -o ep2 $(LIBS)
	
clean: 
	rm -rf ep2