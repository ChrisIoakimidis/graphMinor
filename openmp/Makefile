# Project 1 PaDS Makefile
CC=gcc
CFLAGS=-fopenmp  
DEPS= myStructs.h mmio.h myMatrix.h graphMinor.h


all: main test

test: test.o myStructs.o mmio.o myMatrix.o graphMinor.o
	$(CC) $(CFLAGS) -o $@ $^

main: main.o myStructs.o mmio.o myMatrix.o graphMinor.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
