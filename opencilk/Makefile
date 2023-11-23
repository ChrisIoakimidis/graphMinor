# Project 1 PaDS Makefile
CC=clang
CFLAGS=-fopencilk  
DEPS= myStructs.h mmio.h myMatrix.h graphMinor.h


all: main test

main: main.o myStructs.o mmio.o myMatrix.o graphMinor.o
	$(CC) $(CFLAGS) -o $@ $^

test: test.o myStructs.o mmio.o myMatrix.o graphMinor.o
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
