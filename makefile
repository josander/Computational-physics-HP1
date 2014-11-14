
CC = gcc
CFLAGS = -O3
LIBS = -lm

HEADERS = initfcc.h alpotential.h
OBJECTS = initfcc.o alpotential.o MD_main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

