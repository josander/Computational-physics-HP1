
CC = gcc
CFLAGS = -O3
LIBS = -lm

HEADERS = initfcc.h alpotential.h hpfunc.h
OBJECTS = initfcc.o alpotential.o hpfunc.o MD_main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

