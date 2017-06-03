CC=gcc
BINARY = test
CFLAGS = -std=c99 -pedantic -Wall -Wextra
DFLAGS = -g -O0 -D_DEBUG -DPRINT_PROGRESS
RFLAGS = -O2

SOURCES = *.c
HEADERS = *.h
OBJECTS = curve.o binfields.o integers.o test.o

debug: $(OBJECTS)
	$(CC) $(CFLAGS) $(DFLAGS) -o $(BINARY) $(OBJECTS)
release: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(BINARY) $(SOURCES) $(RFLAGS) 
clean:
	-rm $(OBJECTS)
	-rm *.stackdump

.PHONY: debug release clean

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(DFLAGS) -c -o $@ $< 
