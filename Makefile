CC=g++
LDFLAGS=-std=c++11 -O3 -lm
SOURCES=src/main.cpp src/da.cpp src/parser.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=da
INCLUDES=src/node.h src/da.h

all: $(SOURCES) bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o:  %.c  ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o bin/$(EXECUTABLE)
