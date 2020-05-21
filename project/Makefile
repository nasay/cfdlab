# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c flag.c main.c helper.c checks.c

# Compiler
# --------
CC=dpcpp

CFLAGS=-Wpedantic -Wall -O3 -fsycl

# Linker flags
# ------------
LDFLAGS = -lm -lpthread -lrt -ldl

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)
debug: CFLAGS += -DDEBUG -g -ftree-vectorizer-verbose=5 -fopt-info-vec-optimized
debug: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
