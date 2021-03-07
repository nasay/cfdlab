# Include files
SOURCES=initLB.c visualLB.c boundary.c collision.c streaming.c computeCellValues.c flag.c main.c helper.c checks.c

# Compiler
# --------
CC=dpcpp

CFLAGS=-Wpedantic -Wall -O3 -fsycl -std=c++17

# Linker flags
# ------------
LDFLAGS = -lm -lpthread -lrt -ldl

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)
debug: CFLAGS += -DDEBUG -g -O0 
debug: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
