HEAD = $(wildcard ./include/*.h)
SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))
CC = gcc 
PROG=./bin/schrodingerSolver

all: $(PROG) 

$(PROG) : $(OBJ)
	$(CC) $^ -o $@ -lgsl -lgslcblas -lm

obj/%.o: src/%.c $(HEAD)
	$(CC) -c $< -Wall -o $@

.PHONY : clean cleandoc doc run

clean:
	rm obj/*.o; rm data/*.dat;

cleandoc:
	rm -r doxygen/html

doc:
	doxygen doxygen/Doxyfile

run:
	$(PROG) && ./gnuplotDraw.bash
