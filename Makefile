HEAD = $(wildcard ./include/*.h)
SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))
CC = gcc 
PROG=./bin/schrodingerSolver

all: $(PROG) 

$(PROG) : $(OBJ)
	$(CC) $^ -o $@ -lgsl -lgslcblas -lm

obj/%.o: src/%.c $(HEAD)
	$(CC) -c $< -o $@

.PHONY : cleanlinux cleanwin doc run

cleanlinux:
	rm obj/*.o; rm graphs/*.jpg; rm data/*.dat;

cleanwin:
	del obj\*.o

cleandoclinux:
	rm doxygen/html

cleandocwin:
	del doxygen\html

doc:
	doxygen doxygen/Doxyfile

run:
	$(PROG) && ./gnuplotDraw.bash