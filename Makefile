HEAD = $(wildcard ./include/*.h)
SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))
CC = gcc 
PROG=./bin/test_gsl

all: $(PROG) 

$(PROG) : $(OBJ)
	$(CC) $^ -o $@ -lgsl -lgslcblas -lm

obj/%.o: src/%.c $(HEAD)
	$(CC) -c $< -o $@

.PHONY : cleanlinux cleanwin doc run

cleanlinux:
	rm obj/*.o

cleanwin:
	del obj\*.o

cleandoclinux:
	del doxygen\html

cleandocwin:
	del doxygen\html

doc:
	doxygen doxygen/Doxyfile

run:
	$(PROG)