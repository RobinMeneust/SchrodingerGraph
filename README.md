# Schrodinger's equation Solver

Project written in C and bash<br>
Solves the time independent Schrodinger equation for 3 different potential cases:
- V = 0 <br>
  <img width="400px" src="https://github.com/RobinMeneust/SchrodingerGraph/blob/main/images/potential_1.jpg?raw=true"/><br>

- V is a step<br>
  <img width="400px" src="https://github.com/RobinMeneust/SchrodingerGraph/blob/main/images/potential_2.jpg?raw=true"/><br>
  
- V is rectangular<br>
  <img width="400px" src="https://github.com/RobinMeneust/SchrodingerGraph/blob/main/images/potential_3.jpg?raw=true"/><br>

## Dependencies
  - GSL - GNU Scientific Library
  - Gnuplot

## Installation

Use `make` to compile the code

## Execution

Use `make run` to execute it and generate the graph<br>
You can also directly run the executable in the bin folder to get the data files and/or launch the script gnuplotDraw.bash to get the graph from those files

## Miscellaneous

### To generates the doxygen documentation
````
make doc
````

### To clean the doxygen documentation folder

````
make cleandoc
````

### To remove the .o and .dat files

````
make clean
````

### Notes
For high values (given by the user) the program may take longer or get stuck
