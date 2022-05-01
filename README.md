# Schrodinger's equation Solver

Project written in C and bash
Solves the time independent Schrodinger equation for 3 different potential cases:
  - V = 0
  ``math
  V(x) = \begin{cases}
  \infty &\text{if } x \leq 0 \\
  0 &\text{if } x \in \lbrack 0;L \rbrack  \\
  \infty &\text{if } x \geq L \\
  \end{cases}
  ``
  - V is rectangular
  ``math
  V(x) = \begin{cases}
  \infty &\text{if } x \leq 0 \\
  0 &\text{if } x \in \lbrack 0;a \rbrack  \\
  V0 &\text{if } x \in \lbrack a;b \rbrack  \\
  0 &\text{if } x \in \lbrack b;L \rbrack  \\
  \infty &\text{if } x \geq L \\
  \end{cases}
  ``
  - V is a step
  ``math
  V(x) = \begin{cases}
  \infty &\text{if } x \leq 0 \\
  0 &\text{if } x \in \lbrack 0;a \rbrack  \\
  V0 &\text{if } x \in \lbrack a;L \rbrack  \\
  \infty &\text{if } x \geq L \\
  \end{cases}
  ``

## Dependencies
  - GSL - GNU Scientific Library
  - Gnuplot

## Installation

Use `make` to compile the code

## Execution

Use `make run` to execute it and generate the graph<br>
You can also directly run the executable in the bin folder and/or launch the script gnuplotDraw.bash to get the graph


## Miscellaneous

### To remove the .o files you can do the following

#### Windows
````
make cleanwin
````

#### Linux
````
make cleanlinux
````


