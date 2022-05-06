///////////////////////////////////
/// Schrodinger's equation Solver
///////////////////////////////////

Project written in C and bash

Solves the time independent Schrodinger equation for 3 different potential cases:

	- V = 0
	- V is a step
	- V is rectangular

///////////////////////////////////
Dependencies:

    GSL - GNU Scientific Library
    Gnuplot

///////////////////////////////////
Installation:

Use make to compile the code

///////////////////////////////////
Execution:

Use make run to execute it and generate the graph
You can also directly run the executable in the bin folder to get the data files and/or launch the script gnuplotDraw.bash to get the graph from those files

///////////////////////////////////
Miscellaneous:

To generates the doxygen documentation

	make doc

To clean the doxygen documentation folder

	make cleandoc

To remove the .o and .dat files

	make clean

Notes

	For high values (given by the user) the program may take longer or get stuck