#!/bin/bash

#############################
#   FUNCTIONS
#############################

# Function used to get the approximated value of the decimal logarithm of the given parameter in scientific notation (e.g 5.2e+8)
decimalLogarithmApprox(){
    echo "${1::1}e+$((${#1} - 1))"
}

displayHelp(){
    echo -e "NAME\n\tscript_syracuse.bash\n\nSYNOPSIS\n\tscript_syracuse.bash -h"
    echo -e "\tscript_syracuse.bash UMIN MAX\n\nDESCRIPTION"
    echo -e "\tGenerates graphs in the jpeg format corresponding to maximum altitude, flight duration and altitude duration."
    echo -e "\tIt also generates a resume named synthese-min-max.txt that provide minimum, maximum and average values for each one"
    echo -e "\tUMIN and UMAX are strictly positive integers and UMIN is lesser than UMAX"
    echo -e "\tUMAX must be lesser than 1,000,000,000,000,000\n\n\t-h\n\t\tdisplay this help and exit.\n"
}

# Check if a directory is available. If it's a file then we have to delete it before and if it doesn't exist we create it
checkDirAvailability(){
    if ! [ -d "$1" ] 
    then
        if [ -e "$1" ]
        then
            echo "ERROR: you have to delete the file '$1' before running this script"
            exit 1
        else
            mkdir "$1"
        fi
    fi
}

# Check if a the given executable exist. If it doesn't then we ask if the user want to create it
checkExecutableAvailability(){
    if ! [ -e "$1" ]
    then
        read -p "The executable syracuse doesn't exist and is required. Do you want to create it ? (y/n)" answer
        if [ "$answer" = "y" ]
        then
            gcc main.c -o syracuse
            if [ $? -eq 0 ]
            then
                echo "syracuse was successfully created"
            else
                echo "ERROR: syracuse could not be created"
                exit 1
            fi
        else
            echo "ERROR: 'syracuse' is required. Please read README.txt and follow the instructions in the installation part"
            exit 1
        fi
    elif ! [ -x "syracuse" ]
    then
        echo "ERROR: 'syracuse' is not an executable. Please read README.txt and follow the instructions in the installation part"
        exit 1
    fi
}

# Get a value in the given file (named after u0), return it and save it in a file whose name is given in the third parameter
# $1: u0     $2: line starting from the last three lines    $3: fileOutput
getDataFromDatFiles(){
    local current_value="$(tail -n 3 output_dir/f${1}.dat | sed -n "${2}p" | cut -d'=' -f2)"
    echo "$1 $current_value" >> "$3"
    echo $current_value
}

max(){
    if [ $1 -gt $2 ]
    then
        echo $1
    else
        echo $2
    fi
}

min(){
    if [ $1 -lt $2 ]
    then
        echo $1
    else
        echo $2
    fi
}

#Adapt the yrange if there is only one y value in the file and thus in the created graph
#$1 : fileInput $2 : max value
adaptYRangeIfSingleValue(){
    if [ $(uniq $1 | wc -l) -eq 1 ]
    then
        # yrange is written in scientific notation when it's a large number
        if [ ${#2} -gt 8 ]
        then
            echo "; set yrange [0:$(decimalLogarithmApprox $(($2 * 4)))]"
        else
            echo "; set yrange [0:$(($2 + 1))]"
        fi
    else
        echo ""
    fi
}

#Adapt the xrange if there is only one value in the created graph
#$1 : U0MIN $2 : U0MAX
adaptXRangeIfSingleValue(){
    if [ ${#1} -gt 8 ]
    then
        echo "; set xrange [$(decimalLogarithmApprox $(($1 / 4))):$(decimalLogarithmApprox $(($2 * 4)))]"
    else
        echo "; set xrange [$1:$(($2 + 1))]"
    fi
}


#############################
#   MAIN PART OF THIS SCRIPT
#############################

checkDirAvailability "graphs"

# We initialize the needed variables
gnuplot_instructions=""
gnuplot_instructionsALL=""

# If there is only one point we don't draw a line, we only give one dot, we also have to specify the range (if we don't do this then it will be adjusted automatically and will display warnings)
# To do so we use the following variables that may be modified. They will be used when giving instructions to gnuplot.

line_style="with lines"

# We get the data needed to create the graphs
gnuplot_instructions="'data/phi.dat' u 1:2 title 'phi(x)' $line_style lc rgb 'blue', 'data/potential.dat' u 1:2 title 'V(x)' $line_style lw 4 lc rgb 'black'"
#gnuplot_instructionsALL="'data/phiALL.dat' u 1:2 title 'phi(x)' $line_style lc rgb 'blue', 'data/potential.dat' u 1:2 title 'V(x)' $line_style lc rgb 'black'"

# We create the graphs
gnuplot -e "reset; set terminal jpeg size 1600, 900; set ylabel 'phi'; set xlabel 'x'; set output 'graphs/schrodingerGraph.jpg'; plot $gnuplot_instructions"
#gnuplot -e "reset; set terminal jpeg size 1600, 900; set ylabel 'phi'; set xlabel 'x'; set output 'graphs/schrodingerGraphALL.jpg'; plot $gnuplot_instructionsALL"
