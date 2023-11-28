# graphMinor
This is my implementation of the graph minor alorithm for the 1st assignment of the Parallel and Distributed Systems class.
## For the project to run:
Inside the "myMatrix.h" file the paths to the matrix market dataset and the cluster vector datasets must be set for the parsing.
To build each project, there is a separate makefile for each API inside the corresponding folder. Use "make test" to create the testing program.
## Usage:
./test matrix-name (Without the extension .mtx)
(THE MATRIX AND THE CLUSTER NAMES MUST MATCH, ex. "tx2010.mtx" for the matrix and "tx2010.txt" for the cluster vector)
The number of threads is defined in graphMinor.h as P
## NOTE:
The road_usa graph was too big for my computer to even parse it, even the .mat file couldn't be loaded inside matlab (the application was killed), so only the matrices that I could parse were used.
