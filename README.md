# lifted-nlopt
lifted energy minimization by NLopt library

# install NLopt
tested on macOS 10.15.1 (Apple Clang 11.0.0) and Ubuntu 18.04.3 LTS (gcc 7.4.0).

## macOS
install NLopt (version 2.6.1) by homebrew

    brew install nlopt

## Ubuntu
    sudo apt-get install libnlopt-dev


# compile
    g++ lifted_test.cpp -lnlopt -lm -O3 -o nloptSolve

`-lnlopt -lm` links to NLopt library and math library.

# how to use

the executable `nloptSolve` asks for 3 options: a path to data file, a path to solver options file, and a path to a file 
to store the result.

    ./nloptSolve [dataFile] [solver_options_file] [result_file]

example:

    ./nloptSolve test/lifted test/lifted_solver_options test/lifted_res
