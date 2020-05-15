# Lifting Simplices to Find Injectivity

![](figure/teaser.png)

[Xingyi Du](https://duxingyi-charles.github.io/) [Noam Aigerman](https://research.adobe.com/person/noam-aigerman/), [Qingnan Zhou](https://research.adobe.com/person/qingnan-zhou/), [Shahar Kovalsky](https://shaharkov.github.io/), [Yajie Yan](https://yajieyan.github.io/), [Danny Kaufman](http://dannykaufman.io/), [Tao Ju](https://www.cse.wustl.edu/~taoju/)<br/>
*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2020)*<br/>


# TLC-QN
TLC (Total Lifted Content) energy minimization by quasi-Newton method.

tested on macOS 10.15.1 (Apple Clang 11.0.0) and Ubuntu 18.04.3 LTS (gcc 7.4.0).

## install NLopt

We use the lbfgs quasi-Newton method implemented in NLopt.

### macOS
install NLopt (version 2.6.1) by homebrew

    brew install nlopt

### Ubuntu
    sudo apt-get install libnlopt-dev


## compile
    g++ lifted_test.cpp -lnlopt -lm -O3 -o nloptSolve

`-lnlopt -lm` links to NLopt library and math library.

## how to use

the executable `nloptSolve` asks for 3 options: a path to data file, a path to solver options file, and a path to the file 
to store the result.

    ./nloptSolve [data_file] [solver_options_file] [result_file]

example:

    ./nloptSolve test/lifted test/lifted_solver_options test/lifted_res

## data format

data_file

    num_restVert dimension_restVert
    ... num_restVert x dimension_restVert matrix ...
    num_initVert dimension_initVert
    ... num_initVert x dimension_initVert matrix ...
    num_simplex simplex_size
    ... num_simplex x simplex_size matrix ...
    num_handles
    ... num_handles x 1 matrix ...
    harmonic OR tutte-uniform
    alpha

result_file
    
    name dims
    data
    ...

    
    
    
    
