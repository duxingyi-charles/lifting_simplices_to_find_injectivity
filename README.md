# Lifting Simplices to Find Injectivity

![](figure/teaser.png)

[Xingyi Du](https://duxingyi-charles.github.io/), [Noam Aigerman](https://research.adobe.com/person/noam-aigerman/), [Qingnan Zhou](https://research.adobe.com/person/qingnan-zhou/), [Shahar Kovalsky](https://shaharkov.github.io/), [Yajie Yan](https://yajieyan.github.io/), [Danny Kaufman](https://research.adobe.com/person/danny-kaufman/), [Tao Ju](https://www.cse.wustl.edu/~taoju/)<br/>
*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2020)*<br/>

## Abstract

Mapping a source mesh into a target domain while preserving local injectivity is an important but highly non-trivial task. Existing methods either require an already-injective starting configuration, which is often not available, or rely on sophisticated solving schemes. We propose a novel energy form, called Total Lifted Content (**TLC**), that is equipped with theoretical properties desirable for injectivity optimization. By lifting the simplices of the mesh into a higher dimension and measuring their contents (2D area or 3D volume) there, **TLC** is smooth over the entire embedding space and its global minima are always injective. The energy is simple to minimize using standard gradient-based solvers. Our method achieved _100_% success rate on an extensive benchmark of embedding problems for triangular and tetrahedral meshes, on which existing methods only have varied success.

## TLC-QN

Here we release TLC-QN, a program that finds injective mapping by minimize TLC (Total Lifted Content) energy using quasi-Newton method.

This program has been tested on macOS 10.15.3 (Apple Clang 11.0.3), Ubuntu 18.04.3 LTS (gcc 7.4.0) and Windows 10 (visual studio 2019).

A similar program based on projected Newton method will be released soon.

## Dependency

We use the lbfgs quasi-Newton method implemented in [NLopt](https://nlopt.readthedocs.io/en/latest/).

### macOS
install NLopt (version 2.6.1) by homebrew

    brew install nlopt

### Ubuntu
    sudo apt-get install libnlopt-dev
    
### Window

Source code, DLLs can be downloaded from [NLopt official site](https://nlopt.readthedocs.io/en/latest/NLopt_on_Windows/).


## Compile

Use the following commonds on Mac/Linux,
    
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

You may need to edit NLopt paths in CMakeLists.txt to proper directories on your machine.

## How to use

The executable `findInjective` asks for 3 arguments: a path to input data file, a path to solver options file, and a path to the file to store the result.

    ./findInjective [input_file] [solver_options_file] [result_file]

example:

    ./findInjective example/input example/solver_options example/result

## data format

input_file

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
 
solver_options_file


result_file
    
    name dims
    data
    ...

    
    
    
    
