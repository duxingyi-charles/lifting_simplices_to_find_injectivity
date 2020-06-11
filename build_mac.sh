#!/usr/bin/env bash

echo "-------- install nlopt ------------"
brew ls --version nlopt || brew install nlopt


mkdir build
cd build

echo "-------- CMake -------------------"
cmake -DCMAKE_BUILD_TYPE=Release ..

echo "-------- build -------------------"
make