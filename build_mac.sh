#!/usr/bin/env bash

echo "-------- install eigen ------------"
brew ls --version eigen || brew install eigen
echo "-------- install suitesparse ------"
brew ls --version eigen || brew install suitesparse
echo "-------- install libomp -----------"
brew ls --version eigen || brew install libomp

mkdir build
cd build

echo "-------- CMake -------------------"
cmake -DCMAKE_BUILD_TYPE=Release ..

echo "-------- build -------------------"
make