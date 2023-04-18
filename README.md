# Cool Project Name Here

## Installation and Running

The code can be built using the traditional CMake build procedure, _i.e._:
```
cmake -S . -B build
cmake --build build -- -jN
```
where N is the number of cores to use to compile.

The code can then be run with: `./build/src/bcnutau`. Currently, the matrix element is just set to 1,
and the cuts are hard coded to be p<sub>T</sub> larger than 30 GeV. The integrator is saved into a file
called "results.chili" to enable quicker reloads of the integration for event generation purposes.
