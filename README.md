# Conic Lift-and-Project (CLaP)

This repository contains code 
[]()

MOI-based framework for separating lift-and-project cuts in mixed-integer conic optimization.

## Overview

This repository is structured as follows:
* `dat/` contains test instances. Toy instances are in `dat/examples/`
* `exp/` contains scripts to run the different experiments.
* `src/` contains source code for generating lift-and-project cuts

## Installation

### Dependencies

1. Download and install Julia, e.g., from [here](https://julialang.org/downloads/).
2. Download and install CPLEX 12.10 and, optionally, Gurobi 9.0 and Mosek 9.2.
Other versions may work, but have not been tested with the present code.
3. Download this repository and install the required Julia packages as follows:
```
julia> ]
pkg> activate .
pkg> instantiate
```
Additional steps may be needed to build the solvers' wrappers, see specific instructions in the corresponding packages.

### Data

Test instances are MISOCP instances from the [Conic benchmark library](http://cblib.zib.de/).
The list of instances that were used is in the file `exp/instances_misocp.txt`.

1. Download test instances. To do so, you can either
    * Run the `dat/download_cblib.sh` (on a Linux machine), which will download the entire CBLIB collection
    * Download only the test instances

2. Convert all instances to standard form 
    ```
    julia --project exp/convert_instances.jl dat/<instance>.cbf dat/<instance_sf>.cbf
    ```
    This last step ensures that all rotated second-order cone constraints are reformulated as second-order cone constraints.


### Building a Julia system image

A system image is not required for the code to work, but will speed up computations.
This is done as follows.

1. Generate precompile statements
    ```bash
    julia --project= --trace-compile=precompile.jl exp/snoop.jl
    ```
    This script excepts some test instances to be located in a `dat/cblib/` directory.

2. Generate the system image
    ```julia
    using PackageCompiler
    PackageCompiler.create_sysimage([:JuMP, :MathOptInterface, :LinearAlgebra, :ArgParse, :TimerOutputs, :Gurobi, :CPLEX, :Logging, :Mosek, :MosekTools, :CLaP]; project=".", sysimage_path="JuliaLandP.so", precompile_statements_file="precompile.jl")
    ```

3. To use the system image, simply add the `--sysimage=JuliaLandP.so` option when calling Julia, for instance:
    ```
    julia --sysimage=JuliaLandP.so --project exp/landp/landp.jl --CGCPSolver Gurobi --TimeLimit 120.0 --Rounds 10 dat/examples/misocp2.cbf
    ```

## Running experiments

See instructions in the following directories:
* `exp/cplex_micp`: Initial MICP solve
* `exp/cplex_root`: CPLEX root node
* `exp/landp`: CGCP-based lift-and-project cut generation