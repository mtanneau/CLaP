# Lift-and-Project experiments

## Overview

## Running the experiments
Shell command to run from this directory

### Single instance

```bash
julia --project=@. script.jl [options] <instance>.cbf
```
for instance:
```bash
julia --project=@. script.jl  --MICPSolver CPLEX --CGCPSolver Gurobi --TimeLimit 120.0 ../../examples/dat/misocp2.cbf
```

### Multiple instances with GNU parallel

## Settings

* `MICPSolver`: MICP solver.
    Possible values are `CPLEX` and `Gurobi`
* `CGCPSolver`: CGCP solver.
    Possible values are `CPLEX`, `Gurobi` and `Mosek`.
* `MaxRounds`
* `TimeLimit`: Time limit, in seconds.
* `Normalization`: normalization condition in CGCP.
    Possible values are
    * `Alpha2`: α-normalization
    * `Conic`: conic normalization
* `