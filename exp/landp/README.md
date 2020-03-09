# Lift-and-Project experiments

## Overview

## Running the experiments
Shell command to run from this directory

### Log files

The log files are located in the `log/` folder.
Each log file should be named as follows:
```
<finst>_<MICPSolver>_<CGCPSolver>_<Normalization>_<Rounds>.cpx
```
for the CPLEX output, and
```
<finst>_<MICPSolver>_<CGCPSolver>_<Normalization>_<Rounds>.log
```
for the log file, where
* `<finst>` is the instance name
* `<MICPSolver>` and `<CGCPSolver>` are the MICP and CGCP solvers' names, respectively
* `<Normalization>` is the normalization
* `<Rounds>` is the (maximum) number of rounds

Note that all options are displayed in the log file, so this convention would not affect the parsing of log files.

### Single instance

```bash
julia --project=@. run_landp.jl [options] <instance>.cbf
```
for instance:
```bash
julia --project=@. run_landp.jl  --MICPSolver CPLEX --CGCPSolver Gurobi --TimeLimit 120.0 --Rounds 10 ../../examples/dat/misocp2.cbf
```

### Multiple instances (with GNU `parallel`)

```bash
cat ../../dat/instances_misocp.txt | grep "flay" | parallel -j1 "julia --project=@. run_landp.jl --CGCPSolver CPLEX --Rounds 10 --Normalization Conic ../../dat/cblib/{}.cbf > log/{}_CPX_CPX_SCN_10.cpx 2> log/{}_CPX_CPX_SCN_10.log
```
will run all `flay` instances using 1 thread.

## Settings

* `MICPSolver`: MICP solver.
    Possible values are `CPLEX` and `Gurobi`
* `CGCPSolver`: CGCP solver.
    Possible values are `CPLEX`, `Gurobi` and `Mosek`.
* `Normalization`: normalization condition in CGCP.
    Possible values are
    * `Alpha2`: α-normalization `|α| ⩽ 1`
    * `Conic`: conic normalization `|ρ| + |μ| + u0 + v0 ⩽ 1`
    * `PureConic`: pure conic normalization `|ρ| + |μ| ⩽ 1`
* `Rounds`: Maximum rounds of cutting planes.
    Once the callback is called that many times, it is de-activated.
* `TimeLimit`: Time limit, in seconds.