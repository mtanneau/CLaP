# Lift-and-Project experiments

## Overview

## Running the experiments
Shell command to run from this directory

### Log files

The log files are located in the `log/` folder.
The naming convention is
```
<finst>_<MICPSolver>_<CGCPSolver>_<Normalization>_<Rounds>.cpx
<finst>_<MICPSolver>_<CGCPSolver>_<Normalization>_<Rounds>.log
```
for the CPLEX output and logging output, respectively, where
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
cat exp/instances_misocp.txt | grep "flay" | parallel -j NJOBS --joblog exp/landp/jobs.log "julia --project=. --sysimage=exp/landp/JuliaLandP.so landp.jl --MICPSolver CPLEX --CGCPSolver Gurobi --Rounds 200 --Normalization Conic --TimeLimit 1200.0 dat/cblib/{}.cbf > exp/landp/log/{}_CPX_GRB_SCN_200.cpx 2>exp/landp/log/{}_CPX_GRB_SCN_200.log"
```

## Settings

* `MICPSolver`: MICP solver.
    Possible values are `CPLEX` and `Gurobi`
    
* `CGCPSolver`: CGCP solver.
    Possible values are `CPLEX`, `Gurobi` and `Mosek`.

* `Normalization`: normalization condition in CGCP.
    Possible values are
    * `Alpha2`: α-normalization `|α| ⩽ 1`
    * `Interior`: `α'γ ⩽ 1` where `γ = x* - x_` and `x*` is an interior point
    * `Trivial`: `u0 + v0 ⩽ 1`
    * `Standard`: `|λ| + |μ| + u0 + v0 ⩽ 1`
    * `Uniform`: `|λ| + |μ| ⩽ 1`

* `Rounds`: Maximum number of cutting plane rounds.
    Note that K*-cut used to refine the outer approximation are not counted as rounds.

* `TimeLimit`: Time limit, in seconds.

To view more information, just run
```bash
julia landp.jl --help
```