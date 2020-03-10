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
cat ../../dat/instances_misocp.txt | grep "flay" | parallel -j1 "julia --project=@. run_landp.jl --MICPSolver CPLEX --CGCPSolver Gurobi --Rounds 200 --Normalization Conic --TimeLimit 7200.0 ../../dat/cblib/{}.cbf > log/{}_CPX_GRB_SCN_200.cpx 2>log/{}_CPX_GRB_SCN_200.log"
```

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

To view more information, just run
```bash
julia run_landp.jl --help
```

## Building a Julia sysimage

```bash
$ julia --project=@. --trace-compile=precomp.jl snoop.jl
```
then
```julia
using PackageCompiler
PackageCompiler.create_sysimage([:JuMP, :MathOptInterface, :LinearAlgebra, :ArgParse, :TimerOutputs, :Gurobi, :CPLEX, :Logging]; project="../..", sysimage_path="JuliaLandP.so", precompile_statements_file="precomp.jl")
```

To use the system image, simply add the `-JJuliaLandP.so` option when calling Julia, for instance:
```
julia -JJuliaLandP.so --project=@. run_landp.jl --CGCPSolver Gurobi --TimeLimit 120.0 --Rounds 10 ../../examples/dat/misocp2.cbf
```