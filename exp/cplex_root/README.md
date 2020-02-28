# Root node experiments with CPLEX

## Overview

The `run_cplex_root.jl` script does the following:
1. Read an instance in CBF format
1. Convert the problem instance to standard form
1. Solve the MICP instance with CPLEX and the following parameters:
    * 1 thread
    * 1hr time-limit
    * Outer-approximation algorithm: `CPXPARAM_MIP_Strategy_MIQCPStrat=2`
    * No presolve
    * No heuristics
    * Root node only `CPXPARAM_MIP_Limits_Nodes=0`
    * Cut generation
        * Aggressive lift-and-project: `CPXPARAM_MIP_Cuts_LiftProj=3`
        * User-specified limit on cuts passes (see below)
        * No limit on the number of added cuts `CPXPARAM_MIP_Limits_CutsFactor=1e30`
        * All other cuts de-activated: `CPXPARAM_MIP_Cuts_XXX=-1`
    
1. Save best solution to `res/<instance>.sol`

## Running one instance

```bash
julia --project=@. run_cplex_root.jl <instance_name>.cbf <N>
```
where `N` is the (integer) number of cut passes, i.e., `CPXPARAM_MIP_Limits_CutPasses` will be set to `N`.

## Running all instances

```bash
cat ../../dat/instances_misocp.txt | parallel -j1 "julia --project=@. run_cplex_root.jl ../../dat/cblib/{}.cbf <N> > log/cpx_rootLaP_<N>_{}.log 2>&1"
```