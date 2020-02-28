# Solve MICP instances with CPLEX

## Overview

The `run_cplex_micp.jl` does the following:
1. Read an instance in CBF format
1. Convert the problem instance to standard form
1. Solve the MICP instance with CPLEX and the following parameters:
    * 1 thread
    * 1hr time-limit
    * all other parameters to default
1. Save best solution to `res/<instance>.sol`

## Solving one instance

```bash
julia --project=@. run_cplex_micp.jl <instance_name>.cbf
```

## Solving all instances

```bash
cat ../../dat/instances_misocp.txt | parallel -j2 "julia --project=@. run_cplex_micp.jl ../../dat/cblib/{}.cbf > log/{}.cpx 2>&1"
```