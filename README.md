# Conic Lift-and-Project (CLaP)

MOI-based framework for separating lift-and-project cuts in mixed-integer conic optimization.
## Building a Julia system image

```bash
julia --project= --trace-compile=precompile.jl exp/snoop.jl
```
then
```julia
using PackageCompiler
PackageCompiler.create_sysimage([:JuMP, :MathOptInterface, :LinearAlgebra, :ArgParse, :TimerOutputs, :Gurobi, :CPLEX, :Logging, :Mosek, :MosekTools, :CLaP]; project=".", sysimage_path="JuliaLandP.so", precompile_statements_file="precompile.jl")
```

To use the system image, simply add the `-JJuliaLandP.so` option when calling Julia, for instance:
```
julia -JJuliaLandP.so --project=@. landp.jl --CGCPSolver Gurobi --TimeLimit 120.0 --Rounds 10 ../../examples/dat/misocp2.cbf
```