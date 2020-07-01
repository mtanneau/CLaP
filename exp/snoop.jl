include("landp/landp.jl")

# Run the experiment on a small instance
# This one contains Second Order Cones
for
    nrm in [:Alpha2, :Trivial, :Interior, :Standard, :Uniform],
    cgcp in [:Gurobi, :Mosek],
    finst in ["nvs03", "clay0203h"]

    main([
        "--MICPSolver", "CPLEX",
        "--CGCPSolver", "$cgcp",
        "--Rounds", "20",
        "--TimeLimit", "30.0",
        "--Normalization", "$nrm",
        joinpath(@__DIR__, "../dat/cblib/$(finst).cbf")
    ])

end