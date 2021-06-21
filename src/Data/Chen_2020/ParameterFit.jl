using LinearAlgebra, Plots, GLM, DataFrames # Polynomials
Data = DataFrame(j0_NMC =[2.8e-1;3.48e-1;4.53e-1;5.44e-1;6.02e-1;], j0_Graph = [2.95e-2;3.54e-2;5.76e-2;8.32e-2;1.30e-1],T = [25;30;40;50;60])

cs_max_pos = 51765
Ce0 = 1000
model = lm(@formula(k_NMC ~ (2*j0_NMC/(cs_max_pos*sqrt(ce0)))), Data)
plot!(df.x, predict(model, df), label="model")
#k_Graph = 2*j0_Graph/(cs_max_pos*sqrt(Ce0))

# jN = fit(T, j0_NMC,2)
# jG = fit(T,j0_Graph,2)

# plt = plot(T,j0_NMC,label = "Data")
# plot!(jN,T[1],T[end], label="Fit")
# savefig("j0_NMC.png")

# plt2 = plot(T,j0_Graph,reuse = false, label = "Data")
# plot!(jG,T[1],T[end], label="Fit")
# savefig("j0_Graph.png")

return k_NMC