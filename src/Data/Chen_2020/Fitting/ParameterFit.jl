using LinearAlgebra, Plots, GLM, DataFrames, Polynomials
Data = DataFrame(j0_NMC =[2.8e-1;3.48e-1;4.53e-1;5.44e-1;6.02e-1;], j0_Graph = [2.95e-2;3.54e-2;5.76e-2;8.32e-2;1.30e-1], R_NMC = [1.8;1.48;1.17;1.01;0.94], R_Graph = [14.72;12.51;7.93;5.66;3.74], T = [25;30;40;50;60], cs_max_pos = [63104;63104;63104;63104;63104], cs_max_neg = [33133;33133;33133;33133;33133], ce0 = [1000;1000;1000;1000;1000])

#k_NMC = tuple()
model = lm(@formula(κ ~ (2*j0_NMC/(cs_max_pos*sqrt(ce0)))), Data)
#model = lm(@formula(k_Graph ~ (2*j0_Graph/(cs_max_neg*sqrt(ce0)))), Data)
plot(Data.j0_NMC, predict(model, Data), label="model")
plot!()
#k_Graph = @. 2*j0_Graph/(cs_max_pos*sqrt(ce0))

# jN = fit(T, j0_NMC,2)
# jG = fit(T,j0_Graph,2)

# plt = plot(T,j0_NMC,label = "Data")
# plot!(jN,T[1],T[end], label="Fit")
# savefig("j0_NMC.png")

# plt2 = plot(T,j0_Graph,reuse = false, label = "Data")
# plot!(jG,T[1],T[end], label="Fit")
# savefig("j0_Graph.png")

return