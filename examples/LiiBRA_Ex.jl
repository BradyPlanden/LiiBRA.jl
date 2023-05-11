using LiiBRA, Plots
# plotly()
# default(show = true)

#---------- Cell Definition -----------------#
Cell = Construct("LG M50")
Cell.RA.Fs = 4.0 # Modify transfer function sampling frequency
Cell.RA.SamplingT = 0.25 # Modify final system sampling period
Cell.Neg.Ds = 2.0e-14 # Modify negative electrode diffusion constant
Cell.Const.T = 298.15 # Modify cell temperature
Ŝ = collect(1.0:-1:0) # List of SOC points for model generation
SOC = 0.75 # Starting SOC

#---------- Generate & Simulate Model -----------------#
A, B, C, D = Realise(Cell, Ŝ)
Results, tₑ = HPPC(Cell, Ŝ, SOC, 4.0, -3.0, A, B, C, D);

#----------- Plotting ---------------------------#
plot(Results.t, Results.Cell_V;
     legend = :topright,
     color = :blue,
     bottom_margin = 5Plots.mm,
     left_margin = 5Plots.mm,
     right_margin = 15Plots.mm,
     ylabel = "Terminal Voltage (V)",
     xlabel = "Time (s)",
     title = "HPPC Voltage",
     label = "Voltage",
     size = (1280, 720))

plot(Results.t, Results.Ce;
     legend = :topright,
     bottom_margin = 5Plots.mm,
     left_margin = 5Plots.mm,
     right_margin = 15Plots.mm,
     ylabel = "Electrolyte Concen. (mol/m³)",
     xlabel = "Time (s)",
     title = "Electrolyte Concentration",
     label = ["Neg. Separator Interface" "Neg. Current Collector" "Pos. Current Collector" "Pos. Separator Interface"],
     size = (1280, 720))

plot(Results.t, Results.Cseₚ;
     legend = :topright,
     bottom_margin = 5Plots.mm,
     left_margin = 5Plots.mm,
     right_margin = 15Plots.mm,
     ylabel = "Concentration (mol/m³)",
     xlabel = "Time (s)",
     title = "Positive Electrode Concentration",
     label = ["Current Collector" "Separator Interface"],
     size = (1280, 720))

plot(Results.t, Results.Cseₙ;
     legend = :topright,
     bottom_margin = 5Plots.mm,
     left_margin = 5Plots.mm,
     right_margin = 15Plots.mm,
     ylabel = "Concentration (mol/m³)",
     xlabel = "Time [s]",
     title = "Negative Electrode Concentration",
     label = ["Current Collector" "Separator Interface"],
     size = (1280, 720))
