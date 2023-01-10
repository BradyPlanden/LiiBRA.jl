using LiiBRA, PGFPlotsX, Colors, LaTeXStrings

#---------- Cell Definition -----------------#
Cell = Construct("LG M50")
Cell.RA.Fs = 2.0 # Modify transfer function sampling frequency
Cell.RA.SamplingT = 0.5 # Modify final system sampling period
Cell.Const.T = 298.15 # Modify cell temperature
Cell.RA.H1 = Cell.RA.H2 = [1:2500; 3000:3500; 4000:4500]
Cell.RA.M = 5
Ŝ = collect(1.0:-1.0:0) # List of SOC points for model generation
SOC = 1. # Starting SOC
α = collect(0.95:-0.1:0.75)

function DegLp(α)
Results = tₑₜ = tuple()
    for ψ ∈ α
        Cell.Neg.θ_100 = ψ
        A, B, C, D = Realise(Cell, Ŝ)
        Result, tₑ = CC(Cell,Ŝ,SOC,2.0,9600,A,B,C,D)
        Results = flatten_(Results,Result)
        tₑₜ = flatten_(tₑₜ,tₑ)
    end
    return Results, tₑₜ
end

Results, tₑₜ = DegLp(α)

#----------- Plotting ---------------------------#
κ = distinguishable_colors(10)
ax1 = @pgf Axis(
    {
        ylabel = "Voltage",
        clip_mode="individual",
    }, 
    Plot({color=κ[10], "thick"},
        Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Cell_V[1:tₑₜ[1]]])
        ), LegendEntry(L"$\theta_n^0$ = 0.95"),
    Plot({color=κ[7], "thick"},
        Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Cell_V[1:tₑₜ[2]]])
        ), LegendEntry(L"$\theta_n^0$ = 0.85"),    
    Plot({color=κ[8], "thick"},
        Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Cell_V[1:tₑₜ[3]]])
        ), LegendEntry(L"$\theta_n^0$ = 0.75"),
    )

    ax2 = @pgf Axis(
    {
        ylabel = "Electrolyte Conc. (mol/m³)",
        clip_mode="individual",
        scaled_y_ticks={"base 10:-3"},

    },
    Plot({color=κ[10], "thick"},
        Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Ce[1:tₑₜ[1],1]])
        ),
    Plot({color=κ[10], "thick"},
        Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Ce[1:tₑₜ[1],2]])
        ),
    Plot({color=κ[10], "thick"},
        Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Ce[1:tₑₜ[1],3]])
        ),
    Plot({color=κ[10], "thick"},
        Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Ce[1:tₑₜ[1],4]])
        ),
    Plot({color=κ[7], "thick"},
        Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Ce[1:tₑₜ[2],1]])
        ),
    Plot({color=κ[7], "thick"},
        Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Ce[1:tₑₜ[2],2]])
        ),
    Plot({color=κ[7], "thick"},
        Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Ce[1:tₑₜ[2],3]])
        ),
    Plot({color=κ[7], "thick"},
        Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Ce[1:tₑₜ[2],4]])
        ),  
    Plot({color=κ[8], "thick"},
        Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Ce[1:tₑₜ[3],1]])
        ),
    Plot({color=κ[8], "thick"},
        Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Ce[1:tₑₜ[3],2]])
        ),
    Plot({color=κ[8], "thick"},
        Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Ce[1:tₑₜ[3],3]])
        ),
    Plot({color=κ[8], "thick"},
        Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Ce[1:tₑₜ[3],4]])
        ),  
    )

    ax3 = @pgf Axis(
    {
        ylabel = "Positive Electrode Conc. (mol/m³)",
        clip_mode="individual",

    },
    Plot({color=κ[10], "thick"},
        Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Cse_Pos[1:tₑₜ[1]]])
        ),
    Plot({color=κ[7], "thick"},
        Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Cse_Pos[1:tₑₜ[2]]])
        ),    
    Plot({color=κ[8], "thick"},
        Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Cse_Pos[1:tₑₜ[3]]])
        ),
    )

    ax4 = @pgf Axis(
        {
            ylabel = "Negative Electrode Conc. (mol/m³)",
            clip_mode="individual",

        },
        Plot({color=κ[10], "thick"},
            Table([:x => Results[1].t[1:tₑₜ[1]], :y =>Results[1].Cse_Neg[1:tₑₜ[1]]])
            ),
        Plot({color=κ[7], "thick"},
            Table([:x => Results[2].t[1:tₑₜ[2]], :y =>Results[2].Cse_Neg[1:tₑₜ[2]]])
            ),    
        Plot({color=κ[8], "thick"},
            Table([:x => Results[3].t[1:tₑₜ[3]], :y =>Results[3].Cse_Neg[1:tₑₜ[3]]])
            ),
        )

@pgf GroupPlot(
    {
        group_style = {group_size = "2 by 2", horizontal_sep = "2cm", vertical_sep = "1.5cm"},
        no_markers,
        legend_pos = "north east",
        xlabel = "Time (s)",

    },
    ax1, ax2, ax3, ax4

)