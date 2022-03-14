# Lithium-ion Battery Realisation Algorithms (LiiBRA)

[![Build Status](https://github.com/BradyPlanden/LiiBRA.jl/workflows/CI/badge.svg)](https://github.com/BradyPlanden/LiiBRA.jl/actions)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

<p align="center">
<img src="LiiBRA.png" width="600" align="center"  />
</p>

## Create and Simulation Reduced Order Lithium-Ion Battery Models
LiBRA provides an open-source implementation of realisation algorithms commonly used for generating reduced-order state-space models.
For more information on the current derived transfer functions used to represent the lithium-ion cell, please refer to [1](https://doi.org/10.1016/j.enconman.2007.03.015), [2](https://doi.org/10.1016/j.jpowsour.2012.07.075).

Install (Julia 1.7 and later)
-----------------------------

```julia
(v1.7) pkg> add LiiBRA
```

(Type `]` to enter package mode.)

## Example Usage

```julia
using LiiBRA, Plots
```

Setup:
```julia
Cell = Construct("LG M50")
SList = collect(1.0:-0.25:0.0)
SOC = 0.75
T = 298.15
```

Realisation:
```julia
A,B,C,D = Realise(Cell,SList,T)
```

HPPC Simulation:
```julia
CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, Cell_SOC, tDra, j0, jL, tDra = HPPC(Cell,SList,SOC,4.0,-3.0,A,B,C,D)
```

Plot Results:
```julia
plotly()
display(plot(tDra,CellV, legend=:topright,color=:blue,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Terminal Voltage [V]", xlabel = "Time [s]"))
display(plot(tDra,Ce, legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Electrolyte Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra,Cse_Pos, legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Pos. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra,Cse_Neg, legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Neg. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
```

## Bug Tracking

Please report any issues using the Github [issue tracker]. All feedback is welcome.

[issue tracker]: https://github.com/BradyPlanden/LiiBRA/issues
