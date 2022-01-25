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

Install (Julia 1.6 and later)
-----------------------------

```julia
(v1.6) pkg> add LiiBRA
```

(Type `]` to enter package mode.)

## Example Usage

```julia
using LiiBRA, Plots
```

Construct cell parameterisation, and define SOC points:
```julia
Cell = Construct("LG_M50") #Alternative "Doyle_94"
SList = collect(0.85:-0.05:0.65)
Init_SOC = 0.75
```

Function Call for Realisation:
```julia
function DRA_loop(Cell, SList)
    A = B = C = D = tuple()
    for i in SList 
        #Arrhenius
        Cell.Const.T = 298.15
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R

        #Set Cell Constants
        Cell.Const.SOC = i #Init_SOC
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)

        #Realisation
        A_DRA, B_DRA, C_DRA, D_DRA = DRA(Cell,Cell.RA.s,Cell.RA.f)

        #Flatten output into Tuples
        A = flatten_(A,A_DRA)
        B = flatten_(B,B_DRA)
        C = flatten_(C,C_DRA)
        D = flatten_(D,D_DRA)
    end
return A, B, C, D
end
```

Call Realisation:
```julia
A, B, C, D = DRA_loop(Cell, SList)
```


Function Call for Simulation:
```julia
function Sim_loop(Cell,SList,Init_SOC,A,B,C,D)
    #Set Experiment
    i = Int64(1/Cell.RA.SamplingT) #Sampling Frequency
    Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i)*0.] #1C HPPC Experiment Current Profile
    Tk = ones(size(Iapp))*298.15 #Cell Temperature
    tDra = 0:(1.0/i):((length(Iapp)-1)/i)
    
    #Simulate Model
    return Cell_V, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, Cell_SOC, tDra = Sim_Model(Cell,Iapp,Tk,SList,Init_SOC,A,B,C,D,tDra)
end
```
Simulate Model:
```julia
Cell_V, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, Cell_SOC, tDra, j0, jL = Sim_loop(Cell,SList,Init_SOC,A,B,C,D)
```

Plot Results:
```julia
plotly()
display(plot(tDra[1:end-1],Cell_V[1:end-1], legend=:topright,color=:blue,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Terminal Voltage [V]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Ce[1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Electrolyte Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Cse_Pos[1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Pos. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Cse_Neg[1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Neg. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
```

## Bug Tracking

Please report any issues using the Github [issue tracker]. All feedback is welcome.

[issue tracker]: https://github.com/BradyPlanden/LiiBRA/issues
