using LiBRA, Plots, BenchmarkTools, Infiltrator, MAT, StatsBase

#---------- Cell Definition -----------------#
Cell = Construct("LG_M50") #Alternative "Doyle_94"
Init_SOC = 0.75

#---------- DRA Loop -----------------#
@inline function DRA_loop(Cell, Init_SOC)
    A = B = C = D = tuple()

    for Temp in 5.0:10.0:45.0
        #Arrhenius
        Cell.Const.T = 273.15+Temp
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R

        #Set Cell Constants
        Cell.Const.SOC = Init_SOC
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)

        #DRA
        A_DRA, B_DRA, C_DRA, D_DRA = DRA(Cell,Cell.RA.s,Cell.RA.f)

        #Flatten output into tuples
        A = flatten_(A,A_DRA)
        B = flatten_(B,B_DRA)
        C = flatten_(C,C_DRA)
        D = flatten_(D,D_DRA)
    end
return A, B, C, D
end


#---------- Generate Model -----------------#
A, B, C, D = DRA_loop(Cell, Init_SOC)

#---------- Sim Loop -----------------#
function Sim_loop(Cell, Init_SOC)
    CellV = Ce = jNeg = jPos = RtotNeg = RtotPos = η0 = ηL = η_neg = η_pos = ϕ_ẽ1 = ϕ_ẽ2 = Uocp_Neg = Uocp_Pos = ϕ_e = Cse_Neg = Cse_Pos = tuple()

    #Set Experiment
    i = Int64(1/Cell.RA.SamplingT) #Sampling Frequency
    Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i)*0.] #1C HPPC Experiment Current Profile
    Tk = ones(size(Iapp))*298.15 #Cell Temperature
    tDra = 0:(1.0/i):((length(Iapp)-1)/i)
    
    #Simulate Model
    CellV_, Ce_, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg_, Uocp_Pos_, ϕ_e, Cse_Neg_, Cse_Pos_ = Sim_Model(Cell,Iapp,Tk,Init_SOC,A,B,C,D)
    
    #Flatten output into tuples
    CellV = flatten_(CellV,CellV_)
    Ce = flatten_(Ce,Ce_)
    Cse_Neg = flatten_(Cse_Neg,Cse_Neg_)
    Cse_Pos = flatten_(Cse_Pos,Cse_Pos_)
    Uocp_Neg = flatten_(Uocp_Neg,Uocp_Neg_)
    Uocp_Pos = flatten_(Uocp_Pos,Uocp_Pos_)

return CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra
end

#---------- Simulate Model -----------------#
CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra = Sim_loop(Cell, Init_SOC)

#----------- Plotting ---------------------------#
plotly()
display(plot(tDra[1:end-1],CellV[1][1:end-1], legend=:topright,color=:blue,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Terminal Voltage [V]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Ce[1][1:end-1,[1,2,19,20]], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Electrolyte Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Cse_Pos[1][1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Pos. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Cse_Neg[1][1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Neg. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))