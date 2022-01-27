using LiiBRA, Plots

#---------- Cell Definition -----------------#
Cell = Construct("LG_M50") #Alternative "Doyle_94"
SList = collect(0.85:-0.05:0.65)
SOC = 0.75

#---------- Generate Models -----------------#
function Realise(Cell, SList)
    A = B = C = D = tuple()
    for i in SList 
        #Arrhenius
        Cell.Const.T = 298.15
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R

        #Set Cell Constants
        Cell.Const.SOC = i
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)

        #Realisation
        Aϕ, Bϕ, Cϕ, Dϕ = DRA(Cell,Cell.RA.s,Cell.RA.f)

        #Flatten output into Tuples
        A = flatten_(A,Aϕ)
        B = flatten_(B,Bϕ)
        C = flatten_(C,Cϕ)
        D = flatten_(D,Dϕ)
    end
return A, B, C, D
end


#---------- Simulate Models -----------------#
function Simulate(Cell,SList,Init_SOC,A,B,C,D)

    #Set Experiment
    i = Int64(1/Cell.RA.SamplingT) #Sampling Frequency
    Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i)*0.] #1C HPPC Experiment Current Profile
    Tk = ones(size(Iapp))*298.15 #Cell Temperature
    tDra = 0:(1.0/i):((length(Iapp)-1)/i)
    
    #Simulate Model
    return Sim_Model(Cell,Iapp,Tk,SList,Init_SOC,A,B,C,D,tDra)
end

#---------- Generate Model -----------------#
A, B, C, D = Realise(Cell, SList)


#---------- Simulate Model -----------------#
Cell_V, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, Cell_SOC, tDra, j0, jL = Simulate(Cell,SList,Init_SOC,A,B,C,D)

#----------- Plotting ---------------------------#
plotly()
display(plot(tDra[1:end-1],Cell_V[1:end-1], legend=:topright,color=:blue,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Terminal Voltage [V]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Ce[1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Electrolyte Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Cse_Pos[1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Pos. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
display(plot(tDra[1:end-1],Cse_Neg[1:end-1,:], legend=:topright,bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Neg. Electrode Concen. [mol/m^3]", xlabel = "Time [s]"))
