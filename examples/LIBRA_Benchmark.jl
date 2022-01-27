using LiiBRA, Plots, BenchmarkTools, MAT, StatsBase

#---------- Cell Definition -----------------#
Cell = Construct("LG_M50") #Alternative "Doyle_94"
SList = collect(0.85:-0.05:0.65)
Init_SOC = 0.75

#---------- Experimental Data -----------------# 75.7
#SList = [97.3  94.8  92.3  89.8  87.3  84.8  82.3  79.8  77.2  75  72.2  69.7  67.2  64.7  62.2  59.7  57.2  54.7  52.2  49.7  47.2  44.7  42.2  39.7  37.2  34.7  32.2  29.7  27.2  24.7  22.2  19.7  17.2  14.7  12.2  9.7  7.2  4.7  2.2  0.0] ./100


#---------- DRA Loop -----------------#
function DRA_loop(Cell, SList)
    A = B = C = D = x = Time = tuple()
    for i in SList 
        #Arrhenius
        Cell.Const.T = 298.15
        Cell.RA.Tlen = 65536
        Cell.RA.H1 = 0:3000
        Cell.RA.H2 = 0:3000
        Cell.RA.Fs = 6
        #Cell.Pos.σ = 1.5
        #Cell.Pos.Ds = 6e-15 
        #Cell.Neg.k_norm = 7.264265E-07
        #Cell.Pos.k_norm = 7.26426E-06
        #Cell.Pos.Rs = 5.5e-6
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
        #x = @benchmark DRA(Cell,Cell.RA.s,Cell.RA.f)
        A = flatten_(A,A_DRA)
        B = flatten_(B,B_DRA)
        C = flatten_(C,C_DRA)
        D = flatten_(D,D_DRA)
        #Time = flatten_(Time,x)

    end
return A, B, C, D#,Time
end



#---------- Generate Model -----------------#
A, B, C, D = DRA_loop(Cell, SList)
#Time = DRA_loop(Cell)

#---------- Simulate Model -----------------#
function Sim_loop(Cell,SList,Init_SOC,A,B,C,D)

    #Set Experiment
    i = Int64(1/Cell.RA.SamplingT) #Sampling Frequency
    #Iapp =[ones(1)*0; ones(900*i)*4.8181]
    Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i)*0.] #1C HPPC Experiment Current Profile
    Tk = ones(size(Iapp))*298.15 #Cell Temperature
    tDra = 0:(1.0/i):((length(Iapp)-1)/i)
    
    #Simulate Model
    return Cell_V, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, Cell_SOC, tDra, jeq_neg, jeq_pos, j0, jL, y = Sim_Model(Cell,Iapp,Tk,SList,Init_SOC,A,B,C,D,tDra)
end

function RMSE(y1,y2)
    return sqrt(sum((y1.-y2).^2)/size(y1,1))
end

function Stats(k)
    for i in 1:k
        Rms_V[i] = rmsd(Cell_V[1:end-1],Pyb_V[1:end-2]) #mol

        Rms_Cn[i] = rmsd(Cse_Neg[1:end-1,end],Pyb_Cn[20,20,1:end-2]) #mol
        Max_Cn[i] = maxad(Cse_Neg[1:end-1,end],Pyb_Cn[20,20,1:end-2]) #mol

        Rms_Cp[i] = rmsd(Cse_Pos[1:end-1,end],Pyb_Cp[20,20,1:end-2]) #mol
        Max_Cp[i] = maxad(Cse_Pos[1:end-1,end],Pyb_Cp[20,20,1:end-2]) #mol
        
    end
return Rms_Cn, Max_Cn, Rms_Cp, Max_Cp, Rms_V
end

file = matopen("/Users/bradyplanden/Documents/Git/LIBRA_Paper/Data/PyBaMM/sol_data.mat")
Pyb_Cn = read(file,"c_n")
Pyb_Cp = read(file,"c_p")
Pyb_T = read(file,"t")
Pyb_V = read(file,"V")
Rms_Cn = Array{Float64}(undef,tuple_len(C),1)
Max_Cn = Array{Float64}(undef,tuple_len(C),1)
Rms_Cp = Array{Float64}(undef,tuple_len(C),1)
Max_Cp = Array{Float64}(undef,tuple_len(C),1)
Rms_V = Array{Float64}(undef,tuple_len(C),1)

#HPPC_Data = HPPC_Data_Import(10,data_all) 
Cell_V, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, Cell_SOC, tDra, jeq_neg, jeq_pos, j0, jL, y = Sim_loop(Cell,SList,Init_SOC,A,B,C,D)
#CellV, time, Uocp_Neg, Uocp_Pos = Sim_loop(Cell, [0.747])
Stats(tuple_len(C))


#----------- Plotting ---------------------------#
# plotly()
# display(plot(tDra,CellV[1:end-1], legend=:topright,color=:blue,label="Voltage",bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Voltage [V]"))
# display(plot!(twinx(),HPPC_Data[:,1],HPPC_Data[:,3],legend=:bottomright,color=:red,label="Current", left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Current [A]"))
# display(plot!(title = "LG M50 - HPPC @ 75% SOC", xlabel = "Time [s]"))
# plot(tDra,CellV[1:end-1])
# plot!(HPPC_Data[:,1],HPPC_Data[:,3])
# display(plot(Ce))
# display(plot(ϕ_e))
# display(plot(Cse_Pos[1:end-1,:]))
# display(plot!(Cse_Neg[1:end-1,:]))