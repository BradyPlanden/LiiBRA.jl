using LIBRA, Plots, BenchmarkTools, Infiltrator, MAT, StatsBase

#---------- Cell Definition -----------------#
Cell = Construct("LG_M50") #Alternative "Doyle_94"


#---------- Experimental Data -----------------# 75.7
SOC_Exp = [97.3  94.8  92.3  89.8  87.3  84.8  82.3  79.8  77.2  75  72.2  69.7  67.2  64.7  62.2  59.7  57.2  54.7  52.2  49.7  47.2  44.7  42.2  39.7  37.2  34.7  32.2  29.7  27.2  24.7  22.2  19.7  17.2  14.7  12.2  9.7  7.2  4.7  2.2  0.0] ./100

#---------- DRA Loop -----------------#
@inline function DRA_loop(Cell)
    A = B = C = D = Time = x =tuple()
    for Temp in 25.0:25.0:25.0
        Cell.Const.T = 273.15+Temp
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
            #for i in 3600*4:((3600*4*8)-(3600*4)):(8*4*3600)
                for k in 2500:2500:2500
                    #@show Cell.RA.Tlen = i
                    @show Cell.RA.H1 = 0:k
                    Cell.RA.H2 = 0:k
                    #Cell.RA.Tlen = i
                    #@show Cell.RA.M = k
                    Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
                    Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
                    Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)
                    #Cell.RA.Tlen = 15000+5000*i
                    #Cell.RA.Fs = (i)
                    #Cell.RA.SamplingT = 1/(i)

                    #for SOC in SOC_Exp[13:13]
                        #@show Cell.Const.SOC = SOC
                        A_DRA, B_DRA, C_DRA, D_DRA = DRA(Cell,Cell.RA.s,Cell.RA.f)
                        #x = @benchmark DRA(Cell,Cell.RA.s,Cell.RA.f)
                          A = flatten_(A,A_DRA)
                          B = flatten_(B,B_DRA)
                          C = flatten_(C,C_DRA)
                          D = flatten_(D,D_DRA)
                        #Time = flatten_(Time,x)
                    #end
                end
            #end
    end
return A, B, C, D#,Time 
end


#---------- Generate Model -----------------#
A, B, C, D = DRA_loop(Cell)
#Time = DRA_loop(Cell)

#---------- Simulate Model -----------------#
function Sim_loop(Cell, SOC_Exp, HPPC_Data)
 CellV = tDra = Ce = jNeg = jPos = RtotNeg = RtotPos = η0 = ηL = η_neg = η_pos = ϕ_ẽ1 = ϕ_ẽ2 = Uocp_Neg = Uocp_Pos = ϕ_e = Cse_Neg = Cse_Pos = tuple()
    for i in Int64(1/Cell.RA.SamplingT):Int64(1/Cell.RA.SamplingT):Int64(1/Cell.RA.SamplingT)
        Tk = ones(100*i+1)*298.15 #Cell Temperature
        Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i+2)*0.] #1C HPPC Experiment Current Profile
        Iapp = [ones(1)*0.; ones(50*i)*4.8181; ones(50*i)*-3.613;] #1C HPPC Experiment Current Profile
        Iapp = [Iapp;Iapp[2:end];Iapp[2:end];Iapp[2:end];Iapp[2:end];Iapp[2:end];Iapp[2:end];Iapp[2:end]]
        Tk = [Tk;Tk;Tk;Tk;Tk;Tk;Tk;Tk]
        time_ = 0:(1.0/i):(length(Iapp)*(1/i))
        #time_ = time_[1:end-1]
        tDra = flatten_(tDra,time_)
        #SOC = 0.2:0.1:0.9
        CellV_, Ce_, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg_, Uocp_Pos_, ϕ_e, Cse_Neg_, Cse_Pos_ = Sim_Model(Cell,Iapp,Tk,SOC_Exp[10],A,B,C,D)
        CellV = flatten_(CellV,CellV_)
        Ce = flatten_(Ce,Ce_)
        Cse_Neg = flatten_(Cse_Neg,Cse_Neg_)
        Cse_Pos = flatten_(Cse_Pos,Cse_Pos_)
        Uocp_Neg = flatten_(Uocp_Neg,Uocp_Neg_)
        Uocp_Pos = flatten_(Uocp_Pos,Uocp_Pos_)
    end
    return CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra
end

function RMSE(y1,y2)
    return sqrt(sum((y1.-y2).^2)/size(y1,1))
end

function Stats(k)
    for i in 1:k
        Rms_V[i] = rmsd(CellV[i][1:end-1],Pyb_V) #mol

        Rms_Cn[i] = rmsd(Cse_Neg[i][1:end-1,end],Pyb_Cn[20,20,:]) #mol
        Max_Cn[i] = maxad(Cse_Neg[i][1:end-1,end],Pyb_Cn[20,20,:]) #mol

        Rms_Cp[i] = rmsd(Cse_Pos[i][1:end-1,end],Pyb_Cp[20,20,:]) #mol
        Max_Cp[i] = maxad(Cse_Pos[i][1:end-1,end],Pyb_Cp[20,20,:]) #mol
        
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

HPPC_Data = HPPC_Data_Import(10,data_all) 
CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra = Sim_loop(Cell, SOC_Exp, Pyb_T)
#CellV, time, Uocp_Neg, Uocp_Pos = Sim_loop(Cell, [0.747])
#Stats(tuple_len(C))


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