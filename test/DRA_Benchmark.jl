using BenchmarkTools, LiiBRA, MAT, StatsBase, Plots  # run PROPACK, tsvd (remove k=) then change to svds and do ARPACK

#include("/Users/bradyplanden/Documents/Git/LIBRA_Paper/Data/Experimental/HPPC/HPPCStepCalc.jl")
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
Cell = Construct("LG M50")
SList = collect(0.9:-0.05:0.55)
SOC = 0.75
T = 298.15

function DRA_Loop(Cell, SList::Array, T::Float64)
    A = B = C = D = Time = x = tuple()
    for i in 1500:3000:4500
        
        #Arrhenius
        Cell.Const.T = T
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R
        
        Cell.RA.H1 = 0:i
        Cell.RA.H2 = 0:i
        Cell.RA.Tlen = 16200
        Cell.RA.Fs = 20
        Cell.RA.M = 9
        Cell.RA.SamplingT = 1/10
        
        #Set Cell Constants
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)
        Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
        Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)
        
        for Cell.Const.SOC in SList
            #Realisation
            A_DRA, B_DRA, C_DRA, D_DRA = DRA(Cell,Cell.RA.s,Cell.RA.f)
            A = flatten_(A,A_DRA)
            B = flatten_(B,B_DRA)
            C = flatten_(C,C_DRA)
            D = flatten_(D,D_DRA)
        end
        x = @benchmark DRA(Cell,Cell.RA.s,Cell.RA.f)
        Time = flatten_(Time,x)
    
    end
    return Time, A, B, C, D
end

#---------- Experimental Data -----------------#

function Sim_loop(Cell, Input, SList, SOC, A, B, C, D)
    CellV = tDra = Ce = jNeg = jPos = RtotNeg = RtotPos = η0 = ηL = η_neg = η_pos = ϕ_ẽ1 = ϕ_ẽ2 = Uocp_Neg = Uocp_Pos = ϕ_e = Cse_Neg = Cse_Pos = tuple()
       for i in Int64(1/Cell.RA.SamplingT):Int64(1/Cell.RA.SamplingT):Int64(1/Cell.RA.SamplingT)
           Tk = ones(length(Input))*298.15 #Cell Temperature
           time_ = 0:(1.0/i):(length(Input)*(1/i))
           tDra = flatten_(tDra,time_)
        
           u = 1
           k = 0
           for i in 1:Int64(tuple_len(C)/length(SList))
                k += length(SList)
                Cell_V, Ce_, jNeg, jPos, Rtot_neg, Rtot_pos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg_, Uocp_Pos_, ϕ_e, Cse_Neg_, Cse_Pos_, Cell_SOC, tDra, jeq_neg, jeq_pos, j0, jL = Sim_Model(Cell,Input,"Power",Tk,SList,SOC,A[u:k],B[u:k],C[u:k],D[u:k],tDra)
                u = k+1

                CellV = flatten_(CellV,Cell_V)
                Ce = flatten_(Ce,Ce_)
                Cse_Neg = flatten_(Cse_Neg,Cse_Neg_)
                Cse_Pos = flatten_(Cse_Pos,Cse_Pos_)
                Uocp_Neg = flatten_(Uocp_Neg,Uocp_Neg_)
                Uocp_Pos = flatten_(Uocp_Pos,Uocp_Pos_)

           end
       end
       return CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra
   end
   
   function Stats()
       for i in 1:tuple_len(CellV)
           Rms_V[i] = rmsd(CellV[i],Pyb_V) #mol
   
           Rms_Cn[i] = rmsd(Cse_Neg[i][:,end],Pyb_Cn[20,20,:]) #mol
           Max_Cn[i] = maxad(Cse_Neg[i][:,end],Pyb_Cn[20,20,:]) #mol
   
           Rms_Cp[i] = rmsd(Cse_Pos[i][:,end],Pyb_Cp[20,20,:]) #mol
           Max_Cp[i] = maxad(Cse_Pos[i][:,end],Pyb_Cp[20,20,:]) #mol
           
       end
   return Rms_Cn, Max_Cn, Rms_Cp, Max_Cp, Rms_V
   end
   
   Time, A, B, C, D = DRA_Loop(Cell, SList, T)
   HPPC_Data = HPPC_Data_Import(10,data_all) 
   WLTP_File = matopen("/Users/bradyplanden/Documents/Git/LiiBRA.jl/test/WLTP_M50_M3.mat")
   WLTP_P = read(WLTP_File,"WLTP_M50_M3")
   #Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i+1)*0.] #1C HPPC Experiment Current Profile
   CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra = Sim_loop(Cell, WLTP_P, SList, SOC, A, B, C, D)
   
   Pyb_File = matopen("/Users/bradyplanden/Documents/Git/LIBRA_Paper/Data/PyBaMM/sol_data.mat")
   Pyb_Cn = read(Pyb_File,"c_n")
   Pyb_Cp = read(Pyb_File,"c_p")
   Pyb_T = read(Pyb_File,"t")
   Pyb_V = read(Pyb_File,"V")
   Rms_Cn = Array{Float64}(undef,tuple_len(CellV),1)
   Max_Cn = Array{Float64}(undef,tuple_len(CellV),1)
   Rms_Cp = Array{Float64}(undef,tuple_len(CellV),1)
   Max_Cp = Array{Float64}(undef,tuple_len(CellV),1)
   Rms_V = Array{Float64}(undef,tuple_len(CellV),1)
   
   Rms_Cn, Max_Cn, Rms_Cp, Max_Cp, Rms_V = Stats()