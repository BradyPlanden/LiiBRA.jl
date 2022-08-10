using BenchmarkTools, LiiBRA, MAT, StatsBase, Plots

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 30
Cell = Construct("LG M50")
SList = collect(0.8)
SOC = 1.0
T = 298.15

function DRA_Loop(Cell, SList::Array, T::Float64)
    A = B = C = D = Time = x = tuple()
    for i in 4000:500:4000
        
        #Arrhenius
        Cell.Const.T = T
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R
        
        Cell.RA.H1 = 1:i
        Cell.RA.H2 = 1:i
        Cell.RA.Tlen = 16200
        Cell.RA.Fs = 4
        Cell.RA.M = 6
        Cell.RA.SamplingT = 1/4
        
        #Set Cell Constants
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)
        Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
        Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)
        
        # for Cell.Const.SOC in SList
        #     #Realisation
        #     A_DRA, B_DRA, C_DRA, D_DRA = CIDRA(Cell)
        #     A = flatten_(A,A_DRA)
        #     B = flatten_(B,B_DRA)
        #     C = flatten_(C,C_DRA)
        #     D = flatten_(D,D_DRA)
        # end
         x = @benchmark CIDRA(Cell)
         Time = flatten_(Time,x)
    
    end
    return Time, A, B, C, D
end

#---------- Experimental Data -----------------#

function Sim_loop(Cell, Input, SList, SOC, A, B, C, D)
    CellV = tDra = Ce = jNeg = jPos = RtotNeg = RtotPos = η0 = ηL = η_neg = η_pos = ϕ_ẽ1 = ϕ_ẽ2 = Uocp_Neg = Uocp_Pos = ϕ_e = Cse_Neg = Cse_Pos = Cell_SOC = SimTime = x = tuple()
       for i in Int64(1/Cell.RA.SamplingT):Int64(1/Cell.RA.SamplingT):Int64(1/Cell.RA.SamplingT)
           Tk = ones(length(Input))*298.15 #Cell Temperature
           time_ = 0:(1.0/i):(length(Input)*(1/i))
           tDra = flatten_(tDra,time_)

        #    u = 1
        #    k = 0
        #    for i in 1:Int64(tuple_len(C)/length(SList))

        #         k += length(SList)
                Cell_V, Ce_, jNeg, jPos, Rtot_neg, Rtot_pos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg_, Uocp_Pos_, ϕ_e, Cse_Neg_, Cse_Pos_, Cell_SOC, jeq_neg, jeq_pos, j0, jL, tDra = Sim_Model(Cell,Input,"Power",Tk,SList,SOC,A,B,C,D,tDra)
        #        u = k+1

                CellV = flatten_(CellV,Cell_V)
                Ce = flatten_(Ce,Ce_)
                Cse_Neg = flatten_(Cse_Neg,Cse_Neg_)
                Cse_Pos = flatten_(Cse_Pos,Cse_Pos_)
                Uocp_Neg = flatten_(Uocp_Neg,Uocp_Neg_)
                Uocp_Pos = flatten_(Uocp_Pos,Uocp_Pos_)

        #   end
           #for i in 1:12
           # @time Sim_Model(Cell,Input,"Power",Tk,SList,SOC,A,B,C,D)
           #end
       end
       return CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra, Cell_SOC, SimTime
   end
   
   function Stats()
       for i in 1:tuple_len(CellV)
           Rms_V[i] = rmsd(CellV[i],Pyb_V) # V
           Max_V[i] = maxad(CellV[i],Pyb_V) # V
   
           Rms_Cn[i] = rmsd(Cse_Neg[i][:,end],Pyb_Cn[20,20,:]) #mol
           Max_Cn[i] = maxad(Cse_Neg[i][:,end],Pyb_Cn[20,20,:]) #mol
   
           Rms_Cp[i] = rmsd(Cse_Pos[i][:,end],Pyb_Cp[20,20,:]) #mol
           Max_Cp[i] = maxad(Cse_Pos[i][:,end],Pyb_Cp[20,20,:]) #mol
           
       end
   return Rms_Cn, Max_Cn, Rms_Cp, Max_Cp, Rms_V, Max_V
   end
   

function Sensitivity(Rng,SList,WLTP_P,SOC)
    k=1
    Output = ones(length(Rng)^2,3)
    for i ∈ Rng
        for j ∈ Rng
            @show Cell.Pos.L = i
            @show Cell.Neg.L = j
            Time, A, B, C, D = DRA_Loop(Cell, SList, T)
            #Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i+1)*0.] #1C HPPC Experiment Current Profile
            CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e, Cse_Neg, Cse_Pos, tDra, Cell_SOC, SimTime = Sim_loop(Cell, WLTP_P, SList, SOC, A, B, C, D)
            Output[k,:] = [mean(CellV[1]), i, j]
            k+=1
        end
    end
    return Output
end

WLTP_File = matopen("test/WLTP/WLTP_M50_M3.mat")
WLTP_P = read(WLTP_File,"P_Models")
Time, A, B, C, D = DRA_Loop(Cell, SList, T)
#Output = Sensitivity(20e-6:10e-6:100e-6,SList,WLTP_P,SOC)

#    Pyb_File = matopen("test/WLTP/sol_data.mat")#
#    Pyb_Cn = read(Pyb_File,"c_n")
#    Pyb_Cp = read(Pyb_File,"c_p")
#    Pyb_T = read(Pyb_File,"t")
#    Pyb_V = read(Pyb_File,"V")
#    Rms_Cn = Array{Float64}(undef,tuple_len(CellV),1)
#    Max_Cn = Array{Float64}(undef,tuple_len(CellV),1)
#    Rms_Cp = Array{Float64}(undef,tuple_len(CellV),1)
#    Max_Cp = Array{Float64}(undef,tuple_len(CellV),1)
#    Rms_V = Array{Float64}(undef,tuple_len(CellV),1)
#    Max_V = Array{Float64}(undef,tuple_len(CellV),1)
   
   #Rms_Cn, Max_Cn, Rms_Cp, Max_Cp, Rms_V = Stats()