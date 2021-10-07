using LIBRA, Plots

#---------- Cell Definition -----------------#
Cell = cell("LG_M50") #Alternative "Doyle_94"


#---------- Experimental Data -----------------#
SOC_Exp = [97.3  94.8  92.3  89.8  87.3  84.8  82.3  79.8  77.2  74.7  72.2  69.7  67.2  64.7  62.2  59.7  57.2  54.7  52.2  49.7  47.2  44.7  42.2  39.7  37.2  34.7  32.2  29.7  27.2  24.7  22.2  19.7  17.2  14.7  12.2  9.7  7.2  4.7  2.2  0.0] ./100

#---------- DRA Loop -----------------#
@inline function DRA_loop(Cell)
    A = B = C = D = tuple()
    for Temp in 25.0:25.0:25.0
        Cell.Const.T = 273.15+Temp
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
            for i in 1:1:1
                #Cell.RA.H1 = 0:i
                #Cell.RA.H2 = 0:i
                #Cell.RA.Tlen = 15000+5000*i
                #Cell.RA.Fs = (i)
                #Cell.RA.SamplingT = 1/(i)
                Nfft = 2^(ceil(log2(Cell.RA.Fs*Cell.RA.Tlen)))
                f = 0:Nfft-1
                s = transpose(((2im.*Cell.RA.Fs)*tan.(pi.*f./Nfft)))
                for SOC in SOC_Exp[13:16]
                    @show Cell.Const.SOC = SOC
                    A_DRA, B_DRA, C_DRA, D_DRA = DRA(Cell,s,f)
                    A = flatten_(A,A_DRA)
                    B = flatten_(B,B_DRA)
                    C = flatten_(C,C_DRA)
                    D = flatten_(D,D_DRA)
                end
            end
    end
return A, B, C, D
end


#---------- Generate Model -----------------#
A, B, C, D = DRA_loop(Cell)


#---------- Simulate Model -----------------#
function Sim_loop(Cell, SOC_Exp)
 CellV = time = Ce = jNeg = jPos = RtotNeg = RtotPos = η0 = ηL = η_neg = η_pos = ϕ_ẽ1 = ϕ_ẽ2 = Uocp_Neg = Uocp_Pos = ϕ_e = Cse_Neg = Cse_Pos = tuple()
    for i in 4:4:4
        Tk = ones(100*i+1)*298.15 #Cell Temperature
        Iapp = [ones(1)*0.; ones(10*i)*4.8181; ones(40*i)*0.; ones(10*i)*-3.613; ones(40*i)*0.] #1C HPPC Experiment Current Profile
        time_ = 0:(1.0/i):100
        time_ = time_[1:end-1]
        time = flatten_(time,time_)
        #SOC = 0.2:0.1:0.9
        CellV_, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg_, Uocp_Pos_, ϕ_e, Cse_Neg, Cse_Pos = Sim_Model(Cell,Iapp,Tk,SOC_Exp[13:16],A,B,C,D)
        CellV = flatten_(CellV,CellV_)
        Uocp_Neg = flatten_(Uocp_Neg,Uocp_Neg_)
        Uocp_Pos = flatten_(Uocp_Pos,Uocp_Pos_)
    end
    return CellV, time, Uocp_Neg, Uocp_Pos
end

CellV, time, Uocp_Neg, Uocp_Pos = Sim_loop(Cell, SOC_Exp)
#----------- Plotting ---------------------------#
# plotly()
# display(plot(time,CellV[1:end-1], legend=:topright,color=:blue,label="Voltage",bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Voltage [V]"))
# display(plot!(twinx(),time,Iapp[1:end-1],legend=:bottomright,color=:red,label="Current", left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Current [A]"))
# display(plot!(title = "LG M50 - HPPC @ 75% SOC", xlabel = "Time [s]"))

# display(plot(Ce))
# display(plot(ϕ_e))
# display(plot(Cse_Pos[1:end-1,:]))
# display(plot!(Cse_Neg[1:end-1,:]))