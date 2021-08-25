using LIBRA, Plots

#Cell Definition
#CellTyp = "Doyle_94"

@inline function Impulse() # Create s Vector 
    Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    #Nfft = ceil(2^(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    f = 0:Nfft-1
    s = transpose(((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)))
    return Nfft,f,s
end

@inline function DRA_loop(CellData)
    A_DRA = B_DRA = C_DRA = D_DRA = Dtt = puls = Hank1 = Hank2 = S = U = V = tf__ = tuple()
    #A_DRA = B_DRA = C_DRA = D_DRA = Dtt = tuple()
    for Temp in 25.0:25.0:25.0
        CellData.Const.T = 273.15+Temp
        Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/R
        CellData.Const.κ = CellData.Const.κf(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
            for i in 1:1
                Nfft = ceil(2^(log2(CellData.RA.Fs*CellData.RA.Tlen)))
                #Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
                f = 0:Nfft-1
                s = transpose(((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)))
                for SOC in 0.8:0.8:0.8
                    CellData.Const.SOC = SOC
                    A, B, C, D, Dtt, puls, Hank1, Hank2, S, U, V, tf__ = DRA(CellData,s,f)
                    #A, B, C, D, Dtt = DRA(CellData,s,f)
                    A_DRA = flatten(A_DRA,A)
                    B_DRA = flatten(B_DRA,B)
                    C_DRA = flatten(C_DRA,C)
                    D_DRA = flatten(D_DRA,D)
                end
            end
    end

return A_DRA, B_DRA, C_DRA, D_DRA, Dtt, puls, Hank1, Hank2, S, U, V, tf__
end

#----------Generate Model -----------------#
#TransferFuns = TransferFun()
#Nfft, f, s = Impulse()
#A_DRA, B_DRA, C_DRA, D_DRA, Dtt = DRA_loop(CellData)
A_DRA, B_DRA, C_DRA, D_DRA, Dtt, puls, Hank1, Hank2, S, U, V, tf__ = DRA_loop(CellData)
#save("$CellTyp.jld", "CellData", CellData, "A_DRA", A_DRA, "B_DRA", B_DRA, "C_DRA", C_DRA, "D_DRA", D_DRA, "Dtt", Dtt) #Switch to jld2


#----------Simulate Model -----------------#
Tk = ones(220)*298.15
Iapp = [ones(20)*0.; ones(20)*9.45; ones(80)*0.; ones(20)*-5; ones(80)*0.]
time = 0:0.5:109
CellV, jNeg, jPos, y, x, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, j0, jL, j0_CC_neg, j0_CC_pos, Uocp_Neg, Uocp_Pos, Cse_Neg, Cse_Pos, Ce = Sim_Model(CellData,Dtt,Iapp,Tk,A_DRA,B_DRA,C_DRA,D_DRA)

plot(time,CellV[1:end-1], legend=:topright,color=:blue,label="Voltage",bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Voltage [V]")
plot!(twinx(),time,Iapp[1:end-1],legend=:bottomright,color=:red,label="Current", left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Current [A]")
plot!(title = "LG M50 - HPPC @ 90% SOC", xlabel = "Time [s]")