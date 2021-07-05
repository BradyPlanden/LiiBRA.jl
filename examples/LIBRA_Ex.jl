using LIBRA, JLD, Plots

#Cell Definition
#CellTyp = "Doyle_94"
#include("Data/Chen_2020/LG_M50.jl")

@inline function Impulse() # Create s Vector 
    Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    #Nfft = ceil(2^(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    f = 0:Nfft-1
    s = transpose(((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)))
    return Nfft,f,s
end

@inline function DRA_loop(CellData)
    A_DRA = B_DRA = C_DRA = D_DRA = Dtt = tf__ = s = tuple()
    for Temp in 25.0:25.0:25.0
        CellData.Const.T = 273.15+Temp
        Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/R
        CellData.Const.κ = CellData.Const.κf(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
            for i in 1:1
                #CellData.RA.Tlen += 2000 
                #Nfft = ceil(2^(log2(CellData.RA.Fs*CellData.RA.Tlen)))
                Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
                f = 0:Nfft-1
                s = transpose(((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)))
                #println("Tlen:",CellData.RA.Tlen)
                for SOC in 1.:1.:1.
                    CellData.Const.SOC = SOC
                    A, B, C, D, Dtt, tf__ = DRA(CellData,s,f)
                    A_DRA = flatten(A_DRA,A)
                    B_DRA = flatten(B_DRA,B)
                    C_DRA = flatten(C_DRA,C)
                    D_DRA = flatten(D_DRA,D)
                end
            end
    end

return A_DRA, B_DRA, C_DRA, D_DRA, Dtt, tf__, s
end

#TransferFuns = TransferFun()
#Nfft, f, s = Impulse()
#A_DRA, B_DRA, C_DRA, D_DRA, Dtt = DRA_loop(CellData,s,f)
A_DRA, B_DRA, C_DRA, D_DRA, Dtt, tf__, s = DRA_loop(CellData)
#save("$CellTyp.jld", "CellData", CellData, "A_DRA", A_DRA, "B_DRA", B_DRA, "C_DRA", C_DRA, "D_DRA", D_DRA, "Dtt", Dtt) #Switch to jld2

# plot(collect(eachrow(puls[10:11,:])),xlim = [-10,4000],ylim=[-0.04,0.04])
#plot(puls[7,:])

Tk = ones(500)*298.15
Iapp = ones(500)*0.1
CellV, jNeg, jPos, y, x, η0, ηL, ϕ_ẽ1, ϕ_ẽ2 = Sim_Model(CellData,Dtt,Iapp,Tk,A_DRA,B_DRA,C_DRA,D_DRA)