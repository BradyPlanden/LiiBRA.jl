using LIBRA, JLD

#Cell Definition
CellTyp = "LG_M50_1"
#include("Data/Chen_2020/LG_M50.jl")

@inline function Impulse() # Create s Vector 
    Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    f = 0:Nfft-1
    s = ((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft))'
    return Nfft,f,s
end

@inline function DRA_loop(CellData,s,f)
    A_DRA = B_DRA = C_DRA = D_DRA = Dtt = tuple()
    for Temp in 5.0:5.0:5.0
        CellData.Const.T = 273.15+Temp
        Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/R
        CellData.Const.κ = CellData.Const.κf(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
            for SOC in 0:0.05:0.0
                CellData.Const.SOC = SOC
                A, B, C, D, Dtt = DRA(CellData,s,f)
                A_DRA = flatten(A_DRA,A)
                B_DRA = flatten(B_DRA,B)
                C_DRA = flatten(C_DRA,C)
                D_DRA = flatten(D_DRA,D)
            end
    end

return A_DRA, B_DRA, C_DRA, D_DRA, Dtt
end

#TransferFuns = TransferFun()
Nfft, f, s = Impulse()
A_DRA, B_DRA, C_DRA, D_DRA, Dtt = DRA_loop(CellData,s,f)
#save("$CellTyp.jld", "CellData", CellData, "A_DRA", A_DRA, "B_DRA", B_DRA, "C_DRA", C_DRA, "D_DRA", D_DRA, "Dtt", Dtt) #Switch to jld2

Tk = ones(100)*25
Iapp = ones(100)
CellV = Sim_Model(CellData,Dtt,Iapp,Tk,A_DRA,B_DRA,C_DRA,D_DRA)
