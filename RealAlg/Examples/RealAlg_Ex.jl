using RealAlg, JLD

function CData() # Create structs
    #CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator(),RealisationAlgorthim())
    #loc = (Number[0, 128e-6, 204e-6, 394e-6],Number[128e-6, 204e-6, 394e-6], Number[0,1], Number[128e-6, 204e-6, 394e-6],Number[1],Number[1],Number[0,1],Number[0,1],Number[0,1],Number[0,1])
    TransferFuns = TransferFun()
    return TransferFuns
end


@inline function Impulse() # Create s Vector 
    Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    f = 0:Nfft-1
    s = ((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft))'
    return Nfft,f,s
end

@inline function DRA_loop(CellData,s,f,TransferFuns)
    A_DRA = B_DRA = C_DRA = D_DRA = Dtt = tuple()

    for Temp in 5.0:5.0:5.0
        CellData.Const.T = 273.15+Temp
        Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/R
        CellData.Const.κ = CellData.Const.κf(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
            for SOC in 0:0.05:0.0
                CellData.Const.SOC = SOC
                A, B, C, D, Dtt = DRA(CellData,s,f,TransferFuns)
                A_DRA = flatten(A_DRA,A)
                B_DRA = flatten(B_DRA,B)
                C_DRA = flatten(C_DRA,C)
                D_DRA = flatten(D_DRA,D)
            end
    end

return A_DRA, B_DRA, C_DRA, D_DRA, Dtt
end


TransferFuns = CData()
Nfft, f, s = Impulse()
A_DRA, B_DRA, C_DRA, D_DRA, Dtt = DRA_loop(CellData,s,f,TransferFuns)
save("RealAlg_$CellTyp.jld", "CellData", CellData, "A_DRA", A_DRA, "B_DRA", B_DRA, "C_DRA", C_DRA, "D_DRA", D_DRA, "Dtt", Dtt)
