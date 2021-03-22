using RealAlg, UnitSystems

CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator(),RealisationAlgorthim())
loc = (Number[0, 128e-6, 204e-6, 394e-6],Number[0,1], Number[0,1], Number[128e-6, 204e-6, 394e-6],Number[1],Number[1],Number[0,1],Number[0,1],Number[0,1],Number[0,1])

TransferFuns = TransferFun()

# Create s Vector 

function Impulse()
    Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
    f = 0:Nfft-1
    s = ((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft))'
    return Nfft,f,s
end

function DRA_loop(CellData,s,f,loc,TransferFuns)
    A_DRA = Array{Float64}(undef,0,CellData.RA.M)
    B_DRA = Array{Float64}(undef,0,CellData.RA.M)
    C_DRA = Array{Float64}(undef,0,CellData.RA.M)
    D_DRA = Array{Float64}(undef,0,CellData.RA.M)

    for Temp in 5.0:5.0:5.0
        CellData.Const.T = 273.15+Temp
        Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/universal(SI2019)
        CellData.Const.κ = Kappa(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
            for SOC in 0:0.05:0.0
                CellData.Const.SOC = SOC
                A, B, C, D = DRA(CellData,s,f,loc,TransferFuns)
            end
    end
end
Nfft, f, s = Impulse()
DRA_loop(CellData,s,f,loc,TransferFuns)