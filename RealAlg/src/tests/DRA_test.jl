using RealAlg, UnitSystems

CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator(),RealisationAlgorthim())
loc = (Number[0, 128e-6, 204e-6, 394e-6],Number[0,1], Number[0,1], Number[128e-6, 204e-6, 394e-6],Number[1],Number[1],Number[0,1],Number[0,1],Number[0,1],Number[0,1])


Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/universal(SI2019)
FCall= FCalls(kappa(CellData.Const.ce0))
κ_dra=FCall.Kap.κ*exp(CellData.Const.Ea_κ*Arr_Factor)
TransferFuns = TransferFun()


const A_DRA = Array{Float64}(undef,0,CellData.RA.M)
const B_DRA = Array{Float64}(undef,0,CellData.RA.M)
const C_DRA = Array{Float64}(undef,0,CellData.RA.M)
const D_DRA = Array{Float64}(undef,0,CellData.RA.M)
for Temperature in 5.0:5.0:5.0
    CellData.Const.T = 273.15+Temperature
    for SOC in 0:0.05:0.05
        CellData.Const.SOC = SOC
        A, B, C, D = DRA(CellData,FCall,loc,TransferFuns)
        println("DRA")
    end
end