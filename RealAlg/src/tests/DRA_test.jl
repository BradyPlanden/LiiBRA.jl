using RealAlg, UnitSystems

CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator(),RealisationAlgorthim())
loc = (Number[0, 128e-6, 204e-6, 394e-6],Number[0,1], Number[0,1], Number[128e-6, 204e-6, 394e-6],Number[1],Number[1],Number[0,1],Number[0,1],Number[0,1],Number[0,1])


Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/universal(SI2019)
FCall= FCalls(kappa(CellData.Const.ce0))
κ_dra=FCall.Kap.κ*exp(CellData.Const.Ea_κ*Arr_Factor)


DRA(CellData,FCall,loc)