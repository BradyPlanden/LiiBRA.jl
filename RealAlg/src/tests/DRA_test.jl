using RealAlg

CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator(),RealisationAlgorthim())
loc = (Number[0, 128e-6, 204e-6, 394e-6],Number[0,1], Number[0,1], Number[128e-6, 204e-6, 394e-6],Number[1],Number[1],Number[0,1],Number[0,1],Number[0,1],Number[0,1])

T=298.15
Arr_Factor = (1/CellData.Const.T-1/T)/R
CellData.const.κ_ra = κ(CellData.const.ce0)*exp(cell.const.Ea_κ*Arr_Factor)

DRA(CellData,loc)