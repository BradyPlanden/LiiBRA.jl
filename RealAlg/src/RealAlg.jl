__precompile__()

module RealAlg

using Roots, UnitSystems, DataFrames, CSV, Parameters

import Base: +,-,*,==,>,>=,<,<=,broadcast,sin,cos,tan,cot,abs,exp,log,log10
export Cell, C_e, Negative, Constants, Geometry, Positive, Seperator 

include("RealAlgTypes.jl")
include("Functions/C_e.jl")

CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator())
const Lpos = CellData.Geo.Lpos
const Lneg = CellData.Geo.Lneg
const Lsep = CellData.Geo.Lsep
const Ltot = CellData.Geo.Ltot
const Lnegsep = Lneg+Lsep

const ϵ1 = CellData.Neg.ϵ_e      # Porosity of negative electrode
const ϵ2 = CellData.Sep.ϵ_e      # Porosity of separator
const ϵ3 = CellData.Pos.ϵ_e      # Porosity of positive electrode
const D1 = CellData.Const.De * ϵ1^CellData.Neg.De_brug # Effective ...
const D2 = CellData.Const.De * ϵ2^CellData.Sep.De_brug # diffusivities ...
const D3 = CellData.Const.De * ϵ3^CellData.Pos.De_brug # of cell regions

end # module