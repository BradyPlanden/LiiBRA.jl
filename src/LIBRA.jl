module LIBRA

using Roots, UnitSystems, Parameters, LinearAlgebra, FFTW, Dierckx, Arpack, Printf, JLD, TSVD #, PROPACK

import Base: +,-,*,==,>,>=,<,<=,broadcast,sin,cos,tan,cot,abs,exp,log,log10
export CellData, C_e, Negative, Constants, Positive, Seperator, Flux, C_se, Phi_s, Phi_e, Phi_se, DRA, RealisationAlgorthim, TransferFun, flatten, R, F, CellTyp, Sim_Model

include("LIBRATypes.jl")
include("Functions/Transfer/C_e.jl")
include("Functions/Transfer/C_se.jl")
include("Functions/Transfer/Flux.jl")
include("Functions/Transfer/Phi_s.jl")
include("Functions/Transfer/Phi_e.jl")
include("Functions/Transfer/Phi_se.jl")
include("Methods/DRA.jl")
include("Functions/Sim_Model.jl")
#include("Data/Chen_2020/LG_M50.jl")
include("Data/Doyle_94/Doyle_94.jl")

const CellTyp = "LG_M50"
const ϵ1,ϵ2,ϵ3 = CellData.Neg.ϵ_e, CellData.Sep.ϵ_e, CellData.Pos.ϵ_e        # Porosities
const D1,D2,D3  = CellData.Const.De*ϵ1^CellData.Neg.De_brug, CellData.Const.De*ϵ2^CellData.Sep.De_brug, CellData.Const.De*ϵ3^CellData.Pos.De_brug   #Effective diffusivities
const Lneg, Lpos, Lnegsep, Ltot = CellData.Neg.L, CellData.Pos.L, CellData.Neg.L+CellData.Sep.L,CellData.Neg.L+CellData.Sep.L+CellData.Pos.L

const F,R = faraday(Metric), universal(SI2019)     # Faraday Constant / Universal Gas Constant
const Debug = 0 #Print Variables for Debugging    

function flatten end
flatten() = ()
flatten(a::Tuple) = Tuple(a)
flatten(a) = (a,)
flatten(a::Tuple, b...) = tuple(a..., flatten(b...)...)
flatten(a, b...) = tuple(a, flatten(b...)...)
flatten_tuple(x::Tuple) = flatten(x...)

end # module