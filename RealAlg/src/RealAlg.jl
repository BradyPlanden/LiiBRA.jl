__precompile__()

module RealAlg

using Roots, UnitSystems, DataFrames, CSV, Parameters, LinearAlgebra, FFTW, DataInterpolations

import Base: +,-,*,==,>,>=,<,<=,broadcast,sin,cos,tan,cot,abs,exp,log,log10
export Cell, C_e, Negative, Constants, Geometry, Positive, Seperator, j, ∂Uocp, C_se, Phi_s, Phi_e, Phi_se, DRA, RealisationAlgorthim, kappa, FCalls

include("RealAlgTypes.jl")
include("Functions/C_e.jl")
include("Functions/C_se.jl")
include("Functions/Flux.jl")
include("Functions/Phi_s.jl")
include("Functions/Phi_e.jl")
include("Functions/Phi_se.jl")
include("Methods/DRA.jl")

CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator(),RealisationAlgorthim())
const Lpos = CellData.Pos.L
const Lneg = CellData.Neg.L
const Lsep = CellData.Sep.L
const Ltot = CellData.Geo.Ltot
const Lnegsep = Lneg+Lsep

const ϵ1 = CellData.Neg.ϵ_e      # Porosity of negative electrode
const ϵ2 = CellData.Sep.ϵ_e      # Porosity of separator
const ϵ3 = CellData.Pos.ϵ_e      # Porosity of positive electrode
const D1 = CellData.Const.De * ϵ1^CellData.Neg.De_brug # Effective ...
const D2 = CellData.Const.De * ϵ2^CellData.Sep.De_brug # diffusivities ...
const D3 = CellData.Const.De * ϵ3^CellData.Pos.De_brug # of cell regions
const F = faraday(Metric)      # Faraday Constant
const R = universal(SI2019)       # Universal Gas Constant
const as_neg = 3*CellData.Neg.ϵ_s/CellData.Neg.Rs # Specific interfacial surf. area
const as_pos = 3*CellData.Pos.ϵ_s/CellData.Pos.Rs # Specific interfacial surf. area

const Debug = 1 #Print Variables for Debugging    

@inline function ∂Uocp(Electrode,θ)
    if Electrode == "Neg"
    ∂Uocp = (-20000*exp(-2000*θ) - 3.96*exp(-3*θ))
    else
    ∂Uocp = (-32.4096*exp(-40*(-0.133875 + θ)) - 0.0135664./((0.998432 - θ).^1.49247)+ 0.0595559*exp(-0.04738*θ.^8).*θ.^7 - 0.823297*(sech(8.60942 - 14.5546*θ)).^2)
    end
end

# @inline function κ_funct(ce)
#     return @. (4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24)
# end

end # module