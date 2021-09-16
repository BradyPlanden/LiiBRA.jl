module LIBRA

using Roots, UnitSystems, Parameters, LinearAlgebra, FFTW, Statistics
using Dierckx, Arpack, Printf, JLD, PROPACK, ControlSystems
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
include("Data/Chen_2020/LG_M50.jl")
#include("Data/Doyle_94/Doyle_94.jl")

const ϵ1,ϵ2,ϵ3 = CellData.Neg.ϵ_e, CellData.Sep.ϵ_e, CellData.Pos.ϵ_e        # Porosities
const D1,D2,D3  = CellData.Const.De*ϵ1^CellData.Neg.De_brug, CellData.Const.De*ϵ2^CellData.Sep.De_brug, CellData.Const.De*ϵ3^CellData.Pos.De_brug   #Effective diffusivities
const Lneg, Lpos, Lnegsep, Ltot = CellData.Neg.L, CellData.Pos.L, CellData.Neg.L+CellData.Sep.L,CellData.Neg.L+CellData.Sep.L+CellData.Pos.L

const F,R = faraday(Metric), universal(SI2019)     # Faraday Constant / Universal Gas Constant
const Debug = 0 #Print Variables for Debugging    

function D_Linear(tfs, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos)
    D = Array{Float64}(undef,0,1)
    Dt = Array{Float64}(undef,0,1)
    for i in 1:size(tfs,1)
        if tfs[i,1] == C_e
            Dt =  zeros(length(tfs[i,3]))
        elseif tfs[i,1] == Phi_e
            Dt = zeros(length(tfs[i,3]))
        elseif tfs[i,1] == C_se
            Dt = zeros(length(tfs[i,3]))
        elseif tfs[i,2] == "Pos"
            σ_eff = σ_eff_Pos
            κ_eff = κ_eff_Pos
            if tfs[i,1] == Phi_se
                Dt = @. CellData.Pos.L/(CellData.Const.CC_A*ν_pos*sinh(ν_pos))*((1/κ_eff)*cosh(ν_pos*tfs[i,3])+(1/σ_eff)*cosh(ν_pos*(tfs[i,3]-1))) # Contribution to D as G->∞
            elseif tfs[i,1] == Flux
                Dt = @. ν_pos*(σ_eff*cosh(ν_pos*tfs[i,3])+κ_eff*cosh(ν_pos*(tfs[i,3]-1)))/(CellData.Pos.as*F*CellData.Pos.L*CellData.Const.CC_A*(κ_eff+σ_eff)*sinh(ν_pos))
            elseif tfs[i,1] == Phi_s
                Dt = @. (-CellData.Pos.L*(κ_eff*(cosh(ν_pos)-cosh(tfs[i,3]-1)*ν_pos))-CellData.Pos.L*(σ_eff*(1-cosh(tfs[i,3]*ν_pos)+tfs[i,3]*ν_pos*sinh(ν_pos))))/(CellData.Const.CC_A*σ_eff*(κ_eff+σ_eff)*ν_pos*sinh(ν_pos)) # Contribution to D as G->∞ 
            end
        elseif tfs[i,2] == "Neg"
            σ_eff = σ_eff_Neg
            κ_eff = κ_eff_Neg
            if tfs[i,1] == Phi_se
                Dt = @. CellData.Neg.L/(CellData.Const.CC_A*ν_neg*sinh(ν_neg))*((1/κ_eff)*cosh(ν_neg*tfs[i,3])+(1/σ_eff)*cosh(ν_neg*(tfs[i,3]-1))) # Contribution to D as G->∞
            elseif tfs[i,1] == Flux
                Dt = @. ν_neg*(σ_eff*cosh(ν_neg*tfs[i,3])+κ_eff*cosh(ν_neg*(tfs[i,3]-1)))/(CellData.Neg.as*F*CellData.Neg.L*CellData.Const.CC_A*(κ_eff+σ_eff)*sinh(ν_neg))
            elseif tfs[i,1] == Phi_s
                Dt = @. (-CellData.Neg.L*(κ_eff*(cosh(ν_neg)-cosh(tfs[i,3]-1)*ν_neg))-CellData.Neg.L*(σ_eff*(1-cosh(tfs[i,3]*ν_neg)+tfs[i,3]*ν_neg*sinh(ν_neg))))/(CellData.Const.CC_A*σ_eff*(κ_eff+σ_eff)*ν_neg*sinh(ν_neg)) # Contribution to D as G->∞ 
            end
        end
        D = [D; Dt]
        #Iterate through D terms and export new D Matrix
        #Inputs: All parameters needed? or some method to capture the numerical values of those parameters and dynamically write them
        #Outputs: Vector containing the updated D terms
        #Logic to export only the requested D term or all terms (multiple dispatch?)
    end
    return D
end

function flatten end
flatten() = ()
flatten(a::Tuple) = Tuple(a)
flatten(a) = (a,)
flatten(a::Tuple, b...) = tuple(a..., flatten(b...)...)
flatten(a, b...) = tuple(a, flatten(b...)...)
flatten_tuple(x::Tuple) = flatten(x...)

end # module