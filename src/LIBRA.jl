module LIBRA

using UnitSystems, Parameters, LinearAlgebra, FFTW
using Dierckx, Arpack, Infiltrator, Statistics
export C_e, Negative, Constants, Positive, Seperator, Flux, C_se, Phi_s, Phi_e, Phi_se, DRA
export RealisationAlgorthim, TransferFun, flatten, R, F, CellDef, Sim_Model, D_Linear, _bisection, cell

include("Functions/Transfer/C_e.jl")
include("Functions/Transfer/C_se.jl")
include("Functions/Transfer/Flux.jl")
include("Functions/Transfer/Phi_s.jl")
include("Functions/Transfer/Phi_e.jl")
include("Functions/Transfer/Phi_se.jl")
include("Methods/DRA.jl")
include("Functions/Sim_Model.jl")

const F,R = faraday(Metric), universal(SI2019) #Faraday Constant / Universal Gas Constant
const Debug = 0 #Print Variables for Debugging    

function D_Linear(CellData, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep)
    D = Array{Float64}(undef,0,1)
    Dt = Array{Float64}(undef,0,1)
    D_ = Array{Float64}(undef,5,1)
    q = Int64(1)
    tfs = CellData.Transfer.tfs
    for i in 1:size(tfs,1)
        if tfs[i,1] == C_e
            Dt =  zeros(length(tfs[i,3]))
        elseif tfs[i,1] == Phi_e
            Dt = Array{Float64}(undef,0,1)
            for pt in tfs[i,3]
                if pt <= CellData.Neg.L+eps()
                D_[q]  =  @. (CellData.Neg.L*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg*pt/CellData.Neg.L)) - pt*ν_neg*sinh(ν_neg) + CellData.Neg.L*(cosh(ν_neg)-cosh(ν_neg*(CellData.Neg.L-pt)/CellData.Neg.L)))/(CellData.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg) #Lee. Eqn. 4.22 @ ∞
                elseif pt <= CellData.Neg.L + CellData.Sep.L + eps()
                D_[q] = @. (CellData.Neg.L - pt)/(CellData.Const.CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/(CellData.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) #Lee. Eqn. 4.23 @ ∞
                else
                D_[q] = @. -CellData.Sep.L/(CellData.Const.CC_A*κ_eff_Sep) + CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg)/(CellData.Const.CC_A*(σ_eff_Neg+ κ_eff_Neg)*ν_neg) + (CellData.Pos.L*(-σ_eff_Pos*cosh(ν_pos) + σ_eff_Pos*cosh((CellData.Const.Ltot-pt)*ν_pos/CellData.Pos.L) +  κ_eff_Pos*(cosh((pt-CellData.Neg.L-CellData.Sep.L)*ν_pos/CellData.Pos.L)-1)) - (pt-CellData.Neg.L-CellData.Sep.L)*κ_eff_Pos*sinh(ν_pos)*ν_pos)/(CellData.Const.CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_pos*sinh(ν_pos))
                end
                Dt = [Dt; D_[q]]
                q = q+1
            end

        elseif tfs[i,1] == C_se
            Dt = zeros(length(tfs[i,3]))
        elseif tfs[i,2] == "Pos"
            σ_eff = σ_eff_Pos
            κ_eff = κ_eff_Pos
            if tfs[i,1] == Phi_se
                Dt = @. -1*CellData.Pos.L/(CellData.Const.CC_A*ν_pos*sinh(ν_pos))*((1/κ_eff)*cosh(ν_pos*tfs[i,3])+(1/σ_eff)*cosh(ν_pos*(tfs[i,3]-1))) # Contribution to D as G->∞
            elseif tfs[i,1] == Flux
                Dt = @. -1*ν_pos*(σ_eff*cosh(ν_pos*tfs[i,3])+κ_eff*cosh(ν_pos*(tfs[i,3]-1)))/(CellData.Pos.as*F*CellData.Pos.L*CellData.Const.CC_A*(κ_eff+σ_eff)*sinh(ν_pos))
            elseif tfs[i,1] == Phi_s
                Dt = @. -1*(-CellData.Pos.L*(κ_eff*(cosh(ν_pos)-cosh(tfs[i,3]-1)*ν_pos))-CellData.Pos.L*(σ_eff*(1-cosh(tfs[i,3]*ν_pos)+tfs[i,3]*ν_pos*sinh(ν_pos))))/(CellData.Const.CC_A*σ_eff*(κ_eff+σ_eff)*ν_pos*sinh(ν_pos)) # Contribution to D as G->∞ 
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


"""
    bisection(f, a, b; fa = f(a), fb = f(b), ftol, wtol)

Bisection algorithm for finding the root ``f(x) ≈ 0`` within the initial bracket
`[a,b]`.

Returns a named tuple

`(x = x, fx = f(x), isroot = ::Bool, iter = ::Int, ismaxiter = ::Bool)`.

Terminates when either

1. `abs(f(x)) < ftol` (`isroot = true`),
2. the width of the bracket is `≤wtol` (`isroot = false`),
3. `maxiter` number of iterations is reached. (`isroot = false, maxiter = true`).

which are tested for in the above order. Therefore, care should be taken not to make `wtol` too large.

"""
function bisection(f,CellData, a::Real, b::Real; fa::Real = f(a), fb::Real = f(b),
                   ftol = √eps(), wtol = 0, maxiter = 100)
    @assert fa * fb ≤ 0 "initial values don't bracket zero"
    @assert isfinite(a) && isfinite(b)
    _bisection(f, CellData, float.(promote(a, b, fa, fb, ftol, wtol))..., maxiter)
end

function _bisection(f,CellData, a, b, fa, fb, ftol, wtol, maxiter)
    iter = 0
    abs(fa) < ftol && return (x = a, fx = fa, isroot = true, iter = iter, ismaxiter = false)
    abs(fb) < ftol && return (x = b, fx = fb, isroot = true, iter = iter, ismaxiter = false)
    while true
        iter += 1
        m = middle(a, b)
        fm = f(CellData,m)
        abs(fm) < ftol && return (x = m, fx = fm, isroot = true, iter = iter, ismaxiter = false)
        abs(b-a) ≤ wtol && return (x = m, fx = fm, isroot = false, iter = iter, ismaxiter = false)
        if fa * fm > 0
            a, fa = m, fm
        else
            b, fb = m, fm
        end
        iter == maxiter && return (x = m, fx = fm, isroot = false, iter = iter, ismaxiter = true)
    end
end

function cell(CellType)
    if CellType == "Doyle_94"
        CellType = string(CellType,".jl")
        include(joinpath(dirname(pathof(LIBRA)), "Data/Doyle_94", CellType))
    elseif CellType == "LG_M50"
        CellType = string(CellType,".jl")
        include(joinpath(dirname(pathof(LIBRA)), "Data/Chen_2020", CellType))
    end
    return CellData
end

end # module