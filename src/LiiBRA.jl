module LiiBRA

using UnitSystems, Parameters, LinearAlgebra, FFTW
using TSVD, Roots, Statistics, Interpolations, JLD2
export C_e, Flux, C_se, Phi_s, Phi_e, Phi_se, CIDRA
export flatten_, R, F, Simulate, D_Linear, Construct
export tuple_len, interp, Spatial!, Realise, HPPC, fh!
export mag!, findnearest, CC, WLTP

include("Functions/C_e.jl")
include("Functions/C_se.jl")
include("Functions/Flux.jl")
include("Functions/Phi_s.jl")
include("Functions/Phi_e.jl")
include("Functions/Phi_se.jl")
include("Functions/Simulate.jl")
include("Methods/CIDRA.jl")

const F, R = faraday(Metric), universal(SI2019) #Faraday Constant / Universal Gas Constant 
findnearest(A, x) = argmin(abs.(A .- x)) # Find Nearest for SOC initialisation

#---------- Generate Model -----------------#
function Realise(Cell, Ŝ::Array)
    A = B = C = D = tuple()
    for i in Ŝ
        # Arrhenius
        Arr_Factor = (1 / Cell.Const.T_ref - 1 / Cell.Const.T) / R

        # Set Cell Constants
        Cell.Const.SOC = i
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0) * exp(Cell.Const.Ea_κ * Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs, Cell.RA.Nfft, Cell.RA.f)
        Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
        Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)

        # Realisation
        Aϕ, Bϕ, Cϕ, Dϕ = CIDRA(Cell)

        # Flatten output into Tuples
        A = flatten_(A, Aϕ)
        B = flatten_(B, Bϕ)
        C = flatten_(C, Cϕ)
        D = flatten_(D, Dϕ)
    end
    return A, B, C, D
end

#---------- Hankel Formation & SVD -----------------#
function fh!(H, len₁, len₂, puls, M, Pulsᴸ)
    """

    Inplace mutation of H, returns U,S,V'

    """
    @views for lp1 in 1:length(len₂), lp2 in 1:length(len₁)
        H[(Pulsᴸ * (lp2 - 1) + 1):(Pulsᴸ * lp2), lp1] .= puls[:, len₂[lp1] + len₁[lp2]]
    end

    Init = convert(Vector{float(eltype(H))}, ones(size(H, 1)))
    U, S, V = tsvd(H, M; initvec = Init)

    @views for lp1 in 1:length(len₂), lp2 in 1:length(len₁)
        H[(Pulsᴸ * (lp2 - 1) + 1):(Pulsᴸ * lp2), lp1] .= puls[:, len₂[lp1] + len₁[lp2] + 1]
    end

    for i in 1:size(U, 2)
        if U[1, i] < float(0)
            U[:, i] = -U[:, i]
            V[:, i] = -V[:, i]
        end
    end

    return U, S, V'
end

#---------- HPPC Simulation -----------------#
function HPPC(Cell, Ŝ::Array, SOC::Float64, λ::Float64, ϕ::Float64, A::Tuple, B::Tuple,
              C::Tuple, D::Tuple)

    # Set Experiment
    i = Int64(1 / Cell.RA.SamplingT) #Sampling Frequency
    Input = [zero(i); ones(10 * i) * λ; zeros(40 * i); ones(10 * i) * ϕ; zeros(40 * i)]
    Tk = ones(size(Input)) * Cell.Const.T #Cell Temperature
    t = 0:(1.0 / i):((length(Input) - 1) / i)

    # Simulate Model
    return Simulate(Cell, Input, "Current", Tk, Ŝ, SOC, A, B, C, D, t)
end

#---------- Constant Current Simulation -----------------#
function CC(Cell, Ŝ::Array, SOC::Float64, λ::Float64, γ, A::Tuple, B::Tuple, C::Tuple,
            D::Tuple)

    # Set Experiment
    i = 1 / Cell.RA.SamplingT
    Input = ones(Int64(γ * i)) * λ #CC Profile
    Tk = ones(size(Input)) * Cell.Const.T
    t = 0:(1.0 / i):((length(Input) - 1) / i)

    # Simulate Model
    return Simulate(Cell, Input, "Current", Tk, Ŝ, SOC, A, B, C, D, t)
end

#---------- WLTP Simulation -----------------#
function WLTP(Cell, Ŝ::Array, SOC::Float64, Cycle::Array, A::Tuple, B::Tuple, C::Tuple,
              D::Tuple)

    # Set Experiment
    i = Int64(1 / Cell.RA.SamplingT)
    Tk = ones(size(Cycle)) * Cell.Const.T
    t = 0:(1.0 / i):((length(Cycle) - 1) / i)

    # Simulate Model
    return Simulate(Cell, Cycle, "Power", Tk, Ŝ, SOC, A, B, C, D, t)
end

#---------- Tuple Flatten -----------------#
"""

    flatten_(a::Tuple, b...) 

    Flattens input Tuple "a" and inserts "b" 

"""
function flatten_ end
flatten_() = ()
flatten_(a::Tuple) = Tuple(a)
flatten_(a) = (a,)
flatten_(a::Tuple, b...) = tuple(a..., flatten_(b...)...)
flatten_(a, b...) = tuple(a, flatten_(b...)...)
flatten_tuple(x::Tuple) = flatten_(x...)

#---------- Cell Contstruct -----------------#
"""
    Construct(::String) 

    Function to create Construct dictionary.

    Currently supports:

    1. Doyle '94 parameterisation
    2. Chen 2020 parameterisation

"""
function Construct(CellType::String)
    if CellType == "Doyle_94"
        CellType = string(CellType, ".jl")
        return include(joinpath(dirname(pathof(LiiBRA)), "Data/Doyle_94", CellType))
    elseif CellType == "LG M50"
        CellType = string("LG_M50.jl")
        return include(joinpath(dirname(pathof(LiiBRA)), "Data/Chen_2020", CellType))
    elseif CellType == "A123"
        CellType = string("A123.jl")
        return include(joinpath(dirname(pathof(LiiBRA)), "Data/Prada_2013", CellType))
    end
end

#---------- Cell Update -----------------#
"""
    Spatial!(::String) 

    Function to update the cell dictionary

"""
function Spatial!(Cell, Sₑ, Sₛ)
    Cell.Transfer.Sₑ = Sₑ
    Cell.Transfer.Sₛ = Sₛ
    Cell.Const.Lnegsep, Cell.Const.Ltot = Cell.Neg.L + Cell.Sep.L,
                                          Cell.Neg.L + Cell.Sep.L + Cell.Pos.L
    Cell.Const.D1 = Cell.Const.De * Cell.Neg.ϵ_e^Cell.Neg.De_brug
    Cell.Const.D2 = Cell.Const.De * Cell.Sep.ϵ_e^Cell.Sep.De_brug
    Cell.Const.D3 = Cell.Const.De * Cell.Pos.ϵ_e^Cell.Pos.De_brug
    Cell.Const.Ce_M = size(Cell.Transfer.Locs(Sₑ, Sₛ)[1], 1)
    Cell.RA.Outs = sum([size(Cell.Transfer.Locs(Sₑ, Sₛ)[i],
                             1)
                        for i in 1:length(Cell.Transfer.tfs)])
end

#---------- Tuple Length -----------------#
"""
    tuple_len(::NTuple) 

    Function to return Tuple length. 

"""
tuple_len(::NTuple{N, Any}) where {N} = N #Tuple Size

#---------- Interpolate -----------------#
"""
    interp(MTup::Tuple,Ŝ::Array,SOC) 

    Function to interpolate Tuple indices

"""
function interp(MTup::Tuple, Ŝ::Array, SOC)
    T1 = 1
    T2 = 1
    for i in 1:(length(Ŝ) - 1)
        if Ŝ[i] > SOC >= Ŝ[i + 1]
            T1 = i
            T2 = i + 1
        elseif Ŝ[i] == SOC
            return MTup[i]
        end
    end
    return M = @. MTup[T2] + (MTup[T1] - MTup[T2]) * (SOC - Ŝ[T2]) / (Ŝ[T1] - Ŝ[T2])
end

#---------- Magnitude of an Array -----------------#
function mag(γ::AbstractArray)
    ψ = Array{Float64}(undef, size(γ))
    for j in 1:size(γ, 2), i in 1:size(γ, 1)
        if real(γ[i, j]) < 0
            ψ[i, j] = -abs(γ[i, j])
        else
            ψ[i, j] = abs(γ[i, j])
        end
    end
    return ψ
end

#---------- Magnitude of a Vector -----------------#
function mag(γ::AbstractVector)
    ψ = Vector{Float64}(undef, length(γ))
    for i in 1:length(γ)
        if real(γ[i]) < 0
            ψ[i] = -abs(γ[i])
        else
            ψ[i] = abs(γ[i])
        end
    end
    return ψ
end

#---------- D Lineaisation -----------------#
"""
    D_Linear(Cell,ν_neg,ν_pos,σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep) 

    Function to linearise D array from input cell type and conductivites. 

"""
function D_Linear(Cell, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep)
    # Spatial Domain
    Sₑ = Cell.Transfer.Sₑ
    Sₛ = Cell.Transfer.Sₛ

    #Initial
    D = Array{Float64}(undef, 0, 1)
    Dt = Array{Float64}(undef, 0, 1)
    D_ = Array{Float64}(undef, size(Cell.Transfer.Locs(Sₑ, Sₛ)[2], 1), 1)
    q = Int64(1)
    tfs = Cell.Transfer.tfs
    Elec = Cell.Transfer.Elec
    Locs = Cell.Transfer.Locs(Sₑ, Sₛ)
    for i in 1:size(tfs, 1)
        if tfs[i] == C_e
            Dt = zeros(length(Locs[i]))
        elseif tfs[i] == Phi_e
            Dt = Array{Float64}(undef, 0, 1)
            for pt in Locs[i]
                if pt <= Cell.Neg.L + eps()
                    D_[q] = @. (Cell.Neg.L * (σ_eff_Neg / κ_eff_Neg) *
                                (1 - cosh(ν_neg * pt / Cell.Neg.L)) -
                                pt * ν_neg * sinh(ν_neg) +
                                Cell.Neg.L *
                                (cosh(ν_neg) - cosh(ν_neg * (Cell.Neg.L - pt) / Cell.Neg.L))) /
                               (Cell.Const.CC_A * (κ_eff_Neg + σ_eff_Neg) * sinh(ν_neg) *
                                ν_neg)
                elseif pt <= Cell.Neg.L + Cell.Sep.L + eps()
                    D_[q] = @. (Cell.Neg.L - pt) / (Cell.Const.CC_A * κ_eff_Sep) +
                               (Cell.Neg.L *
                                ((1 - σ_eff_Neg / κ_eff_Neg) * tanh(ν_neg / 2) - ν_neg)) /
                               (Cell.Const.CC_A * (κ_eff_Neg + σ_eff_Neg) * ν_neg)
                else
                    D_[q] = @. -Cell.Sep.L / (Cell.Const.CC_A * κ_eff_Sep) +
                               Cell.Neg.L *
                               ((1 - σ_eff_Neg / κ_eff_Neg) * tanh(ν_neg / 2) - ν_neg) /
                               (Cell.Const.CC_A * (σ_eff_Neg + κ_eff_Neg) * ν_neg) +
                               (Cell.Pos.L * (-σ_eff_Pos * cosh(ν_pos) +
                                 σ_eff_Pos *
                                 cosh((Cell.Const.Ltot - pt) * ν_pos / Cell.Pos.L) +
                                 κ_eff_Pos *
                                 (cosh((pt - Cell.Neg.L - Cell.Sep.L) * ν_pos / Cell.Pos.L) -
                                  1)) -
                                (pt - Cell.Neg.L - Cell.Sep.L) * κ_eff_Pos * sinh(ν_pos) *
                                ν_pos) /
                               (Cell.Const.CC_A * κ_eff_Pos * (κ_eff_Pos + σ_eff_Pos) *
                                ν_pos * sinh(ν_pos))
                end
                Dt = [Dt; D_[q]]
                q = q + 1
            end

        elseif tfs[i] == C_se
            Dt = zeros(length(Locs[i]))
        elseif Elec[i] == "Pos"
            σ_eff = σ_eff_Pos
            κ_eff = κ_eff_Pos
            if tfs[i] == Phi_se
                Dt = @. -1 * Cell.Pos.L / (Cell.Const.CC_A * ν_pos * sinh(ν_pos)) *
                        ((1 / κ_eff) * cosh(ν_pos * Locs[i]) +
                         (1 / σ_eff) * cosh(ν_pos * (Locs[i] - 1)))
            elseif tfs[i] == Flux
                Dt = @. -1 * ν_pos *
                        (σ_eff * cosh(ν_pos * Locs[i]) +
                         κ_eff * cosh(ν_pos * (Locs[i] - 1))) /
                        (Cell.Pos.as * F * Cell.Pos.L * Cell.Const.CC_A * (κ_eff + σ_eff) *
                         sinh(ν_pos))
            elseif tfs[i] == Phi_s
                Dt = @. -1 *
                        (-Cell.Pos.L * (κ_eff * (cosh(ν_pos) - cosh(Locs[i] - 1) * ν_pos)) -
                         Cell.Pos.L * (σ_eff *
                          (1 - cosh(Locs[i] * ν_pos) + Locs[i] * ν_pos * sinh(ν_pos)))) /
                        (Cell.Const.CC_A * σ_eff * (κ_eff + σ_eff) * ν_pos * sinh(ν_pos))
            end
        elseif Elec[i] == "Neg"
            σ_eff = σ_eff_Neg
            κ_eff = κ_eff_Neg
            if tfs[i] == Phi_se
                Dt = @. Cell.Neg.L / (Cell.Const.CC_A * ν_neg * sinh(ν_neg)) *
                        ((1 / κ_eff) * cosh(ν_neg * Locs[i]) +
                         (1 / σ_eff) * cosh(ν_neg * (Locs[i] - 1)))
            elseif tfs[i] == Flux
                Dt = @. ν_neg * (σ_eff * cosh(ν_neg * Locs[i]) +
                         κ_eff * cosh(ν_neg * (Locs[i] - 1))) /
                        (Cell.Neg.as * F * Cell.Neg.L * Cell.Const.CC_A * (κ_eff + σ_eff) *
                         sinh(ν_neg))
            elseif tfs[i] == Phi_s
                Dt = @. (-Cell.Neg.L * (κ_eff * (cosh(ν_neg) - cosh(Locs[i] - 1) * ν_neg)) -
                         Cell.Neg.L * (σ_eff *
                          (1 - cosh(Locs[i] * ν_neg) + Locs[i] * ν_neg * sinh(ν_neg)))) /
                        (Cell.Const.CC_A * σ_eff * (κ_eff + σ_eff) * ν_neg * sinh(ν_neg))
            end
        end
        D = [D; Dt]
    end
    return D
end

end # module
