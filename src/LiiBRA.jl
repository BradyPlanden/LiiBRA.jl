module LiiBRA

using UnitSystems, Parameters, LinearAlgebra, FFTW
using TSVD, Roots, Statistics, Interpolations
export C_e, Flux, C_se, Phi_s, Phi_e, Phi_se, DRA
export flatten_, R, F, Sim_Model, D_Linear, Construct, tuple_len, interp
export Realise, HPPC, fh!, mag!

include("Functions/Transfer/C_e.jl")
include("Functions/Transfer/C_se.jl")
include("Functions/Transfer/Flux.jl")
include("Functions/Transfer/Phi_s.jl")
include("Functions/Transfer/Phi_e.jl")
include("Functions/Transfer/Phi_se.jl")
include("Methods/DRA.jl")
include("Functions/Sim_Model.jl")

const F,R = faraday(Metric), universal(SI2019) #Faraday Constant / Universal Gas Constant 


#---------- Generate Model -----------------#
function Realise(Cell, SList::Array, T::Float64)
    A = B = C = D = tuple()
    for i in SList
        #Arrhenius
        Cell.Const.T = T
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R

        #Set Cell Constants
        Cell.Const.SOC = i
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)
        Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
        Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)

        #Realisation
        Aϕ, Bϕ, Cϕ, Dϕ = DRA(Cell)

        #Flatten output into Tuples
        A = flatten_(A,Aϕ)
        B = flatten_(B,Bϕ)
        C = flatten_(C,Cϕ)
        D = flatten_(D,Dϕ)
    end
return A, B, C, D
end



#---------- Simulate Model -----------------#
function HPPC(Cell,SList::Array,SOC::Float64,λ::Float64,ϕ::Float64,A::Tuple,B::Tuple,C::Tuple,D::Tuple)

    #Set Experiment
    i = Int64(1/Cell.RA.SamplingT) #Sampling Frequency
    Iapp = [ones(1)*0.; ones(10*i)*λ; ones(40*i)*0.; ones(10*i)*ϕ; ones(40*i)*0.] #1C HPPC Experiment Current Profile
    Tk = ones(size(Iapp))*Cell.Const.T #Cell Temperature
    t = 0:(1.0/i):((length(Iapp)-1)/i)
    
    #Simulate Model
    return Sim_Model(Cell,Iapp,Tk,SList,SOC,A,B,C,D,t)
end

#---------- Hankel Formation & SVD -----------------#
function fh!(H,Hlen1,Hlen2,puls,M,Puls_L)
    for lp1 in 1:length(Hlen2), lp2 in 1:length(Hlen1)
            H[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] .= @view puls[:,Hlen2[lp1]+Hlen1[lp2]]
    end

    #Init = abs.(convert(Vector{float(eltype(H))}, randn(size(H,1))))
    Init = convert(Vector{float(eltype(H))}, ones(size(H,1)))
    U,S,V = tsvd(H, M; initvec=Init, maxiter=1e4)
    println("---")  

    for lp1 in 1:length(Hlen2), lp2 in 1:length(Hlen1)
           H[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] .= @view puls[:,Hlen2[lp1]+Hlen1[lp2]+1]
    end
    return U,S,V'
end

"""
    D_Linear(Cell,ν_neg,ν_pos,σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep) 

    Function to linearise D array from input cell type and conductivites. 

"""
function D_Linear(Cell, ν_neg, ν_pos, σ_eff_Neg, κ_eff_Neg, σ_eff_Pos, κ_eff_Pos, κ_eff_Sep)
    D = Array{Float64}(undef,0,1)
    Dt = Array{Float64}(undef,0,1)
    D_ = Array{Float64}(undef,size(Cell.Transfer.Locs[2],1),1)
    q = Int64(1)
    tfs = Cell.Transfer.tfs
    Elec = Cell.Transfer.Elec
    Locs = Cell.Transfer.Locs
    for i in 1:size(tfs,1)
        if tfs[i] == C_e
            Dt =  zeros(length(Locs[i]))
        elseif tfs[i] == Phi_e
            Dt = Array{Float64}(undef,0,1)
            for pt in Locs[i]
                if pt <= Cell.Neg.L+eps()
                D_[q]  =  @. (Cell.Neg.L*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg*pt/Cell.Neg.L)) - pt*ν_neg*sinh(ν_neg) + Cell.Neg.L*(cosh(ν_neg)-cosh(ν_neg*(Cell.Neg.L-pt)/Cell.Neg.L)))/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg) #Lee. Eqn. 4.22 @ ∞
                elseif pt <= Cell.Neg.L + Cell.Sep.L + eps()
                D_[q] = @. (Cell.Neg.L - pt)/(Cell.Const.CC_A*κ_eff_Sep) + (Cell.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) #Lee. Eqn. 4.23 @ ∞
                else
                D_[q] = @. -Cell.Sep.L/(Cell.Const.CC_A*κ_eff_Sep) + Cell.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg)/(Cell.Const.CC_A*(σ_eff_Neg+ κ_eff_Neg)*ν_neg) + (Cell.Pos.L*(-σ_eff_Pos*cosh(ν_pos) + σ_eff_Pos*cosh((Cell.Const.Ltot-pt)*ν_pos/Cell.Pos.L) +  κ_eff_Pos*(cosh((pt-Cell.Neg.L-Cell.Sep.L)*ν_pos/Cell.Pos.L)-1)) - (pt-Cell.Neg.L-Cell.Sep.L)*κ_eff_Pos*sinh(ν_pos)*ν_pos)/(Cell.Const.CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_pos*sinh(ν_pos))
                end
                Dt = [Dt; D_[q]]
                q = q+1
            end

        elseif tfs[i] == C_se
            Dt = zeros(length(Locs[i]))
        elseif Elec[i] == "Pos"
            σ_eff = σ_eff_Pos
            κ_eff = κ_eff_Pos
            if tfs[i] == Phi_se
                Dt = @. -1*Cell.Pos.L/(Cell.Const.CC_A*ν_pos*sinh(ν_pos))*((1/κ_eff)*cosh(ν_pos*Locs[i])+(1/σ_eff)*cosh(ν_pos*(Locs[i]-1))) # Contribution to D as G->∞
            elseif tfs[i] == Flux
                Dt = @. -1*ν_pos*(σ_eff*cosh(ν_pos*Locs[i])+κ_eff*cosh(ν_pos*(Locs[i]-1)))/(Cell.Pos.as*F*Cell.Pos.L*Cell.Const.CC_A*(κ_eff+σ_eff)*sinh(ν_pos))
            elseif tfs[i] == Phi_s
                Dt = @. -1*(-Cell.Pos.L*(κ_eff*(cosh(ν_pos)-cosh(Locs[i]-1)*ν_pos))-Cell.Pos.L*(σ_eff*(1-cosh(Locs[i]*ν_pos)+Locs[i]*ν_pos*sinh(ν_pos))))/(Cell.Const.CC_A*σ_eff*(κ_eff+σ_eff)*ν_pos*sinh(ν_pos)) # Contribution to D as G->∞ 
            end
        elseif Elec[i] == "Neg"
            σ_eff = σ_eff_Neg
            κ_eff = κ_eff_Neg
            if tfs[i] == Phi_se
                Dt = @. Cell.Neg.L/(Cell.Const.CC_A*ν_neg*sinh(ν_neg))*((1/κ_eff)*cosh(ν_neg*Locs[i])+(1/σ_eff)*cosh(ν_neg*(Locs[i]-1))) # Contribution to D as G->∞
            elseif tfs[i] == Flux
                Dt = @. ν_neg*(σ_eff*cosh(ν_neg*Locs[i])+κ_eff*cosh(ν_neg*(Locs[i]-1)))/(Cell.Neg.as*F*Cell.Neg.L*Cell.Const.CC_A*(κ_eff+σ_eff)*sinh(ν_neg))
            elseif tfs[i] == Phi_s
                Dt = @. (-Cell.Neg.L*(κ_eff*(cosh(ν_neg)-cosh(Locs[i]-1)*ν_neg))-Cell.Neg.L*(σ_eff*(1-cosh(Locs[i]*ν_neg)+Locs[i]*ν_neg*sinh(ν_neg))))/(Cell.Const.CC_A*σ_eff*(κ_eff+σ_eff)*ν_neg*sinh(ν_neg)) # Contribution to D as G->∞ 
            end
        end
        D = [D; Dt]
    end
    return D
end


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


"""
    Construct(::String) 

    Function to create Construct dictionary.

    Currently supports:

    1. Doyle '94 parameterisation
    2. Chen 2020 parameterisation

"""
function Construct(CellType::String)
    if CellType == "Doyle 94"
        CellType = string(CellType,".jl")
        include(joinpath(dirname(pathof(LiiBRA)), "Data/Doyle_94", CellType))
    elseif CellType == "LG M50"
        CellType = string("LG_M50.jl")
        include(joinpath(dirname(pathof(LiiBRA)), "Data/Chen_2020", CellType))
    end
    return Cell
end


"""
    tuple_len(::NTuple) 

    Function to return Tuple length. 

"""
tuple_len(::NTuple{N, Any}) where {N} = N #Tuple Size



"""
    interp(MTup::Tuple,SList::Array,SOC) 

    Function to interpolate Tuple indices

"""
function interp(MTup::Tuple,SList::Array,SOC)
    T1 = 0
    T2 = 0
    for i in 1:length(SList)-1
        if SList[i] > SOC >= SList[i+1]
            T2 = i
            T1 = i+1
        end
    end
    return M =  @. MTup[T1]+(MTup[T2]-MTup[T1])*(SOC-SList[T1])/(SList[T2]-SList[T1])
end

function mag!(γ::AbstractArray)
    ψ = Array{Float64}(undef,size(γ))
    for j in 1:size(γ,2), i in 1:size(γ,1)
        if real(γ[i,j]) < 0
            ψ[i,j] = -abs(γ[i,j])
        else
            ψ[i,j] = abs(γ[i,j])
        end
    end
    return ψ
end

function mag!(γ::AbstractVector)
    ψ = Vector{Float64}(undef,length(γ))
    for i in 1:length(γ)
        if real(γ[i]) < 0
            ψ[i] = -abs(γ[i])
        else
            ψ[i] = abs(γ[i])
        end
    end
    return ψ
end

end # module