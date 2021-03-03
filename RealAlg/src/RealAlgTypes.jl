using Parameters

@with_kw struct Constants
    T::Float64 = 298.15
    t_plus::Float64 = 0.363
    De::Float64 = 7.5e-11
    κ::Float64 = 1
    Init_SOC::Float64 = 0.5
    ce0::Float64 = 2000
    dln::Float64 = 3.0
    Ea_κ::Float64 = 0.0
end

@with_kw struct Geometry
    Lsep::Float64 = 7.6e-5
    Ltot::Float64 = 0.000394
    CC_A::Float64 = 1.0
end

@with_kw struct Negative
    L::Float64 = 0.000128
    Rs::Float64 = 1.25e-5
    Ds::Float64 = 3.9e-14
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
    ϵ_s::Float64 = 0.471
    ϵ_e::Float64 = 0.357
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
    σ::Float64 = 3.8
    σ_brug::Float64 = 1.0
    θ_max::Float64 = 0.762
    θ_min::Float64 = 0.17
    cs_max::Float64 = 22860
    α::Float64 = 0.5
    k_norm::Float64 = 2.21e-5
    RFilm::Float64 = 0.0
    k_ref::Float64 = 2.16e-11

end

@with_kw struct Positive
    L::Float64 = 0.00019
    Rs::Float64 = 1.25e-5
    Ds::Float64 = 1.0e-13
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
    ϵ_s::Float64 = 0.471
    ϵ_e::Float64 = 0.357
    De_brug::Float64 = 1.6
    κ_brug::Float64 = 1.5
    σ::Float64 = 100
    σ_brug::Float64 = 1.0
    θ_max::Float64 = 0.53
    θ_min::Float64 = 0.05
    cs_max::Float64 = 26390
    α::Float64 = 0.5
    k_norm::Float64 = 2.28e-05
    RFilm::Float64 = 1.94e-11
    
end

@with_kw struct Seperator
    ϵ_e::Float64 = 0.724
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
end

@with_kw struct RealisationAlgorthim
    Ts::Float64 = 0.5
    M::Int64 = 4
    Tlen::Float64 = 1.5
end

# struct kappa{T<:Number}
#     ce::T
#     function kappa(ce::T) where T
#         return 4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24
#     end
#     new{typeof(ce)}(ce)
# end

struct kappa{T<:Number}
    κ::T
    function kappa(ce::T) where T
        κ  =  4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24
        new{typeof(κ)}(κ)
   end
end


@with_kw struct Cell
    Const::Constants
    Geo::Geometry
    Neg::Negative
    Pos::Positive
    Sep::Seperator
    RA::RealisationAlgorthim
    
end

@with_kw struct FCalls
    Kap::kappa
end
