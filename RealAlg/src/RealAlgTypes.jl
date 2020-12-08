struct Constants
    T
    t_plus
    De
    κ
    Init_SOC
    ce0
end

struct Geometry
    Lneg
    Lpos
    Lsep
    Ltot
    CC_A
end

struct Negative
    Rs::Int64
    Ds
    ϵ_s
    ϵ_e
    De_brug
    κ_brug
    σ
    σ_brug
    θ_max
    θ_min
    cs_max
    α
    k_norm
    RFilm
end

struct Positive
    Rs::Int64
    Ds
    ϵ_s
    ϵ_e
    De_brug
    κ_brug
    σ
    σ_brug
    θ_max
    θ_min
    cs_max
    α
    k_norm
    RFilm
end

struct Seperator
    ϵ_e
    De_brug
end

struct Cell
    Const::Constants
    Geo::Geometry
    Neg::Negative
    Pos::Positive
    Sep::Seperator
end