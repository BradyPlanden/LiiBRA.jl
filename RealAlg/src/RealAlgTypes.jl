struct Constants
    F
    R
    T
    t_plus
    De
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
    eps_s
    eps_e
    brug_De
end

struct Positive
    Rs::Int64
    Ds
    eps_s
    eps_e
    brug_De
end

struct Seperator
    eps_e
    brug_De
end

struct Cell
    Geo::Geometry
    Const::Constants
    Neg::Negative
    Pos::Positive
    Sep::Seperator
end