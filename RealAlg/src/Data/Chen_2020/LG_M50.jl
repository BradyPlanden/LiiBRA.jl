using Parameters

@with_kw mutable struct Constants
    T::Float64 = 273.15
    T_ref::Float64 = 298.15
        t_plus::Float64 = 0.2594 #Inital Transference Number
        tpf::Function = ce -> -0.1287*ce^3+0.4106*ce^2-0.4717*ce+0.4492 #Transference Number Function
        De::Float64 = 1.0e-11   #Inital Electrolyte Diffusivity
        Def::Function = ce -> 8.794e-11*ce^2-3.972e-10*ce+4.862e-10 #Electrolyte Diffusivity Function
    SOC::Float64 = 0.0
    ce0::Float64 = 2000
    dln::Float64 = 3.0
    Ea_κ::Float64 = 0.0
    Ea_De::Float64 = 0.0
        Lnegsep::Float64 = 0.0  #Init Value
        CC_A::Float64 = 0.1027  #Electrode Plate Area 
    κ::Float64 = 1.0
    κf::Function = ce -> 4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24
    Uocp::Function = (Electrode, θ) ->
        if Electrode == "Neg"
            Uocp = 1.97938*2.7182818284^(-39.3631*θ) + 0.2482 - 0.0909*tanh(29.8538*(θ - 0.1234)) - 0.04478*tanh(14.9159*(θ - 0.2769)) - 0.0205*tanh(30.4444*(θ - 0.6103))
        else
            Uocp = -0.8090*ce + 4.4875 - 0.0428*tanh(18.5138*(ce - 0.5542)) - 17.7326*tanh(15.7890*(ce - 0.3117)) + 17.5842*tanh(15.9308*(ce - 0.3120))
        end

    ∂Uocp::Function = (Electrode,θ) -> 
        if Electrode == "Neg"
            ∂Uocp = 0.667934002(tanh(14.9159θ - 4.13021271)^2) + 2.71371042(tanh(29.8538θ - 3.68395892)^2) + 0.6241102(tanh(30.4444θ - 18.58021732)^2) - 4.005754622 - (77.91453287630758(2.7182818284^(-39.3631θ)))
        else
            ∂Uocp = 279.9800214(tanh(15.789t - 4.9214313)^2) + 0.79239064(tanh(18.5138t - 10.26034796)^2) - 1.4510386800000594 - (280.13037336(tanh(15.9308t - 4.9704096)^2))
        end
end

@with_kw mutable struct Negative
        L::Float64 = 85.2e-6
        Rs::Float64 = 5.86e-6   # Particle radius [m]
        Ds::Float64 = 1.74e-15   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
        ϵ_s::Float64 = 0.75     #Active Material Volume Fraction
        ϵ_e::Float64 = 0.25    # Porosity of negative electrode
        De_brug::Float64 = 1.5
        κ_brug::Float64 = 1.5
        σ::Float64 = 215
        σ_brug::Float64 = 1.5
        θ_100::Float64 = 0.9014
        θ_0::Float64 = 0.0279
        cs_max::Float64 = 29583
        α::Float64 = 0.5
        k_norm::Float64 = 6.48e-7
    RFilm::Float64 = 0.0
    k_ref::Float64 = 2.16e-11
    D1::Float64 = 1.0   #Init Value
    D1f::Function = De -> De * ϵ_e^De_brug #Effective Diffusivity
    as::Float64 = 3*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Positive
        L::Float64 = 75.6e-6
        Rs::Float64 = 5.22e-6    # Particle radius [m]
        Ds::Float64 = 1.48e-15   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
        ϵ_s::Float64 = 66.5     #Active Material Volume Fraction
        ϵ_e::Float64 = 0.335    # Porosity of positive electrode
        De_brug::Float64 = 1.5
        κ_brug::Float64 = 1.5
        σ::Float64 = 0.18
        σ_brug::Float64 = 1.5
        θ_100::Float64 = 0.2661
        θ_0::Float64 = 0.9084
        cs_max::Float64 = 51765
        α::Float64 = 0.5
        k_norm::Float64 = 3.42e-6
    RFilm::Float64 = 0
    D3::Float64 = 1.0   #Init Value
    D3f::Function = De -> De * ϵ_e^De_brug
    as::Float64 = 3*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Seperator
        L::Float64 = 12e-6
        ϵ_e::Float64 = 0.47    # Porosity of separator
        De_brug::Float64 = 1.5
        κ_brug::Float64 = 1.5
        D2::Float64 = 1.0   #Init Value
        D2f::Function = De -> De * ϵ_e^De_brug
end

@with_kw mutable struct RealisationAlgorthim
    Fs::Float64 = 2
    SamplingT::Float64 = 1
    M::Int64 = 5 # Model order
    N::Int64 = 1 # Inputs
    Tlen::Int64 = 65536 #1048576 #2097152 #262144 #32768 #24
    H1::Array{Int64,1} = 0:2000 #4000 #4612
    H2::Array{Int64,1} = 0:2000 #4000 #4612
    Outs::Int64 = 21
end

@with_kw mutable struct Cell
    Const::Constants
    Neg::Negative
    Pos::Positive
    Sep::Seperator
    RA::RealisationAlgorthim
end

CellData = Cell(Constants(),Negative(),Positive(),Seperator(),RealisationAlgorthim())