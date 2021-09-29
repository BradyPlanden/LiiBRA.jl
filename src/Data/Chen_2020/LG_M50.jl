using Parameters

@with_kw mutable struct Constants
    CellTyp::String = "Doyle_94"
    T::Float64 = 298.15 # Cell Temperature
    T_ref::Float64 = 298.15 # Reference Temperature
    t_plus::Float64 = 0.2594 # Inital Transference Number
    tpf::Function = ce -> -0.1287*(ce/1000)^3+0.4106*(ce/1000)^2-0.4717*(ce/1000)+0.4492 # Transference Number Function - Requires ce in dm-3
    De::Float64 = 1.769e-10   # Inital Electrolyte Diffusivity
    Def::Function = ce -> 8.794e-11*(ce/1000)^2-3.972e-10*(ce/1000)+4.862e-10 # Electrolyte Diffusivity Function - Requires ce in dm-3
    SOC::Float64 = 0.5 # Initial State of Charge
    ce0::Float64 = 1000 # Initial Electrolyte Concentration
    #dln::Float64 = 3.
    Ea_κ = 0.
    Ea_De::Float64 = 0
    CC_A::Float64 =  0.1027  #Electrode Plate Area 
    κ::Float64 = 0.9487
    κf::Function = ce -> 0.1297*(ce/1000)^3-2.51*(ce/1000)^1.5+3.329*(ce/1000) #Requires ce in dm-3
    Uocp::Function = (Electrode, θ) ->
        if Electrode == "Neg"
            Uocp = @. 1.97938*2.7182818284*exp(-39.3631*θ) + 0.2482 - 0.0909*tanh(29.8538*(θ - 0.1234)) - 0.04478*tanh(14.9159*(θ - 0.2769)) - 0.0205*tanh(30.4444*(θ - 0.6103))
        else
            Uocp = @. -0.8090*θ + 4.4875 - 0.0428*tanh(18.5138*(θ - 0.5542)) - 17.7326*tanh(15.7890*(θ - 0.3117)) + 17.5842*tanh(15.9308*(θ - 0.3120))
        end

    ∂Uocp::Function = (Electrode,θ) -> 
        if Electrode == "Neg"
            #∂Uocp = @. 0.667934002(tanh(14.9159*θ - 4.13021271)^2) + 2.71371042(tanh(29.8538*θ - 3.68395892)^2) + 0.6241102(tanh(30.4444*θ - 18.58021732)^2) - 4.005754622 - (77.91453287630758(2.7182818284^(-39.3631*θ)))
            ∂Uocp = @. -0.62411*((sech(18.5802 - 30.4444*θ))^2 + 4.34813*(sech(3.68396 - 29.8538*θ))^2 + 1.07022*(sech(4.13021 - 14.9159*θ))^2) - 211.794*exp(-39.3631*θ)
        else
            ∂Uocp = @. 279.9800214*(tanh(15.789*θ - 4.9214313)^2) + 0.79239064*(tanh(18.5138*θ - 10.26034796)^2) - 1.4510386800000594 - (280.13037336*(tanh(15.9308*θ - 4.9704096)^2))
        end
    Ce_M::Int64 = 6
    D1::Float64 = 1.0
    D2::Float64 = 1.0
    D3::Float64 = 1.0
    Ltot::Float64 = 0.
    Lnegsep::Float64 = 0.
end

@with_kw mutable struct Negative
    L::Float64 = 85.2e-6    # Negative Electrode Length
    Rs::Float64 = 5.86e-6   # Particle Radius [m]
    Ds::Float64 = 3.3e-14   # Solid Diffusivity [m^2/s]
    Ea_σ::Float64 = 0     # Activation Energy Solid Conductivity
    Ea_Ds::Float64 = 0    # Activation Energy Solid Diffusivity
    ϵ_s::Float64 = 0.75     # Active Material Volume Fraction
    ϵ_e::Float64 = 0.25     # Porosity of Negative Electrode
    De_brug::Float64 = 1.5  # Bruggeman Diffusion Exponent
    κ_brug::Float64 = 1.5   # Bruggeman Electrolyte Conductivity Exponent
    σ::Float64 = 215    # Solid Phase Conductivity
    σ_brug::Float64 = 1.5   # Bruggeman Solid Conductivity Exponent
    θ_100::Float64 = 0.9014 # Theta @ 100% Lithium Concentration
    θ_0::Float64 = 0.0279   # Theta @ 0% Lithium Concentration
    cs_max::Float64 = 33133 # Max Electrode Concentration
    α::Float64 = 0.5    # Alpha Factor
    k_norm::Float64 = 7.264265E-06 #6.48e-6 #7.226781E-06 #6.8973799e-13 #4.1580e-8 #2.12e-10 #Initial Reaction Rate
    Ea_κ::Float64 = 35000   # Activation Energy
    RFilm::Float64 = 0 # Film Resistance
    D1::Float64 = 1.   # Init Value
    D1f::Function = De -> De * ϵ_e^De_brug # Effective Diffusivity
    as::Float64 = 3.0*ϵ_s/Rs # Specific Interfacial Surf. Area
end

@with_kw mutable struct Positive
    L::Float64 = 75.6e-6    # Positive Electrode Length
    Rs::Float64 = 5.22e-6    # Particle radius [m]
    Ds::Float64 = 4e-15   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0
    Ea_Ds::Float64 = 0
    ϵ_s::Float64 = 0.665    # Active Material Volume Fraction
    ϵ_e::Float64 = 0.335    # Porosity of positive electrode
    De_brug::Float64 = 1.5  # Bruggeman Diffusivity Exponent
    κ_brug::Float64 = 1.5   # Bruggeman Electrolyte Conductivity Exponent
    σ::Float64 = 0.18   # Solid Phase Conductivity
    σ_brug::Float64 = 1.5   # Bruggeman Solid Conductivity Exponent
    θ_100::Float64 = 0.27 # Theta @ 100% Lithium Concentration
    θ_0::Float64 = 0.9084   # Theta @ 0% Lithium Concentration
    cs_max::Float64 = 63104 # Max Electrode Concentration
    α::Float64 = 0.5    # Alpha Factor
    k_norm::Float64 =  7.26426E-05 #9E-05 #3.640283886203905e-12 #3.5954e-7 #3.42e-6 #1.12e-9  #Initial Reaction Rate
    Ea_κ::Float64 = 17800   # Activation Energy
    RFilm::Float64 = 0 # Film Resistance
    D3::Float64 = 1   # Init Value
    D3f::Function = De -> De * ϵ_e^De_brug
    as::Float64 = 3.0*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Seperator
    L::Float64 = 12e-6  # Seperator Length
    ϵ_e::Float64 = 0.47    # Porosity of separator
    De_brug::Float64 = 1.5  # Bruggeman Diffusivity Factor
    κ_brug::Float64 = 1.5   # Bruggeman Electrolyte Conductivity Factor
    D2::Float64 = 1   # Init Value
    D2f::Function = De -> De * ϵ_e^De_brug
end

@with_kw mutable struct RealisationAlgorthim
    Fs::Float64 = 4    # Sampling Frequency of Transfer Functions [Hz]
    SamplingT::Float64 = 0.25     # Final Model Sampling Time [s]
    M::Int64 = 10    # Model Order
    N::Int64 = 1    # Number of Inputs
    Tlen::Int64 = 21600 #28800 #43200 #65536 #131072 #1048576 #2097152 #262144 #32768 #24    #Transfer Function Response Length [s] (Change to min)
    H1::Array{Int64,1} = 0:2500 #4000 #4612     # Hankel Dimensions 1
    H2::Array{Int64,1} = 0:2500 #4000 #4612     # Hankel Dimensions 2
    Outs::Int64 = 31    # Number of Outputs
end

@with_kw mutable struct TransferFun
    tfs =   [[C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se] ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] [Number[0, 4.26E-05, 8.52E-05, 9.72E-05, 1.35E-4, 1.728E-04], Number[4.26E-05, 8.52E-05, 9.72E-05, 1.35E-4, 1.728E-04], Number[0,0.5,1], Number[1],Number[0,0.5,1],Number[0,0.5,1],Number[0,0.5,1],Number[1],Number[0,0.5,1],Number[0,0.5,1]]]

end
@with_kw mutable struct Cell
    Const::Constants
    Neg::Negative
    Pos::Positive
    Sep::Seperator
    RA::RealisationAlgorthim
    Transfer::TransferFun
end

CellData = Cell(Constants(),Negative(),Positive(),Seperator(),RealisationAlgorthim(),TransferFun())
CellData.Const.Lnegsep, CellData.Const.Ltot = CellData.Neg.L+CellData.Sep.L,CellData.Neg.L+CellData.Sep.L+CellData.Pos.L
CellData.Const.D1 = CellData.Const.De*CellData.Neg.ϵ_e^CellData.Neg.De_brug
CellData.Const.D2 = CellData.Const.De*CellData.Sep.ϵ_e^CellData.Sep.De_brug
CellData.Const.D3 = CellData.Const.De*CellData.Pos.ϵ_e^CellData.Pos.De_brug
CellData.Const.Ce_M = size(CellData.Transfer.tfs[1,3],1)