using Parameters

@with_kw mutable struct Constants
    CellTyp::String = "LG M50"
    T::Float64 = 298.15 # Cell Temperature (K)
    T_ref::Float64 = 298.15 # Reference Temperature (K)
    t_plus::Float64 = 0.2594 # Inital Transference Number
    tpf::Function = ce -> -0.1287*(ce/1000)^3+0.4106*(ce/1000)^2-0.4717*(ce/1000)+0.4492 # Transference Number Function - Requires ce in dm-3
    De::Float64 = 1.769e-10   # Inital Electrolyte Diffusivity (J mol⁻¹)
    Def::Function = ce -> 8.794e-11*(ce/1000)^2-3.972e-10*(ce/1000)+4.862e-10 # Electrolyte Diffusivity Function - Requires ce in dm-3
    SOC::Float64 = 0.8 # Initial State of Charge
    ce0::Float64 = 1000 # Initial Electrolyte Concentration (mol m⁻³)
    Ea_κ::Float64 = 0.
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
    Ce_M::Int64 = 1     #(Rewritten)
    D1::Float64 = 1.0
    D2::Float64 = 1.0
    D3::Float64 = 1.0
    Ltot::Float64 = 0.
    Lnegsep::Float64 = 0.
    CeRootRange::Float64 = 10
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
    θ_100::Float64 = 0.910612 #0.9014 # Theta @ 100% Lithium Concentration
    θ_0::Float64 = 0.0263473 #0.0279 # Theta @ 0% Lithium Concentration
    cs_max::Float64 = 33133 # Max Electrode Concentration
    α::Float64 = 0.5    # Alpha Factor
    k_norm::Float64 = 6.716047e-12 #6.48e-7 #2.12e-10#1E-5 #7.226781e-7 # Reaction Rate
    Ea_κ::Float64 = 35000   # Activation Energy
    RFilm::Float64 = 0 # Film Resistance - Ωm²
    D1::Float64 = 1.   # Init Value
    D1f::Function = De -> De * ϵ_e^De_brug # Effective Diffusivity
    as::Float64 = 3.0*ϵ_s/Rs # Specific Interfacial Surf. Area
    β::Array{ComplexF64} = [0 - 0.0im] # Init Value
    β!::Function = (s) -> @. Rs*sqrt(s/Ds) # Negative β
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
    σ::Float64 = 0.18 #0.847  # Solid Phase Conductivity
    σ_brug::Float64 = 1.5   # Bruggeman Solid Conductivity Exponent
    θ_100::Float64 = 0.263849 #0.27 # Theta @ 100% Lithium Concentration
    θ_0::Float64 = 0.853974 #0.9084 # Theta @ 0% Lithium Concentration
    cs_max::Float64 = 63104 # Max Electrode Concentration
    α::Float64 = 0.5    # Alpha Factor
    k_norm::Float64 = 3.54458e-11 #3.42e-6 #1.12e-9#4E-05 #7.264272e-6 # Reaction Rate
    Ea_κ::Float64 = 17800   # Activation Energy
    RFilm::Float64 = 0 # Film Resistance - Ωm²
    D3::Float64 = 1   # Init Value
    D3f::Function = De -> De * ϵ_e^De_brug
    as::Float64 = 3.0*ϵ_s/Rs # Specific interfacial surf. area
    β::Array{ComplexF64} = [0 - 0.0im] # Init Value
    β!::Function = (s) -> @. Rs*sqrt(s/Ds) # Positive β
end

@with_kw mutable struct Seperator
    L::Float64 = 12e-6  # Seperator Length
    ϵ_e::Float64 = 0.47    # Porosity of separator
    De_brug::Float64 = 1.5  # Bruggeman Diffusivity Factor
    κ_brug::Float64 = 1.5   # Bruggeman Electrolyte Conductivity Factor
    D2::Float64 = 1   # Init Value
    D2f::Function = De -> De * ϵ_e^De_brug
end

@with_kw mutable struct Realisation
    Fs::Float64 = 4    # Sampling Frequency of Transfer Functions [Hz]
    SamplingT::Float64 = 0.25     # Final Model Sampling Time [s]
    M::Int64 = 4    # Model Order
    N::Int64 = 1    # Number of Inputs
    Tlen::Int64 = 16200 #Transfer Function Response Length [s] (Change to min)
    H1::Array{Int64,1} = 1:2500 #4000 #4612     # Hankel Dimensions 1
    H2::Array{Int64,1} = 1:2500 #4000 #4612     # Hankel Dimensions 2
    Outs::Int64 = 1    # Number of Outputs (Rewritten)
    Nfft::Int64 = 0
    f::UnitRange{Int64} = 0:1
    s::Array{ComplexF64} = [0 - 0.0im]
    Nfft!::Function = (Fs,Tlen) -> ceil(2^(log2(Fs*Tlen))) #2^(ceil(log2(Fs*Tlen))) Old way (constains to 2^ values)
    f!::Function = (Nfft) -> 0:Nfft-1
    s!::Function = (Fs,Nfft,f) -> transpose(((2im.*Fs)*tan.(pi.*f./Nfft)))
end

@with_kw mutable struct TransferFun
    #Electrolyte - Solid

    #6 - 4
    #tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    #Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    #Locs::Vector{Vector{Float64}} = [Float64[0,4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1]]
   
    #4 - 4
    #tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    #Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    #Locs::Vector{Vector{Float64}} =  [Float64[0.0, 8.52e-5,9.72e-5,0.0001728], Float64[8.52e-5,9.72e-5,0.0001728], Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1]]

    #8 - 4
    #tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    #Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    #Locs::Vector{Vector{Float64}} = [Float64[0,2.84e-5,5.68e-5,8.52e-5,9.72e-5,1.224e-4,1.476e-4,1.728e-4], Float64[2.84e-5,5.68e-5,8.52e-5,9.72e-5,1.224e-4,1.476e-4,1.728e-4], Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1]]
    
    # 9 - 6
    #tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    #Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    #Locs::Vector{Vector{Float64}} = [Float64[0,2.84e-5,5.68e-5,8.52e-5,9.12e-5,9.72e-5,0.0001224,0.0001476,0.0001728], Float64[2.84e-5,5.68e-5,8.52e-5,9.12e-5,9.72e-5,0.0001224,0.0001476,0.0001728], Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1],Float64[1],Float64[0,0.333,0.666,1],Float64[0,0.333,0.666,1]]
    
    # 4 - 2
    tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    Locs::Vector{Vector{Float64}} = [Float64[0, 8.52E-05, 9.72E-05, 1.728E-04], Float64[8.52E-05, 9.72E-05, 1.728E-04], Float64[0,1], Float64[1], Float64[0,1], Float64[0,1], Float64[0,1], Float64[1], Float64[0,1], Float64[0,1]]

    # 4 - 6
    #tfs =   [[C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se] ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] [Float64[0.00,4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[0,0.2,0.4,0.6,0.8,1], Float64[1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[0,0.2,0.4,0.6,0.8,1]]]

    # 6 - 2
    #tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    #Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    #Locs::Vector{Vector{Float64}} = [Float64[0,4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[0,1], Float64[1], Float64[0,1], Float64[0,1], Float64[0,1], Float64[1], Float64[0,1], Float64[0,1]]

    # 6 - 6
    #tfs::Vector{Function} = [C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se]
    #Elec::Vector{String} = ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] 
    #Locs::Vector{Vector{Float64}} = [Float64[0.00,4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[4.26e-5,8.52e-5,9.72e-5,1.35e-04,1.728e-4], Float64[0,0.2,0.4,0.6,0.8,1], Float64[1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[1],Float64[0,0.2,0.4,0.6,0.8,1],Float64[0,0.2,0.4,0.6,0.8,1]]

end

@with_kw mutable struct Params
    Const::Constants
    Neg::Negative
    Pos::Positive
    Sep::Seperator
    RA::Realisation
    Transfer::TransferFun
end

Cell = Params(Constants(),Negative(),Positive(),Seperator(),Realisation(),TransferFun())
Cell.Const.Lnegsep, Cell.Const.Ltot = Cell.Neg.L+Cell.Sep.L,Cell.Neg.L+Cell.Sep.L+Cell.Pos.L
Cell.Const.D1 = Cell.Const.De*Cell.Neg.ϵ_e^Cell.Neg.De_brug
Cell.Const.D2 = Cell.Const.De*Cell.Sep.ϵ_e^Cell.Sep.De_brug
Cell.Const.D3 = Cell.Const.De*Cell.Pos.ϵ_e^Cell.Pos.De_brug
Cell.Const.Ce_M = size(Cell.Transfer.Locs[1],1)
Cell.RA.Outs = sum([size(Cell.Transfer.Locs[i],1) for i in 1:length(Cell.Transfer.tfs)])

