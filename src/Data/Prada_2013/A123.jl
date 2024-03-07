using Parameters
"""
Parameters for an LFP cell, from the paper :footcite:t:`Prada2013`
"""

@with_kw mutable struct Constants
    CellTyp::String = "A123"                                                    # LFP, Prada (2013) A123 Systems ANR26650 m1
    T::Float64 = 298.15                                                         # Cell Temperature (K)
    T_ref::Float64 = 298.15                                                     # Reference Temperature (K)
    t_plus::Float64 = 0.36                                                      # Inital Transference Number
    tpf::Function = ce -> -0.1287 * (ce / 1000)^3 + 0.4106 * (ce / 1000)^2 -
                          0.4717 * (ce / 1000) + 0.4492                         # Transference Number Function 
    De::Float64 = 2.0e-10                                                       # Inital Electrolyte Diffusivity (J mol⁻¹)
    Def::Function = ce -> 8.794e-11 * (ce / 1000)^2 - 3.972e-10 *
                                                      (ce / 1000) + 4.862e-10   # CHECK Electrolyte Diffusivity Function 
    SOC::Float64 = 1.0                                                          # Initial State of Charge
    ce0::Float64 = 1200                                                         # Initial Electrolyte Concentration (mol m⁻³)
    Vmin::Float64 = 2.0                                                         # Min Cell Voltage
    Ea_κ::Float64 = 0.0 # check 
    Ea_De::Float64 = 0 # check
    CC_A::Float64 = 0.6*0.3                                                     # Electrode Plate Area (m²)                                           
    κ::Float64 = 0.9487                                                         # Electrolyte Conductivity
    # κf::Function = ce -> 0.1297 * (ce / 1000)^3 - 2.51 * (ce / 1000)^1.5 +
    #                      3.329 * (ce / 1000) 
    κf::Function = ce -> (4.1253e-4 + 5.007 * (ce / 1e6) -
                          4721.2 * (ce / 1e6)^2 +
                          1.5094e6 * (ce / 1e6)^3 -
                          1.6018e8 * (ce / 1e6)^4) * 1e3     # Electrolyte Conductivity Function(S/m)
    Uocp::Function = (Electrode, θ) -> if Electrode == "Neg"
        Uocp = @. 1.97938 * 2.7182818284 * exp(-39.3631 * θ) + 0.2482 -
                  0.0909 * tanh(29.8538 * (θ - 0.1234)) -
                  0.04478 * tanh(14.9159 * (θ - 0.2769)) -
                  0.0205 * tanh(30.4444 * (θ - 0.6103))
    else
        Uocp = @. 3.4077 - 0.020269 * θ + 0.5 * exp(-150*θ) - 0.9 * exp(-30*(1-θ))
    end
    ∂Uocp::Function = (Electrode, θ) -> if Electrode == "Neg"
        ∂Uocp = @. -0.62411 * ((sech(18.5802 - 30.4444 * θ))^2 +
                    4.34813 * (sech(3.68396 - 29.8538 * θ))^2 +
                    1.07022 * (sech(4.13021 - 14.9159 * θ))^2) - 211.794 * exp(-39.3631 * θ)
    else
        ∂Uocp = @. - 0.020269 - 75 * exp(-150 * θ) - 2.52656e-12 * exp(30 * θ)
    end
    Ce_M::Int64 = 1                             # Over-written
    D1::Float64 = 1.0
    D2::Float64 = 1.0
    D3::Float64 = 1.0
    Ltot::Float64 = 0.0
    Lnegsep::Float64 = 0.0
    CeRootRange::Float64 = 20
    Debug::Bool = false
end

# CHECK
@with_kw mutable struct Negative
    L::Float64 = 3.4e-5                         # Negative Electrode Length
    Rs::Float64 = 5e-6                          # Particle Radius [m]
    Ds::Float64 = 3e-15                         # Solid Diffusivity [m^2/s]
    Ea_σ::Float64 = 0                           # Activation Energy Solid Conductivity
    Ea_Ds::Float64 = 0                          # Activation Energy Solid Diffusivity
    ϵ_s::Float64 = 0.58                         # Active Material Volume Fraction
    ϵ_e::Float64 = 0.36                         # Porosity of Negative Electrode
    De_brug::Float64 = 1.5                      # Bruggeman Diffusion Exponent
    κ_brug::Float64 = 1.5                       # Bruggeman Electrolyte Conductivity Exponent
    σ::Float64 = 215                            # Solid Phase Conductivity
    σ_brug::Float64 = 1.5                       # Bruggeman Solid Conductivity Exponent
    θ_100::Float64 = 0.8100363695199028         # Theta @ 100% Lithium Concentration
    θ_0::Float64 = 0.017619411644243348         # Theta @ 0% Lithium Concentration
    cs_max::Float64 = 30555                     # Max Electrode Concentration
    α::Float64 = 0.5                            # Alpha Factor
    k_norm::Float64 = 6.716047e-12              # check Reaction Rate
    Ea_κ::Float64 = 35000                       # check Activation Energy
    RFilm::Float64 = 0                          # Film Resistance - Ωm²
    D1::Float64 = 1.0                           # Init Value
    D1f::Function = De -> De * ϵ_e^De_brug      # Effective Diffusivity
    as::Float64 = 3.0 * ϵ_s / Rs                # Specific Interfacial Surf. Area
    β::Array{ComplexF64} = [0 - 0.0im]          # Init Value
    β!::Function = (s) -> @. Rs * sqrt(s / Ds)  # Negative β

end

@with_kw mutable struct Positive
    L::Float64 = 8e-5                           # Positive Electrode Length
    Rs::Float64 = 5e-08                         # Particle radius [m]
    Ds::Float64 = 5.9e-18                       # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0                           # Activation Energy Solid Conductivity CHECK
    Ea_Ds::Float64 = 0                          # Activation Energy Solid Diffusivity CHECK
    ϵ_s::Float64 = 0.374                        # Active Material Volume Fraction
    ϵ_e::Float64 = 0.426                        # Porosity of positive electrode
    De_brug::Float64 = 1.5                      # Bruggeman Diffusivity Exponent
    κ_brug::Float64 = 1.5                       # Bruggeman Electrolyte Conductivity Exponent
    σ::Float64 = 0.33795074                     # Solid Phase Conductivity
    σ_brug::Float64 = 1.5                       # Bruggeman Solid Conductivity Exponent
    θ_100::Float64 = 0.003767884398705004       # Theta @ 100% Lithium Concentration
    θ_0::Float64 = 0.703500714165226            # Theta @ 0% Lithium Concentration
    cs_max::Float64 = 22806                     # Max Electrode Concentration
    α::Float64 = 0.5                            # Alpha Factor
    k_norm::Float64 = 3.54458e-11               # check Reaction Rate CHECK
    Ea_κ::Float64 = 17800                       # check Activation Energy CHECK
    RFilm::Float64 = 0                          # Film Resistance - Ωm²
    D3::Float64 = 1                             # Init Value CHECK
    D3f::Function = De -> De * ϵ_e^De_brug      # Effective Diffusivity
    as::Float64 = 3.0 * ϵ_s / Rs                # Specific interfacial surf. area
    β::Array{ComplexF64} = [0 - 0.0im]          # Init Value
    β!::Function = (s) -> @. Rs * sqrt(s / Ds)  # Positive β
end

@with_kw mutable struct Separator
    L::Float64 = 2.5e-5                         # Separator Length
    ϵ_e::Float64 = 0.45                         # Porosity of separator
    De_brug::Float64 = 1.5                      # Bruggeman Diffusivity Factor
    κ_brug::Float64 = 1.5                       # Bruggeman Electrolyte Conductivity Factor
    D2::Float64 = 1                             # Init Value
    D2f::Function = De -> De * ϵ_e^De_brug      # Effective Diffusivity
end

@with_kw mutable struct Realisation
    Fs::Float64 = 1                             # Sampling Frequency of Transfer Functions [Hz]
    SamplingT::Float64 = 1                      # Final Model Sampling Time [s]
    M::Int64 = 4                                # Model Order
    N::Int64 = 1                                # Number of Inputs
    Tlen::Int64 = 16200                         # Transfer Function Response Length [s]
    H1::Array{Int64, 1} = [1:2500; 3000:3500]   # Hankel Dimensions 1
    H2::Array{Int64, 1} = [1:2500; 3000:3500]   # Hankel Dimensions 2
    Outs::Int64 = 1                             # Number of Outputs (Rewritten)
    Nfft::Int64 = 0
    f::UnitRange{Int64} = 0:1
    s::Array{ComplexF64} = [0 - 0.0im]
    Nfft!::Function = (Fs, Tlen) -> ceil(2^(log2(Fs * Tlen)))
    f!::Function = (Nfft) -> 0:(Nfft - 1)
    s!::Function = (Fs, Nfft, f) -> transpose(((2im .* Fs) * tan.(pi .* f ./ Nfft)))
end

@with_kw mutable struct TransferFun
    # Number of Spatial Points in Electrolyte and Solid
    Sₑ::Int = 4
    Sₛ::Int = 2

    tfs::Vector{Function} = [
        C_e,
        Phi_e,
        C_se,
        Phi_s,
        Phi_se,
        Flux,
        C_se,
        Phi_s,
        Flux,
        Phi_se,
    ]

    Elec::Vector{String} = [
        "Na",
        "Na",
        "Pos",
        "Pos",
        "Pos",
        "Pos",
        "Neg",
        "Neg",
        "Neg",
        "Neg",
    ]

    Locs::Function = (Sₑ, Sₛ) -> if Sₑ == 6 && Sₛ == 4
        return [
            Float64[0, 4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
        ]
    elseif Sₑ == 4 && Sₛ == 4
        return [
            Float64[0.0, 8.52e-5, 9.72e-5, 1.728e-4],
            Float64[8.52e-5, 9.72e-5, 1.728e-4],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
        ]
    elseif Sₑ == 8 && Sₛ == 4
        return [
            Float64[0, 2.84e-5, 5.68e-5, 8.52e-5, 9.72e-5, 1.224e-4, 1.476e-4, 1.728e-4],
            Float64[2.84e-5, 5.68e-5, 8.52e-5, 9.72e-5, 1.224e-4, 1.476e-4, 1.728e-4],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
        ]
    elseif Sₑ == 9 && Sₛ == 6
        return [
            Float64[0, 2.84e-5, 5.68e-5, 8.52e-5, 9.12e-5, 9.72e-5, 0.0001224, 0.0001476,
                    1.728e-4],
            Float64[2.84e-5, 5.68e-5, 8.52e-5, 9.12e-5, 9.72e-5, 0.0001224, 0.0001476,
                    1.728e-4],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[1],
            Float64[0, 1 / 3, 2 / 3, 1],
            Float64[0, 1 / 3, 2 / 3, 1],
        ]
    elseif Sₑ == 4 && Sₛ == 2
        [
            Float64[0, 8.52E-05, 9.72E-05, 1.728E-04],
            Float64[8.52E-05, 9.72E-05, 1.728E-04],
            Float64[0, 1],
            Float64[1],
            Float64[0, 1],
            Float64[0, 1],
            Float64[0, 1],
            Float64[1],
            Float64[0, 1],
            Float64[0, 1],
        ]
    elseif Sₑ == 4 && Sₛ == 6
        return [
            Float64[0.00, 4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
        ]
    elseif Sₑ == 6 && Sₛ == 2
        return [
            Float64[0, 4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[0, 1],
            Float64[1],
            Float64[0, 1],
            Float64[0, 1],
            Float64[0, 1],
            Float64[1],
            Float64[0, 1],
            Float64[0, 1],
        ]
    elseif Sₑ == 6 && Sₛ == 6
        return [
            Float64[0.00, 4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[4.26e-5, 8.52e-5, 9.72e-5, 1.35e-04, 1.728e-4],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
            Float64[0, 0.2, 0.4, 0.6, 0.8, 1],
        ]
    end
end

@with_kw mutable struct Params
    Const::Constants
    Neg::Negative
    Pos::Positive
    Sep::Separator
    RA::Realisation
    Transfer::TransferFun
end

return Params(Constants(), Negative(), Positive(), Separator(), Realisation(),
              TransferFun())