using Parameters

@with_kw mutable struct Constants
    T::Float64 = 298.15
    T_ref::Float64 = 298.15
    t_plus::Float64 = 0.363
    De::Float64 = 7.5e-11
    SOC::Float64 = 1.
    ce0::Float64 = 2000
    dln::Float64 = 3.0
    Ea_κ::Float64 = 0.
    Ea_De::Float64 = 0.
    CC_A::Float64 = 1.0 
    k_ref::Float64 = 2.16e-11
    κ::Float64 = 1.0
    κf::Function = ce -> 4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24
    Uocp::Function = (Electrode, θ) ->
        if Electrode == "Neg"
            Uocp = @. -0.16+1.32*exp(-3*θ)+10*exp(-2000*θ)
        else
            Uocp = @. 4.19829+0.0565661*tanh(-14.5546*θ+8.60942)-0.0275479*(1/(0.998432-θ)^0.492465-1.90111)-0.157123*exp(-0.04738*θ^8)+0.810239*exp(-40*(θ-0.133875))
        end

    ∂Uocp::Function = (Electrode,θ) -> 
        if Electrode == "Neg"
            ∂Uocp = @. -20000*exp(-2000*θ)-3.96*exp(-3*θ)
        else
            ∂Uocp = @. -32.4096*exp(-40*(-0.133875+θ))-0.0135664/((0.998432-θ)^1.49247)+0.0595559*exp(-0.04738*θ^8)*θ^7-0.823297*(sech(8.60942-14.5546*θ))^2
        end
    Ce_M::Int64 = 4
end
@with_kw mutable struct Negative
    L::Float64 = 1.28e-4
    Rs::Float64 = 1.25e-5   # Particle radius [m]
    Ds::Float64 = 3.9e-14   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.
    Ea_Ds::Float64 = 0.
    ϵ_s::Float64 = 0.471
    ϵ_e::Float64 = 0.357    # Porosity of negative electrode
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
    σ_brug::Float64 = 1.
    σ::Float64 = 100
    θ_100::Float64 = 0.53
    θ_0::Float64 = 0.05
    cs_max::Float64 = 26390
    α::Float64 = 0.5
    k_norm::Float64 = 2.28421166145003e-05
    RFilm::Float64 = 0.
    D1::Float64 = 1.   #Init Value
    D1f::Function = De -> De * ϵ_e^De_brug #Effective Diffusivity
    as::Float64 = 3*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Positive
    L::Float64 = 1.9e-4
    Rs::Float64 = 8.5e-6    # Particle radius [m]
    Ds::Float64 = 1.0e-13   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0     # Activation Energy Solid Conductivity
    Ea_Ds::Float64 = 0.0    # Activation Energy Solid Diffusivity
    ϵ_s::Float64 = 0.297    # Solid Phase Volume Fraction
    ϵ_e::Float64 = 0.444    # Porosity of positive electrode
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
    σ_brug::Float64 = 1.0
    σ::Float64 = 3.8
    θ_100::Float64 = 0.17
    θ_0::Float64 = 0.762004800037954
    cs_max::Float64 = 22860
    α::Float64 = 0.5
    k_norm::Float64 = 2.20728263615611e-05
    RFilm::Float64 = 0
    D3::Float64 = 1.0   #Init Value
    D3f::Function = De -> De * ϵ_e^De_brug
    as::Float64 = 3*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Seperator
    L::Float64 = 7.60e-5
    ϵ_e::Float64 = 0.724    # Porosity of separator
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
    D2::Float64 = 1.0   #Init Value
    D2f::Function = De -> De * ϵ_e^De_brug
end

@with_kw mutable struct RealisationAlgorthim
    Fs::Float64 = 2     #Pulse Sampling Frequency
    SamplingT::Float64 = 1    #Final Sampling Time Step
    M::Int64 = 10 # Model order
    Tlen::Int64 =  131072 #65536 #36000 #65536 #14400 #131072 #1048576 #2097152 #262144 #32768 #24
    H1::Array{Int64,1} = 0:3000 #4000 #4612
    H2::Array{Int64,1} = 0:3000 #4000 #4612
    N::Int64 = 1 # Number of Inputs
    Outs::Int64 = 21    #Number of Outputs
end

@with_kw mutable struct TransferFun
    tfs =   [[C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se] ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] [Number[0, 128e-6, 204e-6, 394e-6], Number[128e-6, 204e-6, 394e-6], Number[0,1], Number[1],Number[0,1],Number[0,1],Number[0,1],Number[1],Number[0,1],Number[0,1]]]
    # tfst_temp = Array{String}(undef,0,1)
    # t1 = Array{String}(undef,0,1)
    # t2 = Array{String}(undef,0,1)
    # tfst::Function = (tfs,tfst_temp,t1,t2,tfst!) ->
    # for i in 1:size(tfs[:,1],1)
    #     for j in 1:size(tfs[i,3],1)
    #         t1 = "$(tfs[i,1])"
    #         t2 = [t2; t1]
    #     end
    #     tfst = [tfst_temp; t2]
    # end
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