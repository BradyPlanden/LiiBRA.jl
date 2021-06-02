using Parameters

@with_kw mutable struct Constants
    T::Float64 = 273.15
    T_ref::Float64 = 298.15
    t_plus::Float64 = 0.363
    De::Float64 = 7.5e-11
    SOC::Float64 = 1.
    ce0::Float64 = 2000
    dln::Float64 = 3.0
    Ea_κ::Float64 = 0.0
    Ea_De::Float64 = 0.0
    Lnegsep::Float64 = 0.0  #Init Value
    CC_A::Float64 = 1.0 
    κ::Float64 = 1.0
    κf::Function = ce -> 4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24
    Uocp::Function = (Electrode, θ) ->
        if Electrode == "Neg"
            Uocp = @. -0.16+1.32*exp(-3*θ)+10*exp(-2000*θ)
        else
            Uocp = @. -20000*exp(-2000*θ)-3.96*exp(-3*θ)
        end

    ∂Uocp::Function = (Electrode,θ) -> 
        if Electrode == "Neg"
            ∂Uocp = @. 4.19829+0.0565661*tanh(-14.5546*θ+8.60942)-0.0275479*(1/(0.998432-θ)^0.492465-1.90111)-0.157123*exp(-0.04738*θ^8)+0.810239*exp(-40*(θ-0.133875))
        else
            ∂Uocp = @. -32.4096*exp(-40*(-0.133875+θ))-0.0135664/((0.998432-θ)^1.49247)+0.0595559*exp(-0.04738*θ^8)*θ^7-0.823297*(sech(8.60942-14.5546*θ))^2
        end
end
@with_kw mutable struct Negative
    L::Float64 = 1.28e-4
    Rs::Float64 = 1.25e-5   # Particle radius [m]
    Ds::Float64 = 3.9e-14   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
    ϵ_s::Float64 = 0.471
    ϵ_e::Float64 = 0.357    # Porosity of negative electrode
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
    σ::Float64 = 100
    σ_brug::Float64 = 1.0
    θ_100::Float64 = 0.53
    θ_0::Float64 = 0.05
    cs_max::Float64 = 26390
    α::Float64 = 0.5
    k_norm::Float64 = 2.2842e-05
    RFilm::Float64 = 0.0
    k_ref::Float64 = 2.16e-11
    D1::Float64 = 1.0   #Init Value
    D1f::Function = De -> De * ϵ_e^De_brug #Effective Diffusivity
    as::Float64 = 3*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Positive
    L::Float64 = 1.9e-4
    Rs::Float64 = 8.5e-6    # Particle radius [m]
    Ds::Float64 = 1.0e-13   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
    ϵ_s::Float64 = 0.297
    ϵ_e::Float64 = 0.444    # Porosity of positive electrode
    De_brug::Float64 = 1.5
    κ_brug::Float64 = 1.5
    σ::Float64 = 3.8
    σ_brug::Float64 = 1.0
    θ_100::Float64 = 0.17
    θ_0::Float64 = 0.762
    cs_max::Float64 = 22860
    α::Float64 = 0.5
    k_norm::Float64 = 2.2073e-5
    RFilm::Float64 = 0
    D3::Float64 = 1.0   #Init Value
    D3f::Function = De -> De * ϵ_e^De_brug
    as::Float64 = 3*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Seperator
    L::Float64 = 7.6e-5
    ϵ_e::Float64 = 0.724    # Porosity of separator
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
    H1::Array{Int64,1} = 0:3000 #4000 #4612
    H2::Array{Int64,1} = 0:3000 #4000 #4612
    Outs::Int64 = 21
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

CellData = Cell(Constants(),Negative(),Positive(),Seperator(),RealisationAlgorthim(), TransferFun())