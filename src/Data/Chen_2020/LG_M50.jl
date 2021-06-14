using Parameters

@with_kw mutable struct Constants
    T::Float64 = 273.15
    T_ref::Float64 = 298.15
        t_plus::Float64 = 0.2594 #Inital Transference Number
        tpf::Function = ce -> -0.1287*ce^3+0.4106*ce^2-0.4717*ce+0.4492 #Transference Number Function
        De::Float64 = 1.0e-11   #Inital Electrolyte Diffusivity
        Def::Function = ce -> 8.794e-11*ce^2-3.972e-10*ce+4.862e-10 #Electrolyte Diffusivity Function
    SOC::Float64 = 1.
    ce0::Float64 = 1000
    dln::Float64 = 3.0
    Ea_De::Float64 = 0.0
        CC_A::Float64 = 0.1027  #Electrode Plate Area 
    κ::Float64 = 1.0
    κf::Function = ce -> 4.1253e-2+500.7*ce*(1e-6)-4.7212e5*ce^2*1e-12+1.5094e8*ce^3*(1e-18)-1.6018e10*ce^4*1e-24
    Uocp::Function = (Electrode, θ) ->
        if Electrode == "Neg"
            Uocp = @. 1.97938*2.7182818284^(-39.3631*θ) + 0.2482 - 0.0909*tanh(29.8538*(θ - 0.1234)) - 0.04478*tanh(14.9159*(θ - 0.2769)) - 0.0205*tanh(30.4444*(θ - 0.6103))
        else
            Uocp = @. -0.8090*θ + 4.4875 - 0.0428*tanh(18.5138*(θ - 0.5542)) - 17.7326*tanh(15.7890*(θ - 0.3117)) + 17.5842*tanh(15.9308*(θ - 0.3120))
        end

    ∂Uocp::Function = (Electrode,θ) -> 
        if Electrode == "Neg"
            ∂Uocp = @. 0.667934002(tanh(14.9159*θ - 4.13021271)^2) + 2.71371042(tanh(29.8538*θ - 3.68395892)^2) + 0.6241102(tanh(30.4444*θ - 18.58021732)^2) - 4.005754622 - (77.91453287630758(2.7182818284^(-39.3631*θ)))
        else
            ∂Uocp = @. 279.9800214(tanh(15.789*θ - 4.9214313)^2) + 0.79239064(tanh(18.5138*θ - 10.26034796)^2) - 1.4510386800000594 - (280.13037336(tanh(15.9308*θ - 4.9704096)^2))
        end
end

@with_kw mutable struct Negative
        L::Float64 = 85.2e-6    #Electrode Length
        Rs::Float64 = 5.86e-6   # Particle radius [m]
        Ds::Float64 = 1.74e-15   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
        ϵ_s::Float64 = 0.75     #Active Material Volume Fraction
        ϵ_e::Float64 = 0.25    # Porosity of negative electrode
        De_brug::Float64 = 1.5
        κ_brug::Float64 = 1.5
        σ::Float64 = 215    #Solid Phase Electronic Conductivity
        σ_brug::Float64 = 1.5
        θ_100::Float64 = 0.9014
        θ_0::Float64 = 0.0279
        cs_max::Float64 = 29583
        α::Float64 = 0.5
        k_norm::Float64 = 6.48e-7
        Ea_κ::Float64 = 35000   #Activation Energy
    RFilm::Float64 = 0.
    D1::Float64 = 1.0   #Init Value
    D1f::Function = De -> De * ϵ_e^De_brug #Effective Diffusivity
    as::Float64 = 3.0*ϵ_s/Rs # Specific interfacial surf. area
end

@with_kw mutable struct Positive
        L::Float64 = 75.6e-6    #Electrode Length
        Rs::Float64 = 5.22e-6    # Particle radius [m]
        Ds::Float64 = 1.48e-15   # Solid diffusivity [m^2/s]
    Ea_σ::Float64 = 0.0
    Ea_Ds::Float64 = 0.0
        ϵ_s::Float64 = 0.665     #Active Material Volume Fraction
        ϵ_e::Float64 = 0.335    # Porosity of positive electrode
        De_brug::Float64 = 1.5
        κ_brug::Float64 = 1.5
        σ::Float64 = 0.18   #Solid Phase Electronic Conductivity
        σ_brug::Float64 = 1.5
        θ_100::Float64 = 0.2661
        θ_0::Float64 = 0.9084
        cs_max::Float64 = 51765
        α::Float64 = 0.5
        k_norm::Float64 = 3.42e-6
        Ea_κ::Float64 = 17800   #Activation Energy
    RFilm::Float64 = 0.
    D3::Float64 = 1.0   #Init Value
    D3f::Function = De -> De * ϵ_e^De_brug
    as::Float64 = 3.0*ϵ_s/Rs #Specific interfacial surf. area
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
    H1::Array{Int64,1} = 0:3000 #4000 #4612
    H2::Array{Int64,1} = 0:3000 #4000 #4612
    Outs::Int64 = 25
end

@with_kw mutable struct TransferFun
    tfs =   [[C_e, Phi_e, C_se, Phi_s, Phi_se, Flux, C_se, Phi_s, Flux, Phi_se] ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] [Number[0, 4.26E-05, 8.52E-05, 9.72E-05, 1.35E-04, 1.73E-04], Number[4.26E-05, 8.52E-05, 9.72E-05, 1.35E-04, 1.73E-04], Number[0,1], Number[1],Number[0,1],Number[0,1],Number[0,1],Number[1],Number[0,1],Number[0,1]]]
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