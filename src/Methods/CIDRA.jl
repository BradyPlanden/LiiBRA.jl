@inline function CIDRA(Cell)
    """ 
    Function for Discrete Realisation Algorithm.

    CIDRA(Cell,s,f)
    
    """

    # Method Check
    if Cell.RA.Fs==(1/Cell.RA.SamplingT)
        println("CIDRA Utilised")
    else
        println("DRA Utilised")
    end

    #Additional Pulse Setup
    tfft = (1/Cell.RA.Fs)*Cell.RA.f
    OrgT = Cell.RA.SamplingT*(0:floor(tfft[end]/Cell.RA.SamplingT))

    #Initialise Loop Variables
    A = Matrix{Float64}(I,Cell.RA.M+1,Cell.RA.M+1)
    B = Vector{Float64}(undef,Cell.RA.M+1)
    C = Matrix{Float64}(undef,Cell.RA.Outs,Cell.RA.M+1)
    D = Vector{Float64}(undef,Cell.RA.Outs)
    puls = Array{Float64}(undef,Cell.RA.Outs,size(OrgT,1)-1)
    C_Aug = Vector{Float64}(undef,Cell.RA.Outs)
    i=Int(1)
    l=Int(1)
    u=Int(0)

    for run in Cell.Transfer.tfs
        tf = Array{ComplexF64}(undef,size(Cell.Transfer.Locs[i],1),size(Cell.RA.s,2))
        Di = Vector{Float64}(undef,size(Cell.Transfer.Locs[i],1))::Vector{Float64}
        res0 = Vector{Float64}(undef,size(Cell.Transfer.Locs[i],1))::Vector{Float64}
        samplingtf = Array{Float64}(undef,size(Cell.Transfer.Locs[i],1),length(OrgT))
        
        if Cell.Transfer.Elec[i] == "Pos"
            run(Cell,Cell.RA.s,Cell.Transfer.Locs[i],"Pos",tf,Di,res0) 
        elseif Cell.Transfer.Elec[i] == "Neg"
            run(Cell,Cell.RA.s,Cell.Transfer.Locs[i],"Neg",tf,Di,res0) 
        else 
            run(Cell,Cell.RA.s,Cell.Transfer.Locs[i],tf,Di,res0) 
        end
        
        # if Cell.RA.Fs==(1/Cell.RA.SamplingT)
        #     f = one(Int)
        #     λ = Int(Cell.RA.Fs*Cell.RA.SamplingT)
        #     for i in 1:λ:size(tfft,2), j in 1:λ:size(tfft,1)
        #         Org[j,f] = tf[j,i]
        #         f+=1
        #     end
        # end

        u += Int(size(Cell.Transfer.Locs[i],1))::Int

        if Cell.RA.Fs==(1/Cell.RA.SamplingT)
            puls[l:u,:] .= real(ifft(tf,2))[:,2:end]
        else
            stpsum = cumsum(real(ifft(tf,2)), dims=2) #cumulative sum of tf response * sample time
            #Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
            for k in 1:size(stpsum,1)
                spl1 = CubicSplineInterpolation(tfft, stpsum[k,:]; bc=Line(OnGrid()), extrapolation_bc=Throw())
                samplingtf[k,:] = spl1(OrgT)
            end
            puls[l:u,:] .= diff(samplingtf, dims=2)
        end

        D[l:u,:] .= Di
        C_Aug[l:u,:] .= res0
        l = u+one(u)::Int
        i += one(i)

    end

    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))
    puls .= puls./SFactor

    #Pre-Allocation for Hankel & SVD
    Puls_L = size(puls,1)
    𝐇 = Array{Float64}(undef,length(Cell.RA.H1)*Puls_L,length(Cell.RA.H2))
    U,S,V = fh!(𝐇,Cell.RA.H1,Cell.RA.H2,puls,Cell.RA.M,Puls_L)

    # Create Observibility and Control Matrices -> Create A, B, and C 
    ʃ = sqrt(Diagonal(S))
    U .= U*ʃ
    V .= ʃ*V
    
    A[2:end,2:end] .= (U\𝐇)/V
    B .= [Cell.RA.SamplingT; V[:,1:Cell.RA.N]]
    C .= [C_Aug SFactor.*U[1:Puls_L,:]]

    #Performance check
    if any(i -> i>1., real(eigvals(A)))
        println("Oscilating System: A has indices of values > 1")
    end
    
    if any(i -> i<0., real(eigvals(A)))
        println("Unstable System: A has indices of negative values")
    end

    return mag!(A), mag!(B), mag!(C), D
 
end