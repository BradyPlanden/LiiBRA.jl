@inline function DRA(Cell,s,f)
    """ 
    Function for Discrete Realisation Algorithm.

    DRA(Cell,s,f)
    
    """

    #Additional Pulse Setup
    tfft = @. (1/Cell.RA.Fs)*f
    OrgT = @. Cell.RA.SamplingT*(0:floor(tfft[end]/Cell.RA.SamplingT))


    #Initialise Loop Variables
    puls = Array{Float64}(undef,Cell.RA.Outs,size(OrgT,1)-1)
    D = Vector{Float64}(undef,Cell.RA.Outs)
    C_Aug = Vector{Float64}(undef,Cell.RA.Outs)
    i=Int64(1)
    l=Int64(1)
    u=Int64(0)


    for run in Cell.Transfer.tfs[:,1]
        tf = Array{ComplexF64}(undef,size(Cell.Transfer.tfs[i,3],1),size(s,2))
        Di = Vector{Float64}(undef,size(Cell.Transfer.tfs[i,3],1))
        res0 = Vector{Float64}(undef,size(Cell.Transfer.tfs[i,3],1))
        samplingtf = Array{Float64}(undef,size(Cell.Transfer.tfs[i,3],1),length(OrgT))

        if Cell.Transfer.tfs[i,2] == "Pos"
            run(Cell,s,Cell.Transfer.tfs[i,3],"Pos",tf,Di,res0) 
        elseif Cell.Transfer.tfs[i,2] == "Neg"
            run(Cell,s,Cell.Transfer.tfs[i,3],"Neg",tf,Di,res0) 
        else 
            run(Cell,s,Cell.Transfer.tfs[i,3],tf,Di,res0) 
        end
        stpsum = (cumsum(Cell.RA.Fs*real(ifft(tf,2)), dims=2).*(1/Cell.RA.Fs)) #cumulative sum of tf response * sample time

        #Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
        for k in 1:size(stpsum,1)
            spl1 = Spline1D(tfft,stpsum[k,:]; k=3)
            samplingtf[k,:] .= evaluate(spl1,OrgT)
        end
        u += size(Cell.Transfer.tfs[i,3],1)
        puls[l:u,:] .= diff(samplingtf, dims=2)
        D[l:u,:] .= Di
        C_Aug[l:u,:] .= res0
        l = u+1
        i += 1
    end

    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))
    puls .= puls./SFactor

    #Pre-Allocation for Hankel & SVD
    Puls_L = size(puls,1)
    Hank = Array{Float64}(undef,length(Cell.RA.H1)*Puls_L,length(Cell.RA.H2))
    U,S,V = fh!(Hank,Cell.RA.H1,Cell.RA.H2,puls,Cell.RA.M,Puls_L)

    # Create Observibility and Control Matrices -> Create A, B, and C 
    S_ = sqrt(diagm(S[1:Cell.RA.M]))
    Observibility = (@view U[:,1:Cell.RA.M])*S_
    Control = S_*(@view V[:,1:Cell.RA.M])'
    
    A = Matrix{Float64}(I,Cell.RA.M+1,Cell.RA.M+1)
    A[2:end,2:end] .= (Observibility\Hank)/Control

    #Performance check
    if any(i -> i>1., real(eigvals(A)))
        println("Oscilating System: A has indices of values > 1")
    end

    if any(i -> i<0., real(eigvals(A)))
        println("Unstable System: A has indices of negative values")
    end
    
    B = [Cell.RA.SamplingT; Control[:,1:Cell.RA.N]]
    C = [C_Aug SFactor.*Observibility[1:size(puls,1),:]]

    #  Final State-Space Form
    #d, Sᵘ = eigen(A,sortby=nothing)
    #d = mag!(d) # taking the magnitude of S and maintaining the sign of the real values
    #S = mag!(S)
    #A = Diagonal(d)
    #C = C*Sᵘ.*(inv(Sᵘ)*B)'
    #B = ones(size(B))

return mag!(A), mag!(B), mag!(C), D   
 
end