@inline function DRA(Cell,s,f)
    """ 
    Function for Discrete Realisation Algorithm.

    DRA(Cell,s,f)
    
    """

    #Additional Pulse Setup
    tfft = @. (1/Cell.RA.Fs)*f
    OrgT = @. Cell.RA.SamplingT*(0:floor(tfft[end]/Cell.RA.SamplingT))
    i=Int64(1)

    #Initialise Loop Variables
    puls = Array{Float64}(undef,0,size(OrgT,1)-1)
    D = Array{Float64}(undef,0,1)
    Dtt = Array{String}(undef,0,1)
    C_Aug = Array{Float64}(undef,0,1)


    for run in Cell.Transfer.tfs[:,1]
        if Cell.Transfer.tfs[i,2] == "Pos"
            tf, Di, res0, Dti = run(Cell,s,Cell.Transfer.tfs[i,3],"Pos") 
        elseif Cell.Transfer.tfs[i,2] == "Neg"
            tf, Di, res0, Dti = run(Cell,s,Cell.Transfer.tfs[i,3],"Neg") 
        else 
            tf, Di, res0, Dti = run(Cell,s,Cell.Transfer.tfs[i,3]) 
        end

        jk = Cell.RA.Fs*real(ifft(tf,2)') # inverse fourier transform tranfser function response (Large Compute)
        stpsum = (cumsum(jk, dims=1).*(1/Cell.RA.Fs))' # cumulative sum of tf response * sample time
        samplingtf = Array{Float64}(undef,size(stpsum,1),length(OrgT))

        # Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
        for Output in 1:size(stpsum,1)
            spl1 = Spline1D(tfft,stpsum[Output,:]; k=3)
            samplingtf[Output,:]= evaluate(spl1,OrgT)
        end

        puls = [puls; diff(samplingtf, dims=2)]
        D = [D; Di]
        Dtt = [Dtt; Dti]
        C_Aug = [C_Aug; res0]
        i += 1

    end

    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))
    puls1 = puls
    puls = puls./SFactor[:,ones(Int64,size(puls,2))]

    #Hankel Formation
    Puls_L = size(puls,1)
    Hank1 = Array{Float64}(undef,length(Cell.RA.H1)*Puls_L,length(Cell.RA.H2))
    Hank2 = Array{Float64}(undef,length(Cell.RA.H1)*Puls_L,length(Cell.RA.H2))

    for lp1 in 1:length(Cell.RA.H2)
        for lp2 in 1:length(Cell.RA.H1)
             Hank1[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = @view puls[:,Cell.RA.H2[lp1]+Cell.RA.H1[lp2]+1]
             Hank2[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = @view puls[:,Cell.RA.H2[lp1]+Cell.RA.H1[lp2]+2]
        end
    end

    #Truncated SVD of Hank1 Matrix
    #T = svds(Hank1; nsv=Cell.RA.M)[1]
    U,S,V = tsvd(Hank1, k=Cell.RA.M)

    # Create Observibility and Control Matrices -> Create A, B, and C 
    S_ = sqrt(diagm(S[1:Cell.RA.M]))
    Observibility = (@view U[:,1:Cell.RA.M])*S_
    Control = S_*(@view V[:,1:Cell.RA.M])'
    
    A = Matrix{Float64}(I,Cell.RA.M+1,Cell.RA.M+1)
    A[2:end,2:end] = (Observibility\Hank2)/Control

    #Error check
    if any(i -> i>1., real(eigvals(A)))
        println("Oscilating System: A has indices of values > 1")
    end

    if any(i -> i<0., real(eigvals(A)))
        println("Unstable System: A has indices of negative values")
    end
    
    B = [Cell.RA.SamplingT; Control[:,1:Cell.RA.N]]
    C = [C_Aug SFactor[:,ones(Int64,Cell.RA.M)].*Observibility[1:size(puls,1),:]]

    #  Final State-Space Form
    d, S = eigen(A)
    A = Diagonal(d)
    B = S\B
    C = (C*S).*B'
    B = ones(size(B))
        
return real(A), real(B), real(C), D
end