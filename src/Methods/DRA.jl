function DRA(CellData,s,f)
    """ 
    Function for Discrete Realisation Algorithm.

    DRA(CellData,s,f)
    
    """

    # Loop Call Transfer Functions ---------------------------------
    tfft = @. (1/CellData.RA.Fs)*f
    OrgT = @. CellData.RA.SamplingT*(0:floor(tfft[end]/CellData.RA.SamplingT))
    i=Int64(1)

    #Initialise Loop Variables
    puls = Array{Float64}(undef,0,size(OrgT,1)-1)
    D = Array{Float64}(undef,0,1)
    Dtt = Array{String}(undef,0,1)
    C_Aug = Array{Float64}(undef,0,1)
    #tf = Vector{Array}
    # stpsum__ = Array{Float64}(undef,0,length(s))
    tf__ = Array{Float64}(undef,0,length(s))
    # jk__ = Array{Float64}(undef,length(s),0)
    # testifft__ = Array{Float64}(undef,length(s),0)
    #D_term = Array{Float64}(undef,0,1)
    #res0 = Array{Float64}(undef,0,1)

    for run in CellData.Transfer.tfs[:,1]
        if CellData.Transfer.tfs[i,2] == "Pos"
            tf, Di, res0, Dti = run(CellData,s,CellData.Transfer.tfs[i,3],"Pos") #high compute line
        elseif CellData.Transfer.tfs[i,2] == "Neg"
            tf, Di, res0, Dti = run(CellData,s,CellData.Transfer.tfs[i,3],"Neg") #high compute line
        else 
            tf, Di, res0, Dti = run(CellData,s,CellData.Transfer.tfs[i,3]) #high compute line
        end

        #tf = reverse!(tf, dims=2)
        #c = plan_ifft(tf,2)
        #jk = CellData.RA.Fs*real((c*tf)') # inverse fourier transform tranfser function response
        jk = CellData.RA.Fs*real(ifft(tf,2)') # inverse fourier transform tranfser function response (Large Compute)
        stpsum = (cumsum(jk, dims=1).*(1/CellData.RA.Fs))' # cumulative sum of tf response * sample time
        samplingtf = Array{Float64}(undef,size(stpsum,1),length(OrgT))
        # Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
        for Output in 1:size(stpsum,1)
            spl1 = Spline1D(tfft,stpsum[Output,:]; k=3) #Large compute - First to improve
            samplingtf[Output,:]= evaluate(spl1,OrgT)
        end
        puls = [puls; diff(samplingtf, dims=2)]
        D = [D; Di]
        Dtt = [Dtt; Dti]
        C_Aug = [C_Aug; res0]
        tf__ = [tf__;tf]
        i += 1

    end

    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))
    puls1 = puls
    puls = puls./SFactor[:,ones(Int64,size(puls,2))]

    #Hankel Formation, perform svd to determine the highest order singular values
    Puls_L = size(puls,1)
    Hank1 = Array{Float64}(undef,length(CellData.RA.H1)*Puls_L,length(CellData.RA.H2))
    Hank2 = Array{Float64}(undef,length(CellData.RA.H1)*Puls_L,length(CellData.RA.H2))

    for lp1 in 1:length(CellData.RA.H2)
        for lp2 in 1:length(CellData.RA.H1)
             Hank1[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = @view puls[:,CellData.RA.H2[lp1]+CellData.RA.H1[lp2]+1] #High compute line
             Hank2[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = @view puls[:,CellData.RA.H2[lp1]+CellData.RA.H1[lp2]+2] #High compute line
        end
    end

    #Truncated SVD of Hank1 Matrix
    T = svds(Hank1; nsv=CellData.RA.M)[1]
    #U,S,V = tsvd(Hank1, k=CellData.RA.M)

    # Create Observibility and Control Matrices -> Create A, B, and C 
    S_ = sqrt(diagm(T.S[1:CellData.RA.M]))
    Observibility = (@view T.U[:,1:CellData.RA.M])*S_
    Control = S_*(@view T.V[:,1:CellData.RA.M])'
    
    A = Matrix{Float64}(I,CellData.RA.M+1,CellData.RA.M+1)
    A[2:end,2:end] = (Observibility\Hank2)/Control #High compute line (Second)

    #Error check
    if any(i -> i>1., real(eigvals(A)))
        println("Oscilating System: A has indices of values > 1")
    end

    if any(i -> i<0., real(eigvals(A)))
        println("Unstable System: A has indices of negative values")
    end
    
    B = [CellData.RA.SamplingT; Control[:,1:CellData.RA.N]]
    C = [C_Aug SFactor[:,ones(Int64,CellData.RA.M)].*Observibility[1:size(puls,1),:]]

    #@infiltrate

    #  Final State-Space Form
    d, S = eigen(A)
    A = Diagonal(d)
    B = S\B
    C = C*S

    # sys = diagonalize(ss(A,B,C,D,CellData.RA.SamplingT))
    # ss_B = sys.B
    # ss_C = sys.C

    #Scale C and tansform B to improve real-time linearisation performance
    # ss_C = ss_C.*ss_B'
    # ss_B = ones(size(ss_B)) 

    C = C.*B'
    B = ones(size(B))
        
return real(A), real(B), real(C), D
end