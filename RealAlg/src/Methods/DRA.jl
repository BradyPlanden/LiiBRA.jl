function DRA(CellData::Cell,FCall::FCalls,L::NTuple{10,Array{Number,1}},TransferFuns)
    """ 
    Discrete Realiastion Algorithm
    # Add License
    # Add Ins and Outs
        # Cell Data
        # Frequency Vector
        # Discretisation Locations
        # Electrode Definition
    """

    # Create s Vector
    Ts = 1/CellData.RA.Fs
    Tlen = CellData.RA.Tlen
    Nfft = 2^(ceil(log2(CellData.RA.Fs*Tlen)))
    f = 0:Nfft-1
    s = (2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)

    # Loop Call Transfer Functions ---------------------------------
    tfft = @. Ts*f
    OrgT = @. (CellData.RA.SamplingT*(0:floor(tfft[end])/CellData.RA.SamplingT))
    #println("tfft:",tfft)
    println("OrgT:",size(OrgT))
    i=Int64(1)

    #Initialise Loop Variables
    puls = Array{Float64}(undef,0,size(OrgT,1)-1)
    D = Array{Float64}(undef,0,1)
    C_Aug = Array{Float64}(undef,0,1)
    DC_Gain = Array{Float64}(undef,Tlen,1)

    for run in TransferFuns.tfs[:,1]
        if TransferFuns.tfs[i,2] == "Pos"
            tf, D_term, res0 = run(CellData,FCall,s,TransferFuns.tfs[i,3],"Pos")
        elseif TransferFuns.tfs[i,2] == "Neg"
            tf, D_term, res0 = run(CellData,FCall,s,TransferFuns.tfs[i,3],"Neg")
        else 
            tf, D_term, res0 = run(CellData,FCall,s,TransferFuns.tfs[i,3])
        end
        #println("D:",D[:,1], "\n")
        jk = CellData.RA.Fs.*real(ifft(tf)) # inverse fourier transform tranfser function response
        stpsum = cumsum(jk, dims=1).*Ts # cumulative sum of tf response * sample time
        nR = size(stpsum,1)
        samplingtf = zeros(nR,length(OrgT))
        dsTf = zeros(nR,length(OrgT))
        println("nR:", nR)
        #println("stpsum:", size(stpsum))
        # Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
            for Output in 1:nR
                spl1 = Spline1D(tfft,stpsum[Output,:]; k=3, s=0.0)
                samplingtf[Output,:]= evaluate(spl1,OrgT)
                dsTf[Output,:] = derivative(spl1,samplingtf[Output,:])
            end

        # println("samplingtf",samplingtf)
         println("dsTf:",size(dsTf))
         println("res0:",res0)
        puls = [puls; dsTf[:,2:end]]
        D = [D; D_term]
        C_Aug = [C_Aug; res0]
        DC_Gain = [DC_Gain; tf[:,1]]
        i = i + 1
        #println("puls:",size(puls))
        if Debug == 1
            #println("jk:",jk, "\n")
            println("stpsum:",size(stpsum))
           #println("tf:",tf, "\n")
            println("D:",D)
            println("nR:",nR)
        end
    end
    println("D:",D)
    println("C_Aug:",size(C_Aug))
    println("puls:",puls[:,1])
    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))
    println("SFactor:",SFactor)
    puls = puls./SFactor
    println("puls_scaled:",puls[:,1])

    # Hankel Formation, perform svd to determine the highest order singular values
    Puls_L = size(puls,1)
    println("Puls_L:",Puls_L)
    Hank1 = zeros(length(CellData.RA.H1)*Puls_L,length(CellData.RA.H2))
    Hank2 = Hank1
    println("Hank1:",size(Hank1))
    println("length H2:",length(CellData.RA.H2))
    for lp1 in 1:length(CellData.RA.H2)
        for lp2 in 1:length(CellData.RA.H1)
            #println("H1:",length(CellData.RA.H1))
            #println("lp1:",lp1)
            #println("lp2:",lp2)
            Hank1[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = puls[:,CellData.RA.H2[lp1]+CellData.RA.H1[lp2]+1]
            Hank1[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = puls[:,CellData.RA.H2[lp1]+CellData.RA.H1[lp2]+2]
        end
    end

    F = svds(Hank1; nsv=CellData.RA.M)[1]
    println("S:",F.S)
    println("F.U:",size(F.U))

    # Create Observibility and Control Matrices -> Create A, B, and C 
    for OrdArg in 1:1#CellData.RA.M
        S_ = sqrt(diagm(F.S))
        println("S_",S_)
        Observibility = F.U*S_
        Control = S_*F.V'
        A = Observibility\Hank2/Control
        eigA = eigen(A)
        println("A:",A)
        println("eigenA:",eigA.values)
        B = Control[:,1:CellData.RA.N]
        C = Observibility[1:CellData.RA.M,:]
        println("C:",size(C))
        # Transform A,B,C matrices to final form
        C = C*SFactor[ones(Int64,size(C,2)),:]
        if C_Aug == 0
            A_Final = A
            B_Final = B
            C_Final = C
        else
            A_Final = diagm([1;diag(A)])
            B_Final = [CellData.RA.SamplingT; B]
            C_Final = [C_Aug; C]
        end
        println("A_Final:",A_Final)
    end

   
    # Create D matrice from H -> ∞

end