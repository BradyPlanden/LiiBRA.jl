@inline function DRA(CellData::Cell,FCall::FCalls,L::NTuple{10,Array{Number,1}},TransferFuns)
    """ 
    Discrete Realiastion Algorithm
    # Add License
    # Add Ins and Outs
        # Cell Data
        # Frequency Vector
        # Discretisation Locations
        # Electrode Definition
    """

    DRA_Debug = 1

    # Create s Vector
    Ts = 1/CellData.RA.Fs
    Tlen = CellData.RA.Tlen
    Nfft = 2^(ceil(log2(CellData.RA.Fs*Tlen)))
    f = 0:Nfft-1
    s = (2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)
    s=s'
    #  println("s:")
    #  display("text/plain", s[:,65531:end])

    # Loop Call Transfer Functions ---------------------------------
    tfft = @. Ts*f
    OrgT = @. CellData.RA.SamplingT*(0:floor(tfft[end]/CellData.RA.SamplingT))
    #println("tfft:",tfft)
    #println("OrgT:",size(OrgT))
    i=Int64(1)

    #Initialise Loop Variables
    puls = Array{Float64}(undef,0,size(OrgT,1)-1)
    stpsum__ = Array{Float64}(undef,0,length(s))
    tf__ = Array{Float64}(undef,0,length(s))
    jk__ = Array{Float64}(undef,length(s),0)
    #testifft__ = Array{Float64}(undef,length(s),0)
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
        #testifft = ifft(tf,2)'
        jk = CellData.RA.Fs.*real(ifft(tf,2)') # inverse fourier transform tranfser function response
        stpsum = (cumsum(jk, dims=1).*Ts)' # cumulative sum of tf response * sample time
        nR = size(stpsum,1)
        samplingtf = zeros(nR,length(OrgT))
        dsTf = zeros(nR,length(OrgT))
        #println("nR:", nR)

        #println("tf:")
        #display("text/plain", tf[:,1:10])
        #println("size jk:",size(jk))
        #println("size stpsum:",size(stpsum))
        # Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
            for Output in 1:nR
                spl1 = Spline1D(tfft,stpsum[Output,:]; k=3)
                samplingtf[Output,:]= evaluate(spl1,OrgT)
                #dsTf[Output,:] = diff(samplingtf[Output,:], dims=2)
            end
            dsTf = [zeros(size(samplingtf,1)) diff(samplingtf, dims=2)]

        # println("samplingtf:")
        # display("text/plain", samplingtf[:,1:10])
        #  println("dsTf:",size(dsTf))
        #  println("stpsum:",size(stpsum))
        puls = [puls; dsTf[:,2:end]]
        D = [D; D_term]
        C_Aug = [C_Aug; res0]
        DC_Gain = [DC_Gain; tf[:,1]]
        i = i + 1
        if Debug == 1
            println("jk:",jk, "\n")
            println("stpsum:",size(stpsum))
            println("tf:",tf, "\n")
            println("D:",D)
            println("nR:",nR)
        end

        if DRA_Debug == 1
            stpsum__ = [stpsum__; stpsum]
            tf__ = [tf__; tf]
            jk__ = [jk__ jk]
            #testifft__ = [testifft__ testifft]
        end
    end
    puls = reverse!(puls, dims=2)
    if DRA_Debug == 1
        # println("tf__:",size(tf__))
        # display("text/plain", tf__)

        #println("testifft__:",size(testifft__))
        #display("text/plain", testifft__)

        # println("jk__:",size(jk__))
        # display("text/plain", jk__)

        # println("stpsum__:")
        # display("text/plain", stpsum__)
    end
    # println("pulsns:")
    # display("text/plain", puls)


    # println("D:")
    # display("text/plain", D)

    # println("C_Aug:",size(C_Aug))
    # println("puls:",puls[:,1])
    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))

     println("SFactor:")
     display("text/plain", SFactor)

    puls = puls./SFactor

     println("puls:")
     display("text/plain", puls)

    # Hankel Formation, perform svd to determine the highest order singular values
    Puls_L = size(puls,1)
    # println("Puls_L:",Puls_L)
    Hank1 = zeros(length(CellData.RA.H1)*Puls_L,length(CellData.RA.H2))
    Hank2 = zeros(length(CellData.RA.H1)*Puls_L,length(CellData.RA.H2))
    # println("Hank1:",size(Hank1))
    # println("length H2:",length(CellData.RA.H2))
    for lp1 in 1:length(CellData.RA.H2)
        for lp2 in 1:length(CellData.RA.H1)
            #println("H1:",length(CellData.RA.H1))
            #println("lp1:",lp1)
            #println("lp2:",lp2)
            Hank1[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = puls[:,CellData.RA.H2[lp1]+CellData.RA.H1[lp2]+1]
            Hank2[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] = puls[:,CellData.RA.H2[lp1]+CellData.RA.H1[lp2]+3]
        end
    end
    #println("test",CellData.RA.H2[100]+CellData.RA.H1[2]+2)

    println("Hank1:")
    display("text/plain", Hank1)

    println("Hank2:")
    display("text/plain", Hank2)

    F = svds(Hank1; nsv=CellData.RA.M)[1]

    println("F.U")
    display("text/plain", F.U)

    println("F.V")
    display("text/plain", F.V)

    println("F.S")
    display("text/plain", F.S)

    # Create Observibility and Control Matrices -> Create A, B, and C 
    # for OrdArg in 1:CellData.RA.M
        S_ = sqrt(diagm(F.S))
         println("S_")
         display("text/plain", S_)
        Observibility = F.U[:,1:CellData.RA.M]*S_

         println("Observibility:")
         display("text/plain", Observibility)

        Control = S_*F.V[:,1:CellData.RA.M]'

         println("Control:")
         display("text/plain", Control)

        A = Observibility\Hank2/Control
        eigA = eigen(A)

        println("eigA:")
        display("text/plain", eigA)

        # println("eigenA:",eigA.values)
        B = Control[:,1:CellData.RA.N]
        C = Observibility[1:size(puls,1),:]

        println("B:")
        display("text/plain", B)

        println("C:")
        display("text/plain", C)

        # Transform A,B,C matrices to final form
        C = C.*SFactor[ones(Int64,size(C,1)),:]
        if C_Aug == 0
            A_Final = A
            B_Final = B
            C_Final = C
        else
            A_Final = diagm([1;eigA.values])
            B_Final = [CellData.RA.SamplingT; B]
            C_Final = [C_Aug C]
        end
        println("C_Final:")
        display("text/plain", C_Final)
  
        println("A_Final:")
        display("text/plain", A_Final)
    #end

    # Create D matrice from H -> ∞
return eigA
end