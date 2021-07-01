@inline function DRA(CellData,s,f)
    """ 
    Discrete Realiastion Algorithm
    # Add License
    # Add Ins and Outs
        # Cell Data
        # Frequency Vector
        # Discretisation Locations
        # Electrode Definition
    """

    DRA_Debug = 0

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
    # tf__ = Array{Float64}(undef,0,length(s))
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

        jk = CellData.RA.Fs.*real(ifft(tf,2)') # inverse fourier transform tranfser function response
        stpsum = (cumsum(jk, dims=1).*(1/CellData.RA.Fs))' # cumulative sum of tf response * sample time
        nR = size(stpsum,1)
        samplingtf = Array{Float64}(undef,nR,length(OrgT))

        # Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
            for Output in 1:nR
                spl1 = Spline1D(tfft,stpsum[Output,:]; k=3) #High compute line
                samplingtf[Output,:]= evaluate(spl1,OrgT)
            end

            dsTf = [Array{Float64}(undef,size(samplingtf,1)) diff(samplingtf, dims=2)]
            puls = [puls; dsTf[:,2:end]]
            D = [D; Di]
            Dtt = [Dtt; Dti]
            C_Aug = [C_Aug; res0]
            i = i+1

        if Debug == 1
            println("jk:",jk, "\n")
            println("stpsum:",size(stpsum))
            println("tf:",tf, "\n")
            println("D:",D)
            println("nR:",nR)
        end

        if DRA_Debug == 1
            # stpsum__ = [stpsum__; stpsum]
            # tf__ = [tf__; tf]
            # jk__ = [jk__ jk]
            # testifft__ = [testifft__ testifft]
        end
    end

    puls = reverse!(puls, dims=2)

    if DRA_Debug == 1
        println("tf__:",size(tf__))
        display("text/plain", tf__)

        println("testifft__:",size(testifft__))
        display("text/plain", testifft__)

        println("jk__:",size(jk__))
        display("text/plain", jk__)

        println("stpsum__:")
        display("text/plain", stpsum__)
    end

    #Scale Transfer Functions in Pulse Response
    SFactor = sqrt.(sum(puls.^2,dims=2))
    puls .= puls./SFactor
   # println("puls:",size(puls))
   
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
    F = svds(Hank1; nsv=CellData.RA.M)[1]

    # Create Observibility and Control Matrices -> Create A, B, and C 
        S_ = sqrt(diagm(F.S))
        Observibility = (@view F.U[:,1:CellData.RA.M])*S_
        Control = S_*(@view F.V[:,1:CellData.RA.M])'
        A = Observibility\Hank2/Control #High compute line

        eigA = eigvals(A)
        E = diagm([1;eigA])
        Ei = E'

        B = @view Control[:,1:CellData.RA.N]
        C = @view Observibility[1:size(puls,1),:]

        # Transform A,B,C matrices to final form
        C = C.*SFactor[:,ones(Int64,size(C,2))]
        if C_Aug == 0
            A_Final = A
            B_Final = diagm([eigA])'*B
            C_Final = C*diagm([eigA])
        else
            A_Final = E
            B_Final = Ei*[CellData.RA.SamplingT; B]
            C_Final = [C_Aug C]*E
        end
        #Scale C and tansform B to improve real-time linearisation performance
        C_Final = C_Final.*B_Final'
        B_Final = ones(size(B_Final))
        
        #  println("A_Final:")
        #  display("text/plain", A_Final)

        #  println("B_Final:")
        #  display("text/plain", B_Final)

        #  println("C_Final:")
        #  display("text/plain", C_Final)

        #  println("D_Final:")
        #  display("text/plain", D)
    #end

return A_Final, B_Final, C_Final, D, Dtt
end