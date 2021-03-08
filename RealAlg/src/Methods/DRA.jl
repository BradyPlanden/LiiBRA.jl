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
#Tlen = 32768
Tlen = 256
Nfft = 2^(ceil(log2(CellData.RA.Fs*Tlen)))
f = 0:Nfft-1
s = (2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)

# Loop Call Transfer Functions ---------------------------------
tfft = @. Ts*f
OrgT = @. (CellData.RA.SamplingT*(0:floor(tfft[end])/CellData.RA.SamplingT))
println("tfft:",tfft)
println("OrgT:",OrgT)
i=Int64(1)
puls = []
D = []
C_Aug = []
DC_Gain = []
 for run in TransferFuns.tfs[:,1]
    if TransferFuns.tfs[i,2] == "Pos"
       tf, D, res0 = run(CellData,FCall,s,TransferFuns.tfs[i,3],"Pos")
    elseif TransferFuns.tfs[i,2] == "Neg"
        tf, D, res0 = run(CellData,FCall,s,TransferFuns.tfs[i,3],"Neg")
    else 
        tf, D, res0 = run(CellData,FCall,s,TransferFuns.tfs[i,3])
    end
    jk = CellData.RA.Fs.*real(ifft(permutedims(tf))) # inverse fourier transform tranfser function response
    stpsum = (cumsum(jk, dims=1).*Ts)' # cumulative sum of tf response * sample time
    nR = size(stpsum,1)
    samplingtf = zeros(nR,length(OrgT))
    dsTf = zeros(nR,length(OrgT))

    println("stpsum:",size(stpsum))
    println("samplingtf:",size(samplingtf))
    for Output in 1:length(nR)
        spl1 = Spline1D(tfft,stpsum[Output,:]; k=3, s=0.0)
        samplingtf[Output,:]= evaluate(spl1,OrgT)
        dsTf = derivative(spl1,samplingtf[Output,:])
    end
    puls = [puls, dsTf[:,2:end]]
    D = [D, D]
    C_Aug = [C_Aug, res0]
   # DC_Gain = [DC_Gain, tf(:,1)]

    i = i + 1
    if Debug == 1
        #println("jk:",jk, "\n")
        println("stpsum:",size(stpsum))
        #println("tf:",tf, "\n")
        println("D:",D)
        println("nR:",nR)
    end
end


# Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
# Create Hankel Matrix, perform svd to determine the highest order singular values
# Create Observibility and Control Matrices -> Create A, B, and C 
# Transform A,B,C matrices to final form
# Create D matrice from H -> ∞

end