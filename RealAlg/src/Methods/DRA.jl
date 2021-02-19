function DRA(CellData::Cell,L::Tuple{Array{Any,1},Array{Any,1},Array{Any,1},Int64})
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
M=4
Fs = 1/CellData.RA.Ts
Tlen = 32768
Nfft = 2^(ceil(log2(Fs*Tlen)))
f = 0:Nfft-1
s = (2im.*Fs)*tan.(pi.*f./Nfft)
println("s:",size(s))
#z = Any[0, 128e-6, 204e-6, 394e-6]

#Transfer functions to call
tfs = [
    C_e(CellData,s,L[1],M),
    C_se(CellData,s,L[3],"Pos"),
    C_se(CellData,s,L[3],"Neg"),
    Phi_e(CellData,s,L[2]),
    Phi_s(CellData,s,L[4],"Pos"),
    Phi_s(CellData,s,L[4],"Neg"),
    Phi_se(CellData,s,L[3],"Pos"),
    Phi_se(CellData,s,L[3],"Neg"),
    j(CellData,s,L[3],"Pos"),
    j(CellData,s,L[3],"Neg")
]

# Call Transfer Functions
Numtf = 9 # Number of Transfer functions
tfft = @. CellData.RA.Ts*f
i = 1
for run in tfs
      tf = fill(0.0,size(L[i]))
      tf = fill(tf,length(s))
      tf = run
      i = i+1
      println("i:",i)
end



# Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
# Create Hankel Matrix, perform svd to determine the highest order singular values
# Create Observibility and Control Matrices -> Create A, B, and C 
# Transform A,B,C matrices to final form
# Create D matrice from H -> ∞

end