function DRA(CellData::Cell)
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
    C_e(CellData,s,Any[0, 128e-6, 204e-6, 394e-6],M),
    C_se(CellData,s,Any[0,1],"Pos"),
    C_se(CellData,s,Any[0,1],"Neg"),
    Phi_e(CellData,s,Any[128e-6, 204e-6, 394e-6]),
    Phi_s(CellData,s,1,"Pos"),
    Phi_s(CellData,s,1,"Neg"),
    Phi_se(CellData,s,Any[0,1],"Pos"),
    Phi_se(CellData,s,Any[0,1],"Neg"),
    Flux(CellData,s,Any[0,1],"Pos"),
    Flux(CellData,s,Any[0,1],"Neg")
]

# Call Transfer Functions
Numtf = 9 # Number of Transfer functions
tfft = @. Ts*(f)
i = 1
for run in tfs
      tf[i] = run
      i = i+1

end



# Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
# Create Hankel Matrix, perform svd to determine the highest order singular values
# Create Observibility and Control Matrices -> Create A, B, and C 
# Transform A,B,C matrices to final form
# Create D matrice from H -> ∞

end