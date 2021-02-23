function DRA(CellData::Cell,L::NTuple{10,Array{Number,1}})
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
Fs = 1/CellData.RA.Ts
Tlen = 128
Nfft = 2^(ceil(log2(Fs*Tlen)))
f = 0:Nfft-1
s = (2im.*Fs)*tan.(pi.*f./Nfft)

# genvar() = (CellData,s,L)
# extract(par, ::typeof(C_e)) = (CellData,s,L[1])
# extract(par, ::typeof(C_se)) = (CellData,s,L[3],)
# extract(par, ::typeof(Phi_e)) = (CellData,s,L[2])
# extract(par, ::typeof(Phi_s)) = (CellData,s,L[4])
# extract(par, ::typeof(Phi_se)) = (CellData,s,L[3])
# extract(par, ::typeof(j)) = (CellData,s,L[3])
#genvar() = (CellData=CellData,s=s,L=L,Def="Neg")

#println("genvar:",typeof(genvar()))
#Transfer functions to call
tfs = [[C_e, C_se, C_se, Phi_e, Phi_s, Phi_s, Phi_se, Phi_se, j, j] ["Na", "Pos", "Neg", "Na", "Pos", "Neg", "Pos", "Neg", "Pos", "Neg"] [Number[0, 128e-6, 204e-6, 394e-6],Number[0,1], Number[0,1], Number[128e-6, 204e-6, 394e-6],Number[1],Number[1],Number[0,1],Number[0,1],Number[0,1],Number[0,1]]]

#k, d = C_e(CellData,s,L[1])
#println("d",d)
#Single calling of each tf test ------------------------------
#= tf, D_term = [
    C_e(CellData,s,L[1]),
    C_se(CellData,s,L[3],"Pos"),
    C_se(CellData,s,L[3],"Neg"),
    Phi_e(CellData,s,L[2]),
    Phi_s(CellData,s,L[4],"Pos"),
    Phi_s(CellData,s,L[4],"Neg"),
    Phi_se(CellData,s,L[3],"Pos"),
    Phi_se(CellData,s,L[3],"Neg"),
    j(CellData,s,L[3],"Pos"),
    j(CellData,s,L[3],"Neg"),
]
println("tf:",size(tf))
println("D_term:",size(D_term)) =#
# Loop Call Transfer Functions ---------------------------------
tfft = CellData.RA.Ts*f
i=Int64(1)
 for run in tfs[:,1]
    if tfs[i,2] == "Pos"
       #var = genvar()
       #par = (extract(var, run), "Pos")
       tf, D = run(CellData,s,tfs[i,3],"Pos")
    elseif tfs[i,2] == "Neg"
        #var = genvar()
        #par = extract(var, run)
        tf, D = run(CellData,s,tfs[i,3],"Neg")
    else 
        tf, D = run(CellData,s,tfs[i,3])
    end
    println("D:",size(D))
    println("tf:",size(tf))
    jk = Fs.*real(ifft(tf)) # inverse fourier transform tranfser function response
    stpsum = cumsum(jk, dims=1).*CellData.RA.Ts # cumulative sum of tf response * sample time
    nR = size(stpsum,2)

    i = i + 1
    #println("jk:",size(jk))
    #println("stpsum:",size(stpsum))
    #println("tf:",typeof(tf))
end

#= #Transfer functions to call
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
i = Number(1)
for run in tfs
      #tf = fill(0.0,size(L[i]))
      #tf = fill(tf,length(s))
      tf = run(CellData,s,L[i])
      i = i+1
      #println("tfs:",)
      #println("tf:",typeof(tf))
end =#


# Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
# Create Hankel Matrix, perform svd to determine the highest order singular values
# Create Observibility and Control Matrices -> Create A, B, and C 
# Transform A,B,C matrices to final form
# Create D matrice from H -> ∞

end