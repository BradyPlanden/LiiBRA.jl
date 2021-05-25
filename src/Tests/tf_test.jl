using RealAlg, CSV, DataFrames

# Cell_Df = CSV.read("/Users/bradyplanden/Documents/Git/RealisationAlgorithms/RealAlg/Data/Sample-Cell.csv", DataFrame)
# Constants_Df = tuple.(Tuple.(eachrow(dropmissing(Cell_Df[:,1:2])))...)
# Geometry_Df = tuple.(Tuple.(eachrow(dropmissing(Cell_Df[:,3:4])))...)
# Negative_Df = tuple.(Tuple.(eachrow(dropmissing(Cell_Df[:,5:6])))...)
# Positive_Df = tuple.(Tuple.(eachrow(dropmissing(Cell_Df[:,7:8])))...)
# Seperator_Df = tuple.(Tuple.(eachrow(dropmissing(Cell_Df[:,9:10])))...)

# CellData = Cell(Constants(1,1.2,1.3,1.4,1.5,1.6),Geometry(1,1.2,1.3,1.4,1.5),Negative(1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,11.11,1.12,1.13,1.14),Positive(1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,11.11,1.12,1.13,1.14),Seperator(1,1.2))
CellData = Cell(Constants(),Geometry(),Negative(),Positive(),Seperator())

M=4
Fs = 2
Tlen = 32768
Nfft = 2^(ceil(log2(Fs*Tlen)))
f = 0:Nfft-1
s = (2im.*Fs)*tan.(pi.*f./Nfft)

# loc = Any[0,1]
loc = Any[0, 128e-6, 204e-6, 394e-6]
Def = "Neg"

test_ce = C_e(CellData,s,M)
println("Ce:",length(test_ce))
test_j = j(CellData,s,loc,Def)
println("j:",typeof(test_j))
test_cse = C_se(CellData,s,loc,Def)
println("Cse:",typeof(test_cse))
test_ϕ = Phi_s(CellData,s,loc,Def)
println("ϕ_s:",typeof(test_ϕ))
test_ϕ_e = Phi_e(CellData,s,loc,Def)
println("ϕ_e:", length(test_ϕ_e))

# c = Cell_Df[:,2]
# c = dropmissing(Cell_Df)
# x = copy(Cell_Df)
# y = Tuple(x)


# c = dropmissing(Cell_Df[:,1:2])
# I1 = [x[1] for x in j]

# q1 = Tuple.(eachrow(dropmissing(Cell_Df[:,1:2])))
# q1 = tuple.(q1...)

