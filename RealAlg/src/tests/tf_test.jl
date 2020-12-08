using RealAlg

Cell_Df = CSV.read("/Users/bradyplanden/Documents/Git/RealisationAlgorithms/RealAlg/Data/Sample-Cell.csv", DataFrame)
CellData = Cell(Constants(1,1.2,1.3,1.4,1.5,1.6),Geometry(1,1.2,1.3,1.4,1.5),Negative(1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,11.11,1.12,1.13,1.14),Positive(1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,11.11,1.12,1.13,1.14),Seperator(1,1.2))

Fs = 2
Tlen = 32768
Nfft = 2^(ceil(log2(Fs*Tlen)))
const f = 0:Nfft-1
const s = (2im.*Fs)*tan.(pi.*f./Nfft)

test_ce = C_e(CellData,s)