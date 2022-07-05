using BenchmarkTools, JLD2, Arpack, TSVD # run PROPACK, tsvd (remove k=) then change to svds and do ARPACK
import PROPACK.tsvd as tsvdpro

H1 = 0:4500
H2 = 0:4500


BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
puls = load("puls_57600.jld2", "puls")
Puls_L = size(puls,1)
Hank = Array{Float64}(undef,length(H1)*Puls_L,length(H2))

#---------- Cell Definition -----------------#
# Cell = Construct("LG M50") #Alternative "Doyle_94"
# Init_SOC = 0.75


# #Arrhenius
# Cell.Const.T = 298.15
# Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R

# #Set Cell Constants
# Cell.Const.SOC = Init_SOC
# Cell.RA.Tlen = 57600
# Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
# Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
# Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
# Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)
# Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
# Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)


function SVD!(Hank,puls,Hlen1,Hlen2)
    A = T = P = tuple()
    Puls_L = size(puls,1)
    for lp1 in 1:length(Hlen2), lp2 in 1:length(Hlen1)
        Hank[Puls_L*(lp2-1)+1:Puls_L*lp2,lp1] .= @view puls[:,Hlen2[lp1]+Hlen1[lp2]+1]
    end

    A = @benchmark svds(Hank; nsv=6)[1]
    T = @benchmark tsvd(Hank, 6)
    P = @benchmark tsvdpro(Hank, k=6)
    return A,T,P
end


# @inline function DRA_loop(Cell)
#     x = Time = tuple()
#     for i in 1500:500:5000
#         Cell.RA.H1 = 0:i
#         Cell.RA.H2 = 0:i
#         x = @benchmark DRA(Cell,Cell.RA.s,Cell.RA.f)
#         Time = flatten_(Time,x)
#     end
#     return Time
# end

A,T,P = SVD!(Hank,puls,H1,H2)
