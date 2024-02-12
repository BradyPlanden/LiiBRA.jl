using BenchmarkTools, JLD2, Arpack, TSVD # run PROPACK, tsvd (remove k=) then change to svds and do ARPACK
import PROPACK.tsvd as tsvdpro

H1 = 0:4500
H2 = 0:4500
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
puls = load("puls_57600.jld2", "puls")
Puls_L = size(puls, 1)
Hank = Array{Float64}(undef, length(H1) * Puls_L, length(H2))

function SVD!(Hank, puls, Hlen1, Hlen2)
    A = T = P = tuple()
    Puls_L = size(puls, 1)
    for lp1 in 1:length(Hlen2), lp2 in 1:length(Hlen1)
        Hank[(Puls_L * (lp2 - 1) + 1):(Puls_L * lp2), lp1] .= @view puls[:,
            Hlen2[lp1] + Hlen1[lp2] + 1]
    end

    A = @benchmark svds(Hank; nsv = 6)[1]
    T = @benchmark tsvd(Hank, 6)
    P = @benchmark tsvdpro(Hank, k = 6)
    return A, T, P
end

A, T, P = SVD!(Hank, puls, H1, H2)
