using LiiBRA, UnitSystems, Test, LinearAlgebra

#Construct cell struct
Cell = Construct("LG_M50")
Cell.RA.Tlen = 128
Cell.RA.Fs = 1
Cell.RA.H1 = 0:Cell.RA.Tlen
Cell.RA.H2 = 0:Cell.RA.Tlen

function TestDRA(Cell)
    A = B = C = D = Time = x = tuple()
    for Temp in 5.0:10.0:45.0
        #Arrhenius
        Cell.Const.T = 273.15+Temp
        Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R

        #Set Cell Constants
        Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
        Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
        Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
        Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)

        #DRA
        A_DRA, B_DRA, C_DRA, D_DRA = DRA(Cell,Cell.RA.s,Cell.RA.f)

        #Flatten output into tuples
        A = flatten_(A,A_DRA)
        B = flatten_(B,B_DRA)
        C = flatten_(C,C_DRA)
        D = flatten_(D,D_DRA)
    end
return A, B, C, D
end

A,B,C,D = TestDRA(Cell)


@test size(A[1],1) == Cell.RA.M+1
@test typeof(A[1]) == Diagonal{Float64, Vector{Float64}}
@test typeof(B[1]) == Matrix{Float64}
@test typeof(C[1]) == Matrix{Float64}
@test typeof(D[1]) == Matrix{Float64}
#@test output eigs and compare to diagonal of A 
#@test sim_model output at multiple SOC?

