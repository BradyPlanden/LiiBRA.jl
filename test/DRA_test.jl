using LiiBRA, UnitSystems, Test, LinearAlgebra, JLD2

#Construct cell struct
Cell = Construct("LG M50")
Cell.RA.Tlen = 128
Cell.RA.Fs = 1
Cell.RA.H1 = 1:Cell.RA.Tlen
Cell.RA.H2 = 1:Cell.RA.Tlen

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
        Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
        Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)

        #DRA
        A_DRA, B_DRA, C_DRA, D_DRA = CIDRA(Cell)

        #Flatten output into tuples
        A = flatten_(A,A_DRA)
        B = flatten_(B,B_DRA)
        C = flatten_(C,C_DRA)
        D = flatten_(D,D_DRA)
    end
return A, B, C, D
end

#Load Test Data
Data = load("CIDRA_Data.jld2")
A,B,C,D = TestDRA(Cell)


@test size(A[1],1) == Cell.RA.M+1
@test typeof(A[1]) == Matrix{Float64}
@test typeof(B[1]) == Vector{Float64}
@test typeof(C[1]) == Matrix{Float64}
@test typeof(D[1]) == Vector{Float64}

@test Data["A0"] ≈ A[1]
@test Data["B0"] ≈ B[1]
@test Data["C0"] ≈ C[1]
@test Data["D0"] ≈ D[1]

#ToDo
#@test output eigs and compare to diagonal of A 
#@test sim_model output at multiple SOC?

