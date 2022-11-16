using BenchmarkTools, LiiBRA
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
Cell = Construct("LG M50")
Ŝ = SOC = 0.8

function DRA_Loop(Cell, Ŝ)

    # Arrhenius
    Cell.Const.T = 298.15
    Arr_Factor = (1/Cell.Const.T_ref-1/Cell.Const.T)/R
    
    # Realisation Variables
    Cell.RA.H1 = [1:2000; 3000:3500]
    Cell.RA.H2 = [1:2000; 3000:3500]
    Cell.RA.Fs = 4
    Cell.RA.Tlen = 16200
    Cell.RA.M = 4
    Cell.RA.SamplingT = 1/4
    
    # Set Cell Constants
    Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0)*exp(Cell.Const.Ea_κ*Arr_Factor)
    Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
    Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
    Cell.RA.s = Cell.RA.s!(Cell.RA.Fs,Cell.RA.Nfft,Cell.RA.f)
    Cell.Neg.β = Cell.Neg.β!(Cell.RA.s)
    Cell.Pos.β = Cell.Pos.β!(Cell.RA.s)
    
    return @benchmark CIDRA(Cell)
end
   
Time = DRA_Loop(Cell, Ŝ)
