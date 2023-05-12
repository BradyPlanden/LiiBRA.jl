using JLD2, Dierckx, FFTW, StatsBase, Infiltrator

#Construct Cell
Cell = Construct("LG M50")

#Arrhenius
Cell.Const.T = 298.15
Arr_Factor = (1 / Cell.Const.T_ref - 1 / Cell.Const.T) / R

#Set Cell Constants
Cell.Const.SOC = 0.8
Cell.Const.κ = Cell.Const.κf(Cell.Const.ce0) * exp(Cell.Const.Ea_κ * Arr_Factor)
Cell.RA.Nfft = Cell.RA.Nfft!(Cell.RA.Fs, Cell.RA.Tlen)
Cell.RA.f = Cell.RA.f!(Cell.RA.Nfft)
Cell.RA.s = Cell.RA.s!(Cell.RA.Fs, Cell.RA.Nfft, Cell.RA.f)

tfft = @. (1 / Cell.RA.Fs) * Cell.RA.f
OrgT = @. Cell.RA.SamplingT * (0:floor(tfft[end] / Cell.RA.SamplingT))
puls = Array{Float64}(undef, Cell.RA.Outs, size(OrgT, 1) - 1)
D = Vector{Float64}(undef, Cell.RA.Outs)
C_Aug = Vector{Float64}(undef, Cell.RA.Outs)

function tf_new!(puls, D, C_Aug, s, f)
    #Initialise Loop Variables
    i = Int64(1)
    l = Int64(1)
    u = Int64(0)

    for run in Cell.Transfer.tfs[:, 1]
        tf = Array{ComplexF64}(undef, size(Cell.Transfer.tfs[i, 3], 1), size(s, 2))
        Di = Vector{Float64}(undef, size(Cell.Transfer.tfs[i, 3], 1))
        res0 = Vector{Float64}(undef, size(Cell.Transfer.tfs[i, 3], 1))

        if Cell.Transfer.tfs[i, 2] == "Pos"
            run(Cell, s, Cell.Transfer.tfs[i, 3], "Pos", tf, Di, res0)
        elseif Cell.Transfer.tfs[i, 2] == "Neg"
            run(Cell, s, Cell.Transfer.tfs[i, 3], "Neg", tf, Di, res0)
        else
            run(Cell, Cell.RA.s, Cell.Transfer.tfs[i, 3], tf, Di, res0)
        end
        jk = Cell.RA.Fs * real(ifft(tf, 2)) # inverse fourier transform tranfser function response (Large Compute)
        stpsum = (cumsum(jk, dims = 2) .* (1 / Cell.RA.Fs)) #cumulative sum of tf response * sample time
        samplingtf = Array{Float64}(undef, size(stpsum, 1), length(OrgT))

        #Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
        for Output in 1:size(stpsum, 1)
            spl1 = Spline1D(tfft, stpsum[Output, :]; k = 3)
            samplingtf[Output, :] .= evaluate(spl1, OrgT)
        end
        u += size(Cell.Transfer.tfs[i, 3], 1)
        #puls = [puls; diff(samplingtf, dims=2)]
        puls[l:u, :] .= diff(samplingtf, dims = 2)
        D[l:u, :] .= Di
        C_Aug[l:u, :] .= res0
        l = u + 1
        i += 1
    end
end

function tf(Cell, s, f)
    tfft = @. (1 / Cell.RA.Fs) * Cell.RA.f
    OrgT = @. Cell.RA.SamplingT * (0:floor(tfft[end] / Cell.RA.SamplingT))
    i = Int64(1)
    puls = Array{Float64}(undef, 0, size(OrgT, 1) - 1)
    D = Array{Float64}(undef, 0, 1)
    Dtt = Array{String}(undef, 0, 1)
    C_Aug = Array{Float64}(undef, 0, 1)
    for run in Cell.Transfer.tfs[:, 1]
        if Cell.Transfer.tfs[i, 2] == "Pos"
            tf, Di, res0, Dti = run(Cell, s, Cell.Transfer.tfs[i, 3], "Pos")
        elseif Cell.Transfer.tfs[i, 2] == "Neg"
            tf, Di, res0, Dti = run(Cell, s, Cell.Transfer.tfs[i, 3], "Neg")
        else
            tf, Di, res0, Dti = run(Cell, s, Cell.Transfer.tfs[i, 3])
        end

        jk = Cell.RA.Fs * real(ifft(tf, 2)') # inverse fourier transform tranfser function response (Large Compute)
        stpsum = (cumsum(jk, dims = 1) .* (1 / Cell.RA.Fs))' # cumulative sum of tf response * sample time
        samplingtf = Array{Float64}(undef, size(stpsum, 1), length(OrgT))

        # Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
        for Output in 1:size(stpsum, 1)
            spl1 = Spline1D(tfft, stpsum[Output, :]; k = 3)
            samplingtf[Output, :] = evaluate(spl1, OrgT)
        end

        puls = [puls; diff(samplingtf, dims = 2)]
        D = [D; Di]
        Dtt = [Dtt; Dti]
        C_Aug = [C_Aug; res0]
        i += 1
    end
    return puls
end

#Main Branch
#puls_ref = tf(Cell,Cell.RA.s,Cell.RA.f)

#Dev Branch
puls_ref = load("puls.jld2", "puls_ref")
tf_new!(puls, D, C_Aug, Cell.RA.s, Cell.RA.f)

maxdiff = Vector{Float64}(undef, size(puls_ref, 1))
rmsediff = Vector{Float64}(undef, size(puls_ref, 1))
for i in 1:size(puls_ref, 1)
    maxdiff[i] = maxad(puls[i, :], puls_ref[i, :])
    rmsediff[i] = rmsd(puls[i, :], puls_ref[i, :])
end
#jldsave("puls.jld2"; puls_ref)
