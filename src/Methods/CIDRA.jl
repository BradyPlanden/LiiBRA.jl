function CIDRA(Cell)
    """ 
    Function for Computationally Informed Discrete Realisation Algorithm.

    CIDRA(Cell)
    """

    # Method
    if Cell.Const.Debug == true
        if Cell.RA.Fs == (1 / Cell.RA.SamplingT)
            println("CIDRA Utilised")
        else
            println("DRA Utilised")
        end
    end

    # Additional Pulse Setup
    tfft = (1 / Cell.RA.Fs) * Cell.RA.f
    OrgT = Cell.RA.SamplingT * (0:floor(tfft[end] / Cell.RA.SamplingT))

    # Initialise Loop Variables
    A = Matrix{Float64}(I, Cell.RA.M + 1, Cell.RA.M + 1)
    B = Vector{Float64}(undef, Cell.RA.M + 1)
    C = Matrix{Float64}(undef, Cell.RA.Outs, Cell.RA.M + 1)
    D = Vector{Float64}(undef, Cell.RA.Outs)
    puls = Array{Float64}(undef, Cell.RA.Outs, size(OrgT, 1) - 1)
    Cₐ = Vector{Float64}(undef, Cell.RA.Outs)
    Sₑ = Cell.Transfer.Sₑ
    Sₛ = Cell.Transfer.Sₛ
    i = Int(1)
    l = Int(1)
    u = Int(0)

    for run in Cell.Transfer.tfs
        tf = Array{ComplexF64}(undef, size(Cell.Transfer.Locs(Sₑ, Sₛ)[i], 1),
                               size(Cell.RA.s, 2))
        Dᵢ = Vector{Float64}(undef, size(Cell.Transfer.Locs(Sₑ, Sₛ)[i], 1))::Vector{Float64}
        res₀ = Vector{Float64}(undef,
                               size(Cell.Transfer.Locs(Sₑ, Sₛ)[i], 1))::Vector{Float64}
        smptf = Array{Float64}(undef, size(Cell.Transfer.Locs(Sₑ, Sₛ)[i], 1),
                               length(OrgT))
        u += Int(size(Cell.Transfer.Locs(Sₑ, Sₛ)[i], 1))

        if Cell.Transfer.Elec[i] == "Pos"
            run(Cell, Cell.RA.s, Cell.Transfer.Locs(Sₑ, Sₛ)[i], "Pos", tf, Dᵢ, res₀)
        elseif Cell.Transfer.Elec[i] == "Neg"
            run(Cell, Cell.RA.s, Cell.Transfer.Locs(Sₑ, Sₛ)[i], "Neg", tf, Dᵢ, res₀)
        else
            run(Cell, Cell.RA.s, Cell.Transfer.Locs(Sₑ, Sₛ)[i], tf, Dᵢ, res₀)
        end

        if Cell.RA.Fs == (1 / Cell.RA.SamplingT)
            puls[l:u, :] .= real(ifft(tf, 2))[:, 2:end]
        else
            #Interpolate H(s) to obtain h_s(s) to obtain discrete-time impulse response
            stpsum = cumsum(real(ifft(tf, 2)), dims = 2)
            for k in 1:size(stpsum, 1)
                spl₁ = CubicSplineInterpolation(tfft, stpsum[k, :]; bc = Line(OnGrid()),
                                                extrapolation_bc = Throw())
                smptf[k, :] = spl₁(OrgT)
            end
            puls[l:u, :] .= diff(smptf, dims = 2)
        end

        D[l:u, :] .= Dᵢ
        Cₐ[l:u, :] .= res₀
        l = u + one(u)
        i += one(i)
    end

    # Scale Transfer Functions in Pulse Response
    S̃ = sqrt.(sum(puls .^ 2, dims = 2))
    puls .= puls ./ S̃

    # Pre-Allocation for Hankel & SVD
    PulsL = size(puls, 1)
    𝐇 = Array{Float64}(undef, length(Cell.RA.H1) * PulsL, length(Cell.RA.H2))
    U, S, V = fh!(𝐇, Cell.RA.H1, Cell.RA.H2, puls, Cell.RA.M, PulsL)

    # Create Observibility and Control Matrices -> Create A, B, and C 
    ʃ = sqrt(Diagonal(S))
    U .= U * ʃ
    V .= ʃ * V

    A[2:end, 2:end] .= (U \ 𝐇) / V
    B .= [Cell.RA.SamplingT; V[:, 1:(Cell.RA.N)]]
    C .= [Cₐ S̃ .* U[1:PulsL, :]]

    # System checks
    if any(i -> i > 1.0, real(eigvals(A)))
        println("Oscilating System: A has indices of values > 1")
    end

    if any(i -> i < 0.0, real(eigvals(A)))
        println("Unstable System: A has indices of negative values")
    end

    if conj(eigvals(A)) != eigvals(A)
        println("System has non-real eigenvalues at SOC:$(Cell.Const.SOC)")
    end

    # Transform the SS system for interpolation
    d, Sᵘ = eigen(A, sortby = nothing)
    A = inv(Sᵘ) * A * Sᵘ
    B = inv(Sᵘ) * B
    C = C * Sᵘ * Diagonal(B)
    B = ones(length(B))

    return mag(A), mag(B), mag(C), D
end
