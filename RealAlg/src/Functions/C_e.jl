function C_e(CellData::Cell,s,z,M)
""" 
Electrolyte Concentration Transfer Function
# Add License
# Add Ins and Outs
    # Locations to be computed
    # Sampling Frequency
"""

T = CellData.Const.T      # Temperature
t_plus = CellData.Const.t_plus  # Transference Number
ζ = (1-t_plus)/F    #Simplifying Variable
Ds_Neg = CellData.Neg.Ds       # Solid diffusivity [m^2/s]
Ds_Pos = CellData.Pos.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Geo.CC_A   # Current-collector area [m^2]
κ_eff_Neg = CellData.Const.κ*ϵ1^CellData.Neg.κ_brug
κ_eff_Pos = CellData.Const.κ*ϵ3^CellData.Pos.κ_brug
σ_eff_Neg = CellData.Neg.σ*ϵ1^CellData.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = CellData.Pos.σ*ϵ3^CellData.Pos.σ_brug #Effective Conductivity Pos

#Defining SOC
θ_neg = CellData.Const.Init_SOC * (CellData.Neg.θ_max-CellData.Neg.θ_min) + CellData.Neg.θ_min 
θ_pos = CellData.Const.Init_SOC * (CellData.Pos.θ_max-CellData.Pos.θ_min) + CellData.Pos.θ_min 

#Beta's
βn = @. Rs_Neg*sqrt(s/Ds_Neg)
βn = permutedims(βn)
βp = @. Rs_Pos*sqrt(s/Ds_Pos)
βp = permutedims(βp)
println("β:",size(βn))

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max_neg = CellData.Neg.cs_max
cs0_neg = cs_max_neg * θ_neg
cs_max_pos = CellData.Pos.cs_max
cs0_pos = cs_max_pos * θ_pos

α_neg = CellData.Neg.α
α_pos = CellData.Pos.α

#Current Flux Density
j0_neg = CellData.Neg.k_norm*(ce0*(cs_max_neg-cs0_neg))^(1-α_neg)*cs0_neg^α_neg
j0_pos = CellData.Pos.k_norm*(ce0*cs_max_pos*cs0_pos)^(1-α_pos)*cs0_pos^α_pos

#Resistances
Rct_neg = R*T/(j0_neg*F)^2
Rtot_neg = Rct_neg + CellData.Neg.RFilm

Rct_pos = R*T/(j0_pos*F)^2
Rtot_pos = Rct_pos + CellData.Pos.RFilm

#∂Uocp_neg = UOCP(θ_neg)
∂Uocp_neg = (-20000*exp(-2000*θ_neg) - 3.96*exp(-3*θ_neg))
#∂Uocp_pos = UOCP(θ_pos)
∂Uocp_pos = (-32.4096*exp(-40*(-0.133875 + θ_pos)) - 0.0135664./((0.998432 - θ_pos).^1.49247)+ 0.0595559*exp(-0.04738*θ_pos.^8).*θ_pos.^7 - 0.823297*(sech(8.60942 - 14.5546*θ_pos)).^2)

s = permutedims(s)
#Condensing Variable
ν_n =  @. Lneg*sqrt((as_neg/σ_eff_Neg+as_neg/κ_eff_Neg)/(Rtot_neg+∂Uocp_neg*(Rs_Neg/(F*Ds_Neg))*(tanh(βn)/(tanh(βn)-βn))))
ν_p =  @. Lpos*sqrt((as_pos/σ_eff_Pos+as_pos/κ_eff_Pos)/(Rtot_pos+∂Uocp_pos*(Rs_Pos/(F*Ds_Pos))*(tanh(βp)/(tanh(βp)-βp))))

λ = roots(M+1)
println("λ:",size(λ))
λ = (λ[1:size(λ,1) .!= 1,: ]) #Delete first element relating to location zero

#Create all k's
in1 = sqrt.(λ.*ϵ1./D1)
in2 = sqrt.(λ*ϵ2./D2)
in3 = sqrt.(λ*ϵ3./D3)

Bound_Neg_1 = Lneg * in1
Bound_Sep_0 = Lneg * in2
Bound_Sep_1 = (Lneg+Lsep) * in2
Bound_Pos_0 = (Lneg+Lsep) * in3
Bound_Pos_1 = Ltot * in3
Bound_Pos_2 = Lpos * in3

#Scaled coefficients
k3_s = cos.(Bound_Neg_1).*cos.(Bound_Sep_0)+D1*in1.*sin.(Bound_Neg_1).*sin.(Bound_Sep_0)./(D2*in2)
k4_s = cos.(Bound_Neg_1).*sin.(Bound_Sep_0) - D1*in1.*cos.(Bound_Sep_0).*sin.(Bound_Neg_1)./(D2*in2);
k5_s = k3_s.*(cos.(Bound_Sep_1).*cos.(Bound_Pos_1)+D2*in2.*sin.(Bound_Sep_1).*sin.(Bound_Pos_1)./(D3*in3))+k4_s.*(sin.(Bound_Sep_1).*cos.(Bound_Pos_1)-D2*in2.*cos.(Bound_Sep_1).*sin.(Bound_Pos_1)./(D3*in3));
k6_s = k3_s.*(cos.(Bound_Sep_1).*sin.(Bound_Pos_1)-D2*in2.*sin.(Bound_Sep_1).*cos.(Bound_Pos_1)./(D3*in3))+k4_s.*(sin.(Bound_Sep_1).*sin.(Bound_Pos_1)+D2*in2.*cos.(Bound_Sep_1).*cos.(Bound_Pos_1)./(D3*in3));

 println("ν_p:",size(ν_p))
 println("λ:",size(λ))
 println("in1:",size(in1))
# println("in2:",in2)
 println("Bound_Neg_1:",size(Bound_Neg_1))
# println("k3_s:",k3_s)
# println("βn:",length(βn))

#Solving for k1:
Int_ψ1 = ϵ1*(2*Bound_Neg_1+sin.(2*Bound_Neg_1)./(4*in1))
Int_ψ2 = ϵ2./(4*in2) .* (2 .* (k3_s.^2 .+ k4_s.^2) .* Lsep .* in2 + 2*k3_s .* k4_s .* cos.(2 .* Bound_Sep_0) .- 2 .* k3_s .* k4_s .* cos.(2 .* Bound_Sep_1) .- (k3_s .- k4_s) .* (k3_s .+ k4_s) .* (sin.(2 .* Bound_Sep_0) .- sin.(2 .* Bound_Sep_1)))
Int_ψ3 = ϵ3./(4*in3) .* (2 .* (k5_s.^2+k6_s.^2) .* Lpos .* in3 + 2*k5_s .* k6_s .* cos.(2 .* Bound_Pos_0) .- 2 .* k5_s .* k6_s .* cos.(2 .* Bound_Pos_1) .- (k5_s .- k6_s) .* (k5_s .+ k6_s) .* (sin.(2 .* Bound_Pos_0) .- sin.(2 .* Bound_Pos_1)))

# println("Int_ψ1:",Int_ψ1)
# println("Int_ψ2:",Int_ψ2)
# println("Int_ψ3:",Int_ψ3)

k1 = @. 1/(sqrt(Int_ψ1+Int_ψ2+Int_ψ3))
k3 = @. k1*k3_s
k4 = @. k1*k4_s
k5 = @. k1*k5_s
k6 = @. k1*k6_s

 println("k1:",size(k1))
# println("ζ:", ζ)

j_Neg = @. ((κ_eff_Neg+σ_eff_Neg)*cosh(ν_n)*ν_n)*(k1*ζ*Bound_Neg_1*sin(Bound_Neg_1))/(ν_n^2*sinh(ν_n)+CC_A*(κ_eff_Neg+σ_eff_Neg)*((Bound_Neg_1^2)))

# println("j_Neg:", length(j_Neg))
println("Bound_Pos_2:", size(Bound_Pos_2))
Hlp1 = σ_eff_Pos+κ_eff_Pos
Hlp2 = @. (Hlp1*cosh(ν_p)*ν_p)
Hlp3 = @. (Bound_Pos_2^2)+ν_p^2*sinh(ν_p)

j_Pos1 = @. (k6*ζ*Bound_Pos_2*cos(Bound_Pos_1)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos2 = @. (k5*ζ*Bound_Pos_2*sin(Bound_Pos_2)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos3 = @. (k6*ζ*Bound_Pos_2*cos(Bound_Pos_0)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos4 = @. (k5*ζ*Bound_Pos_2*sin(Bound_Pos_1)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos5 = @. (k5*ζ*σ_eff_Pos*cos(Bound_Pos_0)*κ_eff_Pos*cos(Bound_Pos_1)*ν_p^2)/(CC_A*Hlp1*(Bound_Pos_2^2 + ν_p^2))
j_Pos6 = @. (k6*ζ*σ_eff_Pos*sin(Bound_Pos_0)*κ_eff_Pos*sin(Bound_Pos_1)*ν_p^2)/(CC_A*Hlp1*(Bound_Pos_2^2 + ν_p^2))

j_Pos = j_Pos1 - j_Pos2 + j_Pos3 - j_Pos4 - j_Pos5 - j_Pos6
println("j_Pos:",size(j_Pos))

C_e =  @. ((j_Neg + j_Pos)/(s+λ))

i=1
ψ = fill(0.0,1,length(z))
println("ψ:",size(ψ))
for x in z
    if x < Lneg
       ψ[i] = k1[i]*cos(in1[i]*x) #negative electrode
    elseif x < Lnegsep
        ψ[i] = @. k3[i]*cos(in2[i]*x)+k4[i]*sin(in2[i]*x) # separator
    else
        ψ[i] = @. k5[i]*cos(in3[i]*x)+k6[i]*sin(in3[i]*x) # postive electrode
    end
i = i+1
end
ψ_dims = ψ.*diagm(0=>fill(1., size(ψ,2)))
println("ψ_dims:",size(ψ_dims))
C_e = ψ_dims*C_e
println("Ce:",size(C_e))
return C_e
end

function roots(roots_n)
    root = Any[0]
    i = 0.00001
    ∇ = 0.00001
    if(roots_n > 1)
        while length(root) <= roots_n-1
            if flambda(i-∇)*flambda(i)<0
            push!(root, find_zero(flambda,(i-∇,i)))
            end
        i = i+∇
        end
    end
    return root
end

function flambda(λ)
    k1 = 1
    sle1 = sqrt(λ*ϵ1/D1)
    sle2 = sqrt(λ*ϵ2/D2)
    sle3 = sqrt(λ*ϵ3/D3)
    k3 = k1*(cos(sle1*Lneg).*cos(sle2*Lneg) + D1*sle1.*sin(sle1*Lneg).*sin(sle2*Lneg)./(D2*sle2))
    k4 = k1*(cos(sle1*Lneg).*sin(sle2*Lneg) - D1*sle1.*cos(sle2*Lneg).*sin(sle1*Lneg)./(D2*sle2))
    k5 = k3*(cos(sle2*Lnegsep).*cos(sle3*Lnegsep) + D2*sle2.*sin(sle2*Lnegsep).*sin(sle3*Lnegsep)./(D3*sle3))+k4*(sin(sle2*Lnegsep).*cos(sle3*Lnegsep) - D2*sle2.*cos(sle2*Lnegsep).*sin(sle3*Lnegsep)./(D3*sle3))
    k6 = k3*(cos(sle2*Lnegsep).*sin(sle3*Lnegsep) - D2*sle2.*sin(sle2*Lnegsep).*cos(sle3*Lnegsep)./(D3*sle3))+k4*(sin(sle2*Lnegsep).*sin(sle3*Lnegsep) + D2*sle2.*cos(sle2*Lnegsep).*cos(sle3*Lnegsep)./(D3*sle3))
    Psiprime = -k5.*sle3.*sin(sle3*Ltot) + k6.*sle3.*cos(sle3*Ltot)
  end