@inline function C_e(CellData::Cell,s,z)
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
σ_eff_Neg = CellData.Neg.σ*CellData.Neg.ϵ_s^CellData.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = CellData.Pos.σ*CellData.Pos.ϵ_s^CellData.Pos.σ_brug #Effective Conductivity Pos

#Defining SOC
θ_neg = CellData.Const.SOC * (CellData.Neg.θ_100-CellData.Neg.θ_0) + CellData.Neg.θ_0 
θ_pos = CellData.Const.SOC * (CellData.Pos.θ_100-CellData.Pos.θ_0) + CellData.Pos.θ_0 

#Beta's
βn = @. CellData.Neg.Rs*sqrt(s/Ds_Neg)
βp = @. CellData.Pos.Rs*sqrt(s/Ds_Pos)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max_neg = CellData.Neg.cs_max
cs0_neg = cs_max_neg * θ_neg
cs_max_pos = CellData.Pos.cs_max
cs0_pos = cs_max_pos * θ_pos

#Current Flux Density
κ_pos = CellData.Pos.k_norm/CellData.Pos.cs_max/ce0^(1-CellData.Pos.α)
κ_neg = CellData.Neg.k_norm/CellData.Neg.cs_max/ce0^(1-CellData.Neg.α)
j0_neg = κ_neg*(ce0*(cs_max_neg-cs0_neg))^(1-CellData.Neg.α)*cs0_neg^CellData.Neg.α
j0_pos = κ_pos*(ce0*(cs_max_pos-cs0_pos))^(1-CellData.Pos.α)*cs0_pos^CellData.Pos.α

#Resistances
Rtot_neg = R*CellData.Const.T/(j0_neg*F^2) + CellData.Neg.RFilm
Rtot_pos = R*CellData.Const.T/(j0_pos*F^2) + CellData.Pos.RFilm

#OCP derivative
∂Uocp_pos = ∂Uocp("Pos",θ_pos)/cs_max_pos
∂Uocp_neg = ∂Uocp("Neg",θ_neg)/cs_max_neg

#Condensing Variable
ν_n =  @. Lneg*sqrt((as_neg/σ_eff_Neg+as_neg/κ_eff_Neg)/(Rtot_neg+∂Uocp_neg*(CellData.Neg.Rs/(F*Ds_Neg))*(tanh(βn)/(tanh(βn)-βn))))
ν_p =  @. Lpos*sqrt((as_pos/σ_eff_Pos+as_pos/κ_eff_Pos)/(Rtot_pos+∂Uocp_pos*(CellData.Pos.Rs/(F*Ds_Pos))*(tanh(βp)/(tanh(βp)-βp))))
# println("ν_n:")
# display("text/plain", ν_n[:,1:10])

# println("ζ:",ζ)
# println("κ_eff_pos:",κ_eff_Pos)
# println("κ_eff_neg:",κ_eff_Neg)
# println("σ_eff_Pos:",σ_eff_Pos)
# println("σ_eff_Neg:",σ_eff_Neg)
# println("∂Uocp_pos:",∂Uocp_pos)
# println("Rtot_pos:",Rtot_pos)
# println("Ds_pos:",Ds_Pos)
# println("j0_pos:",j0_pos)

λ = roots(CellData.RA.M+1)
λ = (λ[1:size(λ,1) .!= 1,: ]) #Delete first element relating to location zero
# println("λ:",size(λ))
# println("ν_n:",size(ν_n))
# println("∂Uocp:",∂Uocp_neg)

#Create all k's
in1 = @. sqrt(λ*ϵ1/D1)
in2 = @. sqrt(λ*ϵ2/D2)
in3 = @. sqrt(λ*ϵ3/D3)
# println("in1",in1)
# println("in2",in2)
# println("in3",in3)

Bound_Neg_1 = Lneg * in1
Bound_Sep_0 = Lneg * in2
Bound_Sep_1 = (Lneg+Lsep) * in2
Bound_Pos_0 = (Lneg+Lsep) * in3
Bound_Pos_1 = Ltot * in3
Bound_Pos_2 = Lpos * in3

# println("Bound_Neg_1:",Bound_Neg_1)
# println("Bound_Sep_0:",Bound_Sep_0)
# println("Bound_Sep_1:",Bound_Sep_1)
# println("Bound_Pos_0:",Bound_Pos_0)
# println("Bound_Pos_1:",Bound_Pos_1)
# println("Bound_Pos_2:",Bound_Pos_2)

#Scaled coefficients
k3_s = cos.(Bound_Neg_1).*cos.(Bound_Sep_0)+D1*in1.*sin.(Bound_Neg_1).*sin.(Bound_Sep_0)./(D2*in2)
k4_s = cos.(Bound_Neg_1).*sin.(Bound_Sep_0) - D1*in1.*cos.(Bound_Sep_0).*sin.(Bound_Neg_1)./(D2*in2);
k5_s = k3_s.*(cos.(Bound_Sep_1).*cos.(Bound_Pos_0)+D2*in2.*sin.(Bound_Sep_1).*sin.(Bound_Pos_0)./(D3*in3))+k4_s.*(sin.(Bound_Sep_1).*cos.(Bound_Pos_0)-D2*in2.*cos.(Bound_Sep_1).*sin.(Bound_Pos_0)./(D3*in3));
k6_s = k3_s.*(cos.(Bound_Sep_1).*sin.(Bound_Pos_0)-D2*in2.*sin.(Bound_Sep_1).*cos.(Bound_Pos_0)./(D3*in3))+k4_s.*(sin.(Bound_Sep_1).*sin.(Bound_Pos_0)+D2*in2.*cos.(Bound_Sep_1).*cos.(Bound_Pos_0)./(D3*in3));
# println("k3_s:",k3_s)
# println("k4_s:",k4_s)
# println("k5_s:",k5_s)
# println("k6_s:",k6_s)

#Solving for k1:
Int_ψ1 = @. ϵ1*(2*Bound_Neg_1+sin(2*Bound_Neg_1))/(4*in1)
Int_ψ2 = @. ϵ2/(4*in2)*(2*(k3_s^2+k4_s^2)*Lsep*in2+2*k3_s*k4_s*cos(2*Bound_Sep_0)-2*k3_s*k4_s*cos(2*Bound_Sep_1)-(k3_s-k4_s)*(k3_s+k4_s)*(sin(2*Bound_Sep_0)-sin(2*Bound_Sep_1)))
Int_ψ3 = ϵ3./(4*in3) .* (2 .* (k5_s.^2+k6_s.^2) .* Lpos .* in3 + 2*k5_s .* k6_s .* cos.(2 .* Bound_Pos_0) .- 2 .* k5_s .* k6_s .* cos.(2 .* Bound_Pos_1) .- (k5_s .- k6_s) .* (k5_s .+ k6_s) .* (sin.(2 .* Bound_Pos_0) .- sin.(2 .* Bound_Pos_1)))


k1 = @. 1/(sqrt(Int_ψ1+Int_ψ2+Int_ψ3))
k3 = @. k1*k3_s
k4 = @. k1*k4_s
k5 = @. k1*k5_s
k6 = @. k1*k6_s

# println("k5:",k5)
# println("k6:",k6)
# println("k4:",k4)
# println("k3:",k3)
# println("k1:",k1)
# println("Int_ψ1:",Int_ψ1)
# println("Int_ψ2:",Int_ψ2)
# println("Int_ψ3:",Int_ψ3)

j_Neg = @. k1*ζ*ν_n.*(Bound_Neg_1*(κ_eff_Neg+σ_eff_Neg*cosh.(ν_n))*sin.(Bound_Neg_1)+(κ_eff_Neg+σ_eff_Neg*cos.(Bound_Neg_1)).*sinh.(ν_n).*ν_n) ./(CC_A*(κ_eff_Neg+σ_eff_Neg)*(Bound_Neg_1^2 +ν_n.^2).*sinh.(ν_n));
#j_Neg = @. (((κ_eff_Neg+σ_eff_Neg*cosh(ν_n))*ν_n)*(k1*ζ*Bound_Neg_1*sin(Bound_Neg_1)))/(sinh(ν_n)*CC_A*(κ_eff_Neg+σ_eff_Neg)*((Bound_Neg_1^2+ν_n^2))) + (((κ_eff_Neg+σ_eff_Neg*cosh(Bound_Neg_1))*ν_n^2)*(k1*ζ))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*((Bound_Neg_1^2+ν_n^2)))
zero_tf_neg = @. k1*ζ*sin(Bound_Neg_1)/(CC_A*Bound_Neg_1)
j_Neg[:,findall(s.==0)] .= zero_tf_neg[:,findall(s.==0)]
# println("size:jn:",size(j_Neg))

#= Hlp1 = @. ((κ_eff_Pos+σ_eff_Pos*cosh(ν_p))*ν_p)
Hlp2 = @. ((σ_eff_Pos+κ_eff_Pos*cosh(ν_p))*ν_p)
Hlp3 = @. (Bound_Pos_2^2+ν_p^2)*sinh(ν_p)

j_Pos1 = @. (k6*ζ*Bound_Pos_2*cos(Bound_Pos_1)*Hlp2)/(CC_A*(σ_eff_Pos+κ_eff_Pos)*Hlp3)
j_Pos2 = @. (k5*ζ*Bound_Pos_2*sin(Bound_Pos_2)*Hlp1)/(CC_A*(σ_eff_Pos+κ_eff_Pos)*Hlp3)
j_Pos3 = @. (k6*ζ*Bound_Pos_2*cos(Bound_Pos_0)*Hlp1)/(CC_A*(σ_eff_Pos+κ_eff_Pos)*Hlp3)
j_Pos4 = @. (k5*ζ*Bound_Pos_2*sin(Bound_Pos_1)*Hlp2)/(CC_A*(σ_eff_Pos+κ_eff_Pos)*Hlp3)
j_Pos5 = @. (k5*ζ*(σ_eff_Pos*cos(Bound_Pos_0)+κ_eff_Pos*cos(Bound_Pos_1))*ν_p^2)/(CC_A*(σ_eff_Pos+κ_eff_Pos)*(Bound_Pos_2^2+ν_p^2))
j_Pos6 = @. (k6*ζ*(σ_eff_Pos*sin(Bound_Pos_0)+κ_eff_Pos*sin(Bound_Pos_1))*ν_p^2)/(CC_A*(σ_eff_Pos+κ_eff_Pos)*(Bound_Pos_2^2+ν_p^2))

j_Pos = j_Pos1 - j_Pos2 + j_Pos3 - j_Pos4 - j_Pos5 - j_Pos6 =#

j_Pos = @. -ζ*ν_p/(CC_A*(κ_eff_Pos+σ_eff_Pos)*(Bound_Pos_2^2+ν_p^2)*sinh(ν_p))*(-k6*Bound_Pos_2*cos(Bound_Pos_1)*(σ_eff_Pos+κ_eff_Pos*cosh(ν_p))+Bound_Pos_2*(κ_eff_Pos+σ_eff_Pos*cosh(ν_p))*(k6*cos(Bound_Pos_0)-k5*sin(Bound_Pos_0))+k5*Bound_Pos_2*(σ_eff_Pos+ κ_eff_Pos*cosh(ν_p))*sin(Bound_Pos_1)+sinh(ν_p)*(k5*σ_eff_Pos*cos(Bound_Pos_0)+k5*κ_eff_Pos*cos(Bound_Pos_1)+k6*σ_eff_Pos*sin(Bound_Pos_0)+k6*κ_eff_Pos*sin(Bound_Pos_1))*ν_p)
zero_tf = @. -ζ*(k6*(cos(Bound_Pos_0)-cos(Bound_Pos_1))+k5*(sin(Bound_Pos_1)-sin(Bound_Pos_0)))/(CC_A*Bound_Pos_2)
j_Pos[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]

C_e =  ((j_Neg + j_Pos)./(s.+λ))

# println("C_e:")
# display("text/plain", C_e[:,1:10])

# println("j_Pos:")
# display("text/plain", j_Pos[:,1:10])

# println("j_Neg:")
# display("text/plain", j_Neg[:,1:10])
# println("zero_tf_neg:",zero_tf_neg)
# println("zero_tf:",zero_tf)
# println("λ:",length(λ))
i=1
ψ = fill(0.0,length(z),length(λ))
for loop in 1:length(λ)
    i=1
    for x in z #Eigen Weighting
        if x < Lneg+eps()
            ψ[i,loop] = k1[loop]*cos(in1[loop]*x) #negative electrode
        elseif x > Lnegsep-eps()
            ψ[i,loop] = k5[loop]*cos(in3[loop]*x)+k6[loop]*sin(in3[loop]*x) # postive electrode
        else
            ψ[i,loop] = k3[loop]*cos(in2[loop]*x)+k4[loop]*sin(in2[loop]*x) # separator
        end
    
    # println("x:",x)
    # println("lneg:",Lneg+eps())
    # println("lnegsep:",Lnegsep-eps())
    i = i+1
    end
end
# println("ψ:")
# display("text/plain", ψ)
Ce_tf = zeros(ComplexF64,CellData.RA.M,length(s))
#= for loop in 1:length(λ)
    for broad in 1:length(z)
    tmp = ψ.*C_e[broad,:]
    Ce_tf = tmp'+Ce_tf
    end
end =#

#ψ_dims = ψ.*diagm(0=>fill(1., size(ψ,2)))
Ce_tf = ψ*C_e
D_term = zeros(length(z))
res0 = zeros(length(z))

# print("z:",z, "\n")
# println("Ce_tf:",Ce_tf[:,1])

# println("βn:",size(βn)) #good
# println("j0_neg:",size(j0_neg)) #good
# println("θ_neg:",θ_neg) #good
# println("Rtot_neg:",Rtot_neg) #good
#  #good
# println("s:",size(s))
# println("j_Neg:",size(j_Neg), "\n") #no 
# println("j_Pos:",j_Pos[4,1:10], "\n") #no
# println("ν_n:",size(ν_n)) #no

return Ce_tf, D_term, res0
end

function roots(roots_n)
    root = Float64[0]
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
    k1 = 1.0
    sle1 = sqrt(λ*ϵ1/D1)
    sle2 = sqrt(λ*ϵ2/D2)
    sle3 = sqrt(λ*ϵ3/D3)
    k3 = k1*(cos(sle1*Lneg).*cos(sle2*Lneg) + D1*sle1.*sin(sle1*Lneg).*sin(sle2*Lneg)./(D2*sle2))
    k4 = k1*(cos(sle1*Lneg).*sin(sle2*Lneg) - D1*sle1.*cos(sle2*Lneg).*sin(sle1*Lneg)./(D2*sle2))
    k5 = k3*(cos(sle2*Lnegsep).*cos(sle3*Lnegsep) + D2*sle2.*sin(sle2*Lnegsep).*sin(sle3*Lnegsep)./(D3*sle3))+k4*(sin(sle2*Lnegsep).*cos(sle3*Lnegsep) - D2*sle2.*cos(sle2*Lnegsep).*sin(sle3*Lnegsep)./(D3*sle3))
    k6 = k3*(cos(sle2*Lnegsep).*sin(sle3*Lnegsep) - D2*sle2.*sin(sle2*Lnegsep).*cos(sle3*Lnegsep)./(D3*sle3))+k4*(sin(sle2*Lnegsep).*sin(sle3*Lnegsep) + D2*sle2.*cos(sle2*Lnegsep).*cos(sle3*Lnegsep)./(D3*sle3))
    Psiprime = -k5.*sle3.*sin(sle3*Ltot) + k6.*sle3.*cos(sle3*Ltot)
end