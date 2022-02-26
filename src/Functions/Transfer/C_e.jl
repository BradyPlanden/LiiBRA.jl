@inline function C_e(Cell,s,z,tf,D,res0)
""" 
Electrolyte Concentration Transfer Function

C_e(Cell,s,z)

"""

ζ = (1-Cell.Const.t_plus)/F    #Simplifying Variable
κ_eff_Neg = Cell.Const.κ*Cell.Neg.ϵ_e^Cell.Neg.κ_brug
κ_eff_Pos = Cell.Const.κ*Cell.Pos.ϵ_e^Cell.Pos.κ_brug
σ_eff_Neg = Cell.Neg.σ*Cell.Neg.ϵ_s^Cell.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = Cell.Pos.σ*Cell.Pos.ϵ_s^Cell.Pos.σ_brug #Effective Conductivity Pos


#Defining SOC
θ_neg = Cell.Const.SOC * (Cell.Neg.θ_100-Cell.Neg.θ_0) + Cell.Neg.θ_0 
θ_pos = Cell.Const.SOC * (Cell.Pos.θ_100-Cell.Pos.θ_0) + Cell.Pos.θ_0 

#Prepare for j0
ce0 = Cell.Const.ce0
cs_max_neg = Cell.Neg.cs_max
cs0_neg = cs_max_neg * θ_neg
cs_max_pos = Cell.Pos.cs_max
cs0_pos = cs_max_pos * θ_pos

#Current Flux Density

if Cell.Const.CellTyp == "Doyle_94"
    κ_pos = Cell.Pos.k_norm/Cell.Pos.cs_max/ce0^(1-Cell.Pos.α)
    κ_neg = Cell.Neg.k_norm/Cell.Neg.cs_max/ce0^(1-Cell.Neg.α)
    j0_neg = κ_neg*(ce0*(cs_max_neg-cs0_neg))^(1-Cell.Neg.α)*cs0_neg^Cell.Neg.α
    j0_pos = κ_pos*(ce0*(cs_max_pos-cs0_pos))^(1-Cell.Pos.α)*cs0_pos^Cell.Pos.α
 else
    j0_neg = Cell.Neg.k_norm*(Cell.Const.ce0*cs0_neg*(Cell.Neg.cs_max-cs0_neg))^(1-Cell.Neg.α)
    j0_pos = Cell.Pos.k_norm*(Cell.Const.ce0*cs0_pos*(Cell.Pos.cs_max-cs0_pos))^(1-Cell.Pos.α)
 end

#Resistances
Rtot_neg = R*Cell.Const.T/(j0_neg*F^2) + Cell.Neg.RFilm
Rtot_pos = R*Cell.Const.T/(j0_pos*F^2) + Cell.Pos.RFilm
#Rtot_neg = R*Cell.Const.T/(j0_neg*Cell.Const.CC_A*F) + Cell.Neg.RFilm
#Rtot_pos = R*Cell.Const.T/(j0_pos*Cell.Const.CC_A*F) + Cell.Pos.RFilm

#OCP derivative
∂Uocp_pos = Cell.Const.∂Uocp("Pos",θ_pos)/cs_max_pos
∂Uocp_neg = Cell.Const.∂Uocp("Neg",θ_neg)/cs_max_neg

#Condensing Variable
ν_n =  @. Cell.Neg.L*sqrt((Cell.Neg.as/σ_eff_Neg+Cell.Neg.as/κ_eff_Neg)/(Rtot_neg+∂Uocp_neg*(Cell.Neg.Rs/(F*Cell.Neg.Ds))*(tanh(Cell.Neg.β)/(tanh(Cell.Neg.β)-Cell.Neg.β))))
ν_p =  @. Cell.Pos.L*sqrt((Cell.Pos.as/σ_eff_Pos+Cell.Pos.as/κ_eff_Pos)/(Rtot_pos+∂Uocp_pos*(Cell.Pos.Rs/(F*Cell.Pos.Ds))*(tanh(Cell.Pos.β)/(tanh(Cell.Pos.β)-Cell.Pos.β))))


R_ce = find_zeros(x->flambda(Cell,x),0.00,Cell.Const.CeRootRange)
if size(R_ce,1) >= Cell.Const.Ce_M+1
    λ = R_ce[2:Cell.Const.Ce_M+1]
else
    throw(DomainError(R_ce, "Ce roots 'R_ce' does not contain enough values, increase Cell.Const.CeRootRange"))
end


#Create all k's
in1 = @. sqrt(λ*Cell.Neg.ϵ_e/Cell.Const.D1)
in2 = @. sqrt(λ*Cell.Sep.ϵ_e/Cell.Const.D2)
in3 = @. sqrt(λ*Cell.Pos.ϵ_e/Cell.Const.D3)

Bound_Neg_1 = Cell.Neg.L * in1
Bound_Sep_0 = Cell.Neg.L * in2
Bound_Sep_1 = (Cell.Neg.L+Cell.Sep.L) * in2
Bound_Pos_0 = (Cell.Neg.L+Cell.Sep.L) * in3
Bound_Pos_1 = (Cell.Neg.L+Cell.Sep.L+Cell.Pos.L) * in3
Bound_Pos_2 = Cell.Pos.L * in3

#Scaled coefficients
k3_s = @. cos(Bound_Neg_1)*cos(Bound_Sep_0) + Cell.Const.D1*in1*sin(Bound_Neg_1)*sin(Bound_Sep_0)/(Cell.Const.D2*in2)
k4_s = @. cos(Bound_Neg_1)*sin(Bound_Sep_0) - Cell.Const.D1*in1*cos(Bound_Sep_0)*sin(Bound_Neg_1)/(Cell.Const.D2*in2)
k5_s = @. k3_s*(cos(Bound_Sep_1)*cos(Bound_Pos_0)+Cell.Const.D2*in2*sin(Bound_Sep_1)*sin(Bound_Pos_0)/(Cell.Const.D3*in3))+k4_s*(sin(Bound_Sep_1)*cos(Bound_Pos_0)-Cell.Const.D2*in2*cos(Bound_Sep_1)*sin(Bound_Pos_0)/(Cell.Const.D3*in3))
k6_s = @. k3_s*(cos(Bound_Sep_1)*sin(Bound_Pos_0)-Cell.Const.D2*in2*sin(Bound_Sep_1)*cos(Bound_Pos_0)/(Cell.Const.D3*in3))+k4_s*(sin(Bound_Sep_1)*sin(Bound_Pos_0)+Cell.Const.D2*in2*cos(Bound_Sep_1)*cos(Bound_Pos_0)/(Cell.Const.D3*in3))

#Solving for k1:
Int_ψ1 = @. Cell.Neg.ϵ_e*(2*Bound_Neg_1+sin(2*Bound_Neg_1))/(4*in1)
Int_ψ2 = @. Cell.Sep.ϵ_e/(4*in2)*(2*(k3_s^2+k4_s^2)*Cell.Sep.L*in2+2*k3_s*k4_s*cos(2*Bound_Sep_0)-2*k3_s*k4_s*cos(2*Bound_Sep_1)-(k3_s-k4_s)*(k3_s+k4_s)*(sin(2*Bound_Sep_0)-sin(2*Bound_Sep_1)))
Int_ψ3 = @. Cell.Pos.ϵ_e/(4*in3)*(2*(k5_s^2+k6_s^2)*Cell.Pos.L*in3+2*k5_s*k6_s*cos(2*Bound_Pos_0)-2*k5_s*k6_s*cos(2*Bound_Pos_1)-(k5_s-k6_s)*(k5_s+k6_s)*(sin(2*Bound_Pos_0)-sin(2*Bound_Pos_1)))

k1 = @. 1/(sqrt(Int_ψ1+Int_ψ2+Int_ψ3))
k3 = @. k1*k3_s
k4 = @. k1*k4_s
k5 = @. k1*k5_s
k6 = @. k1*k6_s

j_Neg = @. k1*ζ*ν_n*(Bound_Neg_1*(κ_eff_Neg+σ_eff_Neg*cosh(ν_n))*sin(Bound_Neg_1)+(κ_eff_Neg+σ_eff_Neg*cos(Bound_Neg_1))*sinh(ν_n)*ν_n)/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*(Bound_Neg_1^2+ν_n^2)*sinh(ν_n))
zero_tf_neg = @. k1*ζ*sin(Bound_Neg_1)/(Cell.Const.CC_A*Bound_Neg_1)
j_Neg[:,findall(s.==0)] .= zero_tf_neg[:,findall(s.==0)]

j_Pos = @. -ζ*ν_p/(Cell.Const.CC_A*(κ_eff_Pos+σ_eff_Pos)*(Bound_Pos_2^2+ν_p^2)*sinh(ν_p))*(-k6*Bound_Pos_2*cos(Bound_Pos_1)*(σ_eff_Pos+κ_eff_Pos*cosh(ν_p))+Bound_Pos_2*(κ_eff_Pos+σ_eff_Pos*cosh(ν_p))*(k6*cos(Bound_Pos_0)-k5*sin(Bound_Pos_0))+k5*Bound_Pos_2*(σ_eff_Pos+ κ_eff_Pos*cosh(ν_p))*sin(Bound_Pos_1)+sinh(ν_p)*(k5*σ_eff_Pos*cos(Bound_Pos_0)+k5*κ_eff_Pos*cos(Bound_Pos_1)+k6*σ_eff_Pos*sin(Bound_Pos_0)+k6*κ_eff_Pos*sin(Bound_Pos_1))*ν_p)
zero_tf = @. -ζ*(k6*(cos(Bound_Pos_0)-cos(Bound_Pos_1))+k5*(sin(Bound_Pos_1)-sin(Bound_Pos_0)))/(Cell.Const.CC_A*Bound_Pos_2)
j_Pos[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]

 
i=Int64(1)
ψ = fill(0.0,length(z),length(λ))
for loop in 1:length(λ)
    i=1
    for x in z #Eigen Weighting
        if x < Cell.Neg.L+eps()
            ψ[i,loop] = k1[loop]*cos(in1[loop]*x) #negative electrode
        elseif x > (Cell.Neg.L+Cell.Sep.L)-eps()
            ψ[i,loop] = k5[loop]*cos(in3[loop]*x)+k6[loop]*sin(in3[loop]*x) # postive electrode
        else
            ψ[i,loop] = k3[loop]*cos(in2[loop]*x)+k4[loop]*sin(in2[loop]*x) # separator
        end
    i = i+1
    end
end

tf .= ψ*((j_Neg .+ j_Pos)./(s.+λ))
D .=  zeros(length(z))
res0 .= zeros(length(z))

end

function flambda(Cell,λ)
    ϵ1 = Cell.Neg.ϵ_e
    ϵ2 = Cell.Sep.ϵ_e
    ϵ3 = Cell.Pos.ϵ_e
    D1 = Cell.Const.D1
    D2 = Cell.Const.D2
    D3 = Cell.Const.D3
    Lneg = Cell.Neg.L
    Lnegsep = Cell.Const.Lnegsep
    Ltot = Cell.Neg.L+Cell.Sep.L+Cell.Pos.L
    k1 = Float64(1)
    sle1 = sqrt(λ*ϵ1/D1)
    sle2 = sqrt(λ*ϵ2/D2)
    sle3 = sqrt(λ*ϵ3/D3)
    k3 = @. k1*(cos(sle1*Lneg)*cos(sle2*Lneg) + D1*sle1*sin(sle1*Lneg)*sin(sle2*Lneg)/(D2*sle2))
    k4 = @. k1*(cos(sle1*Lneg)*sin(sle2*Lneg) - D1*sle1*cos(sle2*Lneg)*sin(sle1*Lneg)/(D2*sle2))
    k5 = @. k3*(cos(sle2*Lnegsep)*cos(sle3*Lnegsep) + D2*sle2*sin(sle2*Lnegsep)*sin(sle3*Lnegsep)/(D3*sle3))+k4*(sin(sle2*Lnegsep)*cos(sle3*Lnegsep) - D2*sle2*cos(sle2*Lnegsep)*sin(sle3*Lnegsep)/(D3*sle3))
    k6 = @. k3*(cos(sle2*Lnegsep)*sin(sle3*Lnegsep) - D2*sle2*sin(sle2*Lnegsep)*cos(sle3*Lnegsep)/(D3*sle3))+k4*(sin(sle2*Lnegsep)*sin(sle3*Lnegsep) + D2*sle2*cos(sle2*Lnegsep)*cos(sle3*Lnegsep)/(D3*sle3))
    Psiprime = @. -k5*sle3*sin(sle3*Ltot) + k6*sle3*cos(sle3*Ltot)
    if λ == 0
        Psiprime = 0
    end
    return Psiprime
end