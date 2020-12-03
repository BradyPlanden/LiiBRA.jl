function C_e(CellData::Cell)
""" 
Electrolyte Concentration Transfer Function
# Add License
# Add Ins and Outs
    # Locations to be computed
    # Sampling Frequency
"""
Lpos = CellData.Geo.Lpos
Lneg = CellData.Geo.Lneg
Lsep = CellData.Geo.Lsep
Ltot = CellData.Geo.Ltot


F = CellData.Const.F      # Faraday Constant
R = CellData.Const.R       # Universal Gas Constant
T = CellData.Const.T      # Temperature
t_plus = CellData.Const.t_plus  # Transference Number
ζ = (1-t_plus)/F    #Simplifying Variable


Rs_Neg = CellData.Neg.Rs       # Particle radius [m]
Rs_Pos = CellData.Pos.Rs       # Particle radius [m]
Ds_Neg = CellData.Neg.Ds       # Solid diffusivity [m^2/s]
Ds_Pos = CellData.Pos.Ds       # Solid diffusivity [m^2/s]
CC_A = CellData.Geo.CC_A   # Current-collector area [m^2]
as_neg = 3*CellData.Neg.ϵ_s/Rs_Neg # Specific interfacial surf. area
as_pos = 3*CellData.Pos.ϵ_s/Rs_Pos # Specific interfacial surf. area
ϵ1 = CellData.Neg.ϵ_e      # Porosity of negative electrode
ϵ2 = CellData.Sep.ϵ_e      # Porosity of separator
ϵ3 = CellData.Pos.ϵ_e      # Porosity of positive electrode
D1 = CellData.Const.De * ϵ1^CellData.Neg.De_brug # Effective ...
D2 = CellData.Const.De * ϵ2^CellData.Sep.De_brug # diffusivities ...
D3 = CellData.Const.De * ϵ3^CellData.Pos.De_brug # of cell regions


κ_eff_Neg = CellData.Const.κ*ϵ1^CellData.Neg.κ_brug
κ_eff_Pos = CellData.Const.κ*ϵ3^CellData.Pos.κ_brug

σ_eff_Neg = CellData.Neg.σ*ϵ1^CellData.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = CellData.Pos.σ*ϵ3^CellData.Pos.σ_brug #Effective Conductivity Pos


#Defining SOC
θ_neg = CellData.Const.Init_SOC * (CellData.Neg.θ_max-CellData.Neg.θ_min) + CellData.Neg.θ_min 
θ_pos = CellData.Const.Init_SOC * (CellData.Pos.θ_max-CellData.Pos.θ_min) + CellData.Pos.θ_min 

#Beta's
βn = Rs_Neg*(s*Ds_Neg)^(1/2)
βp = Rs_Pos*(s*Ds_Pos)^(1/2)


#Prepare for j0
ce0 = CellData.Const.ce0
cs_max_neg = CellData.Neg.cs_max
cs0_neg = cs_max_neg * θ_neg

cs_max_pos = CellData.Pos.cs_max
cs0_pos = cs_max_pos * θ_pos

α_neg = CellData.Neg.α
α_pos = CellData.Pos.α

#Current Flux Density
j0_neg = CellData.Neg.k_norm*(ce0*cs_max_neg*cs0_neg)^(1-α_neg)*cs0_neg^α_neg
j0_pos = CellData.Neg.k_norm*(ce0*cs_max_pos*cs0_pos)^(1-α_pos)*cs0_pos^α_pos

#Resistances
Rct_neg = R*T/(j0_neg*F)^2
Rtot_neg = Rct_neg + CellData.Const.Rfilm_neg

Rct_pos = R*T/(j0_pos*F)^2
Rtot_pos = Rct_pos + CellData.Const.Rfilm_pos


∂Uocp_neg = UOCP(θ_neg)
∂Uocp_pos = UOCP(θ_pos)


#Create all k's
in1 = sqrt(λ*ϵ1/D1)
in2 = sqrt(λ*ϵ2/D2)
in3 = sqrt(λ*ϵ3/D3)

Bound_Neg_1 = Lneg * in1
Bound_Sep_0 = Lneg * in2
Bound_Sep_1 = (Lneg+Lsep) * in2
Bound_Pos_0 = (Lneg+Lsep) * in3
Bound_Pos_1 = Ltot * in3



#Condensing Variable
ν_n = Lneg*(as_neg/σ_eff_Neg+as_neg/κ_eff_Neg)^(1/2)/(Rtot_neg+∂Uocp_neg*(Rs_Neg/(F*Ds_Neg))*(tanh(βn)/(tanh(βn)-βn)))
ν_p = Lpos*(as_pos/σ_eff_Pos+as_pos/κ_eff_Pos)^(1/2)/(Rtot_pos+∂Uocp_pos*(Rs_Pos/(F*Ds_Pos))*(tanh(βp)/(tanh(βp)-βp)))


Lneg⋆ = Lneg * (ϵ1 * λ_k / Ds_Neg)^1/2
Lpos⋆ = Lpos * (ϵ3 * λ_k / Ds_Pos)^1/2
L⋆ = Ltot * (ϵ3 * λ_k / Ds_Pos)^1/2
Lnm⋆ = (Lneg+Lsep) * (ϵ3 * λ_k / Ds_Pos)^1/2


j_Neg = k1*ζ*Lneg⋆*sin(Lneg⋆)*(κ_eff_Neg+σ_eff_Neg*cosh(ν_n)*ν_s)/(CC_A*(κ_eff_Neg+σ_eff_Neg)*(Lneg⋆^2+ν_n^2*sinh(ν_n)))

Hlp1 = σ_eff_Pos+κ_eff_Pos
Hlp2 = (Hlp1*cosh(ν_p)*ν_p)
Hlp3 = (Lpos⋆^2 + ν_p^2)*sinh(ν_p)

j_Pos1 = (k6*ζ*Lpos⋆*cos(L⋆)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos2 = (k5*ζ*Lpos⋆*sin(Lpos⋆)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos3 = (k6*ζ*Lpos⋆*cos(Lnm⋆)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos4 = (k5*ζ*Lpos⋆*sin(L⋆)*Hlp2)/(CC_A*Hlp1*Hlp3)
j_Pos5 = (k5*ζ*σ_eff_Pos*cos(Lnm⋆)*κ_eff_Pos*cos(L⋆)*ν_p^2)/(CC_A*Hlp1*(Lpos⋆^2 + ν_p^2))
j_Pos6 = (k6*ζ*σ_eff_Pos*sin(Lnm⋆)*κ_eff_Pos*sin(L⋆)*ν_p^2)/(CC_A*Hlp1*(Lpos⋆^2 + ν_p^2))


j_Pos = j_Pos1 - j_Pos2 + j_Pos3 - j_Pos4 - j_Pos5 - j_Pos6


C_e = (1/(s+λ_k))* (j_Neg + j_Pos)

return Ltot

end

