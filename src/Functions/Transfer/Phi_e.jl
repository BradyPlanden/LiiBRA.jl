@inline function Phi_e(CellData,s,z)
   """ 
   Electrolyte Potential Transfer Function
   # Add License
   # Add Ins and Outs
       # Cell Data 
       # Frequency Vector 
       # Discretisation Locations
       # Electrode Definition
   """


CC_A = CellData.Const.CC_A   # Current-collector area [m^2]
De = CellData.Const.De # Electrolyte Diffusivity
κ_eff_Neg = CellData.Const.κ*CellData.Neg.ϵ_e^CellData.Neg.κ_brug
κ_eff_Sep = CellData.Const.κ*CellData.Sep.ϵ_e^CellData.Sep.κ_brug
κ_eff_Pos = CellData.Const.κ*CellData.Pos.ϵ_e^CellData.Pos.κ_brug
σ_eff_Neg = CellData.Neg.σ*CellData.Neg.ϵ_s^CellData.Neg.σ_brug #Effective Conductivity Neg
σ_eff_Pos = CellData.Pos.σ*CellData.Pos.ϵ_s^CellData.Pos.σ_brug #Effective Conductivity Pos
#dln = CellData.Const.dln  #Electrolyte activity coefficient term (Rod. 17)
#κ_D_eff = (2*R*CellData.Const.T/F)*κ_eff_Neg*(1-CellData.Const.t_plus)*(1+dln) #Diffision Effective Electrolyte Conductivity

#Defining SOC
θ_neg = CellData.Const.SOC * (CellData.Neg.θ_100-CellData.Neg.θ_0) + CellData.Neg.θ_0 
θ_pos = CellData.Const.SOC * (CellData.Pos.θ_100-CellData.Pos.θ_0) + CellData.Pos.θ_0

#Beta's
βn = @. CellData.Neg.Rs*sqrt(s/CellData.Neg.Ds)
βp = @. CellData.Pos.Rs*sqrt(s/CellData.Pos.Ds)

#Prepare for j0
ce0 = CellData.Const.ce0
cs_max_neg = CellData.Neg.cs_max
cs0_neg = cs_max_neg * θ_neg
cs_max_pos = CellData.Pos.cs_max
cs0_pos = cs_max_pos * θ_pos
α_neg = CellData.Neg.α
α_pos = CellData.Pos.α

#Current Flux Density
if CellData.Const.CellTyp == "Doyle_94"
   κ_neg = CellData.Neg.k_norm/CellData.Neg.cs_max/ce0^(1-α_neg)
   κ_pos = CellData.Pos.k_norm/CellData.Pos.cs_max/ce0^(1-α_pos)
   j0_neg = κ_neg*(ce0*(cs_max_neg-cs0_neg))^(1-α_neg)*cs0_neg^α_neg
   j0_pos = κ_pos*(ce0*(cs_max_pos-cs0_pos))^(1-α_pos)*cs0_pos^α_pos
else
   j0_neg = CellData.Neg.k_norm*(ce0*(cs0_neg/cs_max_neg*(1-cs0_neg/cs_max_neg)))^(1-CellData.Neg.α)
   j0_pos = CellData.Pos.k_norm*(ce0*(cs0_pos/cs_max_pos*(1-cs0_pos/cs_max_pos)))^(1-CellData.Pos.α)
end

#Resistances 
Rtot_neg = R*CellData.Const.T/(j0_neg*F^2) + CellData.Neg.RFilm
Rtot_pos = R*CellData.Const.T/(j0_pos*F^2) + CellData.Pos.RFilm

#OCP derivative
∂Uocp_neg = CellData.Const.∂Uocp("Neg",θ_neg)/cs_max_neg
∂Uocp_pos = CellData.Const.∂Uocp("Pos",θ_pos)/cs_max_pos

ν_neg = @. CellData.Neg.L*sqrt((CellData.Neg.as/σ_eff_Neg+CellData.Neg.as/κ_eff_Neg)/(Rtot_neg+∂Uocp_neg*(CellData.Neg.Rs/(F*CellData.Neg.Ds))*(tanh(βn)/(tanh(βn)-βn)))) #Condensing Variable - eq. 4.13
ν_neg_∞ = @. CellData.Neg.L*sqrt(CellData.Neg.as*((1/κ_eff_Neg)+(1/σ_eff_Neg))/(Rtot_neg))
ν_pos = @. CellData.Pos.L*sqrt((CellData.Pos.as/σ_eff_Pos+CellData.Pos.as/κ_eff_Pos)/(Rtot_pos+∂Uocp_pos*(CellData.Pos.Rs/(F*CellData.Pos.Ds ))*(tanh(βp)/(tanh(βp)-βp)))) #Condensing Variable - eq. 4.13
ν_pos_∞ = @. CellData.Pos.L*sqrt(CellData.Pos.as*((1/κ_eff_Pos)+(1/σ_eff_Pos))/(Rtot_pos))


# if tf = Rod17 #Consider splitting into different function calls.

#    for i < length(z)
#    pt = z[i]

#       if pt <= CellData.Neg.L
#          ϕ_tf = @. CellData.Neg.L*(σ_eff*(1-cosh(ν_neg*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (CellData.Neg.L*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-pt))-pt*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
#          zero_tf = @. -(CellData.Neg.L*pt^2)/(2*as*κ_eff) - (κ_d*CellData.Neg.L)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
#          #D  =  @. CellData.Neg.L*(σ_eff*(1-cosh(ν_∞*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (CellData.Neg.L*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞)
#          ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

#       elseif pt <= CellData.Neg.L + CellData.Sep.L
#          ϕ_tf1 = @. -pt*CellData.Sep.L/(as*κ_eff)
#          ϕ_tf2 = @. -(κ_D_eff/κ_eff_Pos*ce0)*(cS1*(e^ΛS1*pt-1)+cS2*(e^ΛS1*pt-1))
#          ϕ_n_tf = @. CellData.Neg.L*(σ_eff*(1-cosh(ν_neg*CellData.Neg.L)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (CellData.Neg.L*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-CellData.Neg.L))-CellData.Neg.L*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
#          ϕ_tf = @. ϕ_n_tf + ϕ_tf1 + ϕ_tf2 
#          zero_tf = @. -(CellData.Neg.L*pt^2)/(2*as*κ_eff) - (κ_d*CellData.Neg.L)(t_plus-1)*pt^2)/(2*A*ce0*De*F*κ_eff)
#          #D  = @. CellData.Sep.L*(σ_eff*(1-cosh(ν_∞*pt)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) + (CellData.Sep.L*(κ_eff*cosh(ν_∞)-cosh(ν_∞*(1-pt))-pt*sinh(ν_∞)*ν_∞))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_∞)*ν_∞) - pt*CellData.Sep.L/(as*κ_eff)
#          ϕ_tf[findall(s.==0),:] .= zero_tf[findall(s.==0),:]

#       else
#          ϕ_s_tf1 = @. -pt*CellData.Sep.L/(as*κ_eff)
#          ϕ_s_tf2 = @. -(κ_D_eff/κ_eff_Pos*ce0)*(cS1*(e^ΛS1*pt-1)+cS2*(e^ΛS1*pt-1))
#          ϕ_n_tf = @. CellData.Neg.L*(σ_eff*(1-cosh(ν_neg*CellData.Neg.L)/(as*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg) + (CellData.Neg.L*(κ_eff*cosh(ν_neg)-cosh(ν_neg*(1-CellData.Neg.L))-CellData.Neg.L*sinh(ν_neg)*ν_neg))/(A*κ_eff*(κ_eff+σ_eff)*sinh(ν_neg)*ν_neg)
#          ϕ_p_tf1 = ((pt-1)*CellData.Pos.L/κ_eff_Pos*CC_A)+(CellData.Pos.as*F*CellData.Pos.L^2/κ_eff_Pos) * (j_Pos_1*(e^(ΛP1*pt)+(1-pt)*ΛP1*e^ΛP1-e^ΛP1)/ΛP1^2 + j_Pos_2*(e^(-ΛP1*pt)+(pt-1)*ΛP1*e^-ΛP1-e^-ΛP1)/ΛP1^2 + j_Pos_1*(e^(ΛP1*pt)+(1-pt)*ΛP1*e^ΛP1-e^ΛP1)/ΛP1^2 + j_Pos_3*(e^(ΛP2*pt)+(1-pt)*ΛP2*e^ΛP2-e^ΛP2)/ΛP2^2 + j_Pos_4*(e^(-ΛP2*pt)+(pt-1)*ΛP2*e^-ΛP2-e^-ΛP2)/ΛP2^2)
#          ϕ_p_tf2 = -(κ_D_eff/κ_eff_Pos*ce0)*(cp1*(e^ΛP1*pt-e^Λp1)+cp2*(e^ΛP1*pt-e^Λp1)+cp3*(e^ΛP2*pt-e^Λp2)+cp4*(e^ΛP2*pt-e^Λp2))
#          ϕ_tf = @. ϕ_n_tf + ϕ_s_tf1 + ϕ_s_tf2 + ϕ_p_tf1 + ϕ_p_tf2

#          #D  = @.  
#    end

ϕ_tf = Array{ComplexF64}(undef,length(z),length(s))
D_term = Array{String}(undef,length(z),1)
D = fill(0.0,length(z))
i=Int64(1)
# Loop Tf's
   for pt in z
         if pt <= CellData.Neg.L+eps()
            ϕ_tf[i,:] = @. (CellData.Neg.L*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg*pt/CellData.Neg.L)) - pt*ν_neg*sinh(ν_neg) + CellData.Neg.L*(cosh(ν_neg)-cosh(ν_neg*(CellData.Neg.L-pt)/CellData.Neg.L)))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg) #Lee. Eqn. 4.22
            zero_tf = @. -(pt^2)/(2*CC_A*κ_eff_Neg*CellData.Neg.L)
            D[i]  =  @. (CellData.Neg.L*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg_∞*pt/CellData.Neg.L)) - pt*ν_neg_∞*sinh(ν_neg_∞) + CellData.Neg.L*(cosh(ν_neg_∞)-cosh(ν_neg_∞*(CellData.Neg.L-pt)/CellData.Neg.L)))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg_∞)*ν_neg_∞) #Lee. Eqn. 4.22 @ ∞
            ϕ_tf[i,findall(s.==0)] .= zero_tf
            D_term[i] = "@. ($(CellData.Neg.L)*($σ_eff_Neg/$κ_eff_Neg)*(1-cosh(ν_neg*$pt/$(CellData.Neg.L))) - $pt*ν_neg*sinh(ν_neg))/($CC_A*($κ_eff_Neg+$σ_eff_Neg)*sinh(ν_neg)*ν_neg) + ($(CellData.Neg.L)*(cosh(ν_neg)-cosh(ν_neg*($(CellData.Neg.L)-$pt)/$(CellData.Neg.L))/($CC_A*$κ_eff_Neg*($κ_eff_Neg+$σ_eff_Neg)*sinh(ν_neg)*ν_neg)))" #Lee. Eqn. 4.22 @ ∞
         
         elseif pt <= CellData.Neg.L + CellData.Sep.L + eps()
            ϕ_tf[i,:] = @. (CellData.Neg.L - pt)/(CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) #Lee. Eqn. 4.23
            zero_tf = @. (2*κ_eff_Neg*CellData.Neg.L-κ_eff_Sep*CellData.Neg.L-2*κ_eff_Neg*pt)/(2*CC_A*κ_eff_Neg*κ_eff_Sep)
            D[i] = @. (CellData.Neg.L - pt)/(CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg_∞/2)-ν_neg_∞))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg_∞) #Lee. Eqn. 4.23 @ ∞
            ϕ_tf[i,findall(s.==0)] .= zero_tf
            D_term[i] = "@. ($(CellData.Neg.L - pt))/($CC_A*$κ_eff_Sep) + ($(CellData.Neg.L)*((1-$σ_eff_Neg/$κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/($CC_A*($κ_eff_Neg+$σ_eff_Neg)*ν_neg)" #Lee. Eqn. 4.23 @ ∞
   
         else
            #ϕ_tf[i,:] = @. -CellData.Sep.L/(CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) - (CellData.Pos.L*(1+(σ_eff_Pos/κ_eff_Pos)*cosh(ν_pos)) + CellData.Pos.L*cosh((CellData.Neg.L+CellData.Sep.L-pt)*ν_pos/CellData.Pos.L) + CellData.Pos.L*((σ_eff_Pos/κ_eff_Pos)*cosh((Ltot-pt)*ν_pos/CellData.Pos.L)))/(CC_A*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos)*ν_pos) + (CellData.Neg.L + CellData.Sep.L - pt)/(CC_A*((κ_eff_Pos+σ_eff_Pos)))  #Lee. Eqn. 4.24
            #ϕ_tf[i,:] = @. (-CellData.Neg.L*((σ_eff_Neg-κ_eff_Neg)*tanh(ν_neg/2)))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) - CellData.Neg.L/(CC_A*(κ_eff_Neg+σ_eff_Neg)) - (CellData.Sep.L/(CC_A*κ_eff_Sep)) - CellData.Pos.L*(1-cosh((pt-1)*ν_pos))/((CC_A*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos)*ν_pos)) - CellData.Pos.L*σ_eff_Pos*(cosh(ν_pos)-cosh(pt*ν_pos))/((CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos)*ν_pos)) - CellData.Pos.L*(pt-1)/(CC_A*(κ_eff_Pos+σ_eff_Pos))
            ϕ_tf[i,:] = @. -CellData.Sep.L/(CC_A*κ_eff_Sep) + CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg)/(CC_A*(σ_eff_Neg+ κ_eff_Neg)*ν_neg) + (CellData.Pos.L*(-σ_eff_Pos*cosh(ν_pos) + σ_eff_Pos*cosh((Ltot-pt)*ν_pos/CellData.Pos.L) +  κ_eff_Pos*(cosh((pt-CellData.Neg.L-CellData.Sep.L)*ν_pos/CellData.Pos.L)-1)) - (pt-CellData.Neg.L-CellData.Sep.L)*κ_eff_Pos*sinh(ν_pos)*ν_pos)/(CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_pos*sinh(ν_pos))
            #zero_tf = @. -CellData.Neg.L/(2*CC_A*κ_eff_Neg) - CellData.Sep.L/(2*CC_A*κ_eff_Sep) + (CellData.Neg.L+CellData.Sep.L-pt)*(CellData.Neg.L+2*CellData.Pos.L+CellData.Sep.L-pt)/(2*CC_A*κ_eff_Pos*CellData.Pos.L)
            zero_tf = @. - CellData.Neg.L/(2*CC_A*κ_eff_Neg) - CellData.Sep.L/(CC_A*κ_eff_Sep) - (CellData.Pos.L-(Ltot-pt)^2/CellData.Pos.L)/(2*CC_A*κ_eff_Pos) 
            D[i] = @. -CellData.Sep.L/(CC_A*κ_eff_Sep) + CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg_∞/2)-ν_neg_∞)/(CC_A*(σ_eff_Neg+ κ_eff_Neg)*ν_neg_∞) + (CellData.Pos.L*(-σ_eff_Pos*cosh(ν_pos_∞) + σ_eff_Pos*cosh((Ltot-pt)*ν_pos_∞/CellData.Pos.L) +  κ_eff_Pos*(cosh((pt-CellData.Neg.L-CellData.Sep.L)*ν_pos_∞/CellData.Pos.L)-1)) - (pt-CellData.Neg.L-CellData.Sep.L)*κ_eff_Pos*sinh(ν_pos_∞)*ν_pos_∞)/(CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_pos_∞*sinh(ν_pos_∞))
            #D[i]  = @. -CellData.Sep.L/(CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg_∞/2)-ν_neg_∞))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg_∞) - (CellData.Pos.L*(1+(σ_eff_Pos/κ_eff_Pos)*cosh(ν_pos_∞)) + CellData.Pos.L*cosh((CellData.Neg.L+CellData.Sep.L-pt)*ν_pos_∞/CellData.Pos.L) + CellData.Pos.L*((σ_eff_Pos/κ_eff_Pos)*cosh((Ltot-pt)*ν_pos_∞/CellData.Pos.L)))/(CC_A*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos_∞)*ν_pos_∞) + (CellData.Neg.L + CellData.Sep.L - pt)/(CC_A*((κ_eff_Pos+σ_eff_Pos))) #Lee. Eqn. 4.24 @ ∞
            ϕ_tf[i,findall(s.==0)] .= zero_tf
            D_term[i] = "@. -$(CellData.Sep.L)/($CC_A*$κ_eff_Sep) + ($(CellData.Neg.L)*((1-$σ_eff_Neg/$κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/($CC_A*($κ_eff_Neg+$σ_eff_Neg)*ν_neg) - ($(CellData.Pos.L)*(1+($σ_eff_Pos/$κ_eff_Pos)*cosh(ν_pos)) + $(CellData.Pos.L)*cosh(($(CellData.Neg.L)+$(CellData.Sep.L)-$pt)*ν_pos/$(CellData.Pos.L)) + $(CellData.Pos.L)*(($σ_eff_Pos/$κ_eff_Pos)*cosh(($Ltot-$pt)*ν_pos/$(CellData.Pos.L))))/($CC_A*($κ_eff_Pos+$σ_eff_Pos)*sinh(ν_pos)*ν_pos) + ($(CellData.Neg.L) + $(CellData.Sep.L) - $pt)/($CC_A*(($κ_eff_Pos+$σ_eff_Pos)))" #Lee. Eqn. 4.24 @ ∞
         end
      i=i+1
 end
 if Debug == 1
      println("D:Phi_e:",size(D))
      println("ϕ_tf:Phi_e:",ϕ_tf[:,2])
      println("D:Phi_e:Pos",D)
      println("z:Phi_e:Pos",z)
      println("ν_∞:Phi_e:Pos",ν_neg_∞)
      println("j0:Phi_e:Pos",j0_neg)
      println("Rtot:Phi_e:Pos",Rtot_neg)
      println("σ_eff:Phi_e:Pos",σ_eff_Neg)
      println("∂Uocp_elc:Phi_e:Pos",∂Uocp_neg)
      println("L:Phi_e:Pos",CellData.Neg.L)
      println("D_term:", D_term)
 end
res0 = zeros(length(z))

return ϕ_tf, D, res0, D_term

end
