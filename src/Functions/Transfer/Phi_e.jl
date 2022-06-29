function Phi_e(Cell,s,z,ϕ_tf,D,res0)
   """ 
   Electrolyte Potential Transfer Function

   Phi_e(CellData,s,z)

   """

   κ_eff_Neg = Cell.Const.κ*Cell.Neg.ϵ_e^Cell.Neg.κ_brug
   κ_eff_Sep = Cell.Const.κ*Cell.Sep.ϵ_e^Cell.Sep.κ_brug
   κ_eff_Pos = Cell.Const.κ*Cell.Pos.ϵ_e^Cell.Pos.κ_brug
   σ_eff_Neg = Cell.Neg.σ*Cell.Neg.ϵ_s^Cell.Neg.σ_brug #Effective Conductivity Neg
   σ_eff_Pos = Cell.Pos.σ*Cell.Pos.ϵ_s^Cell.Pos.σ_brug #Effective Conductivity Pos

   #Defining SOC
   θ_neg = Cell.Const.SOC * (Cell.Neg.θ_100-Cell.Neg.θ_0) + Cell.Neg.θ_0 
   θ_pos = Cell.Const.SOC * (Cell.Pos.θ_100-Cell.Pos.θ_0) + Cell.Pos.θ_0

   #Prepare for j0
   cs0_neg = Cell.Neg.cs_max * θ_neg
   cs0_pos = Cell.Pos.cs_max * θ_pos

   #Current Flux Density
   if Cell.Const.CellTyp == "Doyle_94"
      κ_neg = Cell.Neg.k_norm/Cell.Neg.cs_max/Cell.Const.ce0^(1-Cell.Neg.α)
      κ_pos = Cell.Pos.k_norm/Cell.Pos.cs_max/Cell.Const.ce0^(1-Cell.Pos.α)
      j0_neg = κ_neg*(Cell.Const.ce0*(Cell.Neg.cs_max-cs0_neg))^(1-Cell.Neg.α)*cs0_neg^Cell.Neg.α
      j0_pos = κ_pos*(Cell.Const.ce0*(Cell.Pos.cs_max-cs0_pos))^(1-Cell.Pos.α)*cs0_pos^Cell.Pos.α
   else
      j0_neg = Cell.Neg.k_norm*(Cell.Const.ce0*cs0_neg*(Cell.Neg.cs_max-cs0_neg))^(1-Cell.Neg.α)
      j0_pos = Cell.Pos.k_norm*(Cell.Const.ce0*cs0_pos*(Cell.Pos.cs_max-cs0_pos))^(1-Cell.Pos.α)
   end

   #Resistances 
   Rtot_neg = R*Cell.Const.T/(j0_neg*F^2) + Cell.Neg.RFilm
   Rtot_pos = R*Cell.Const.T/(j0_pos*F^2) + Cell.Pos.RFilm

   #OCP derivative
   ∂Uocp_neg = Cell.Const.∂Uocp("Neg",θ_neg)/Cell.Neg.cs_max
   ∂Uocp_pos = Cell.Const.∂Uocp("Pos",θ_pos)/Cell.Pos.cs_max

   ν_n = @. Cell.Neg.L*sqrt((Cell.Neg.as/σ_eff_Neg+Cell.Neg.as/κ_eff_Neg)/(Rtot_neg+∂Uocp_neg*(Cell.Neg.Rs/(F*Cell.Neg.Ds))*(tanh(Cell.Neg.β)/(tanh(Cell.Neg.β)-Cell.Neg.β)))) #Condensing Variable - eq. 4.13
   ν_n_∞ = @. Cell.Neg.L*sqrt(Cell.Neg.as*((1/κ_eff_Neg)+(1/σ_eff_Neg))/(Rtot_neg))
   ν_p = @. Cell.Pos.L*sqrt((Cell.Pos.as/σ_eff_Pos+Cell.Pos.as/κ_eff_Pos)/(Rtot_pos+∂Uocp_pos*(Cell.Pos.Rs/(F*Cell.Pos.Ds ))*(tanh(Cell.Pos.β)/(tanh(Cell.Pos.β)-Cell.Pos.β)))) #Condensing Variable - eq. 4.13
   ν_p_∞ = @. Cell.Pos.L*sqrt(Cell.Pos.as*((1/κ_eff_Pos)+(1/σ_eff_Pos))/(Rtot_pos))


   i=Int64(1)
   # Loop Tf's
   for pt in z
      if pt <= Cell.Neg.L+eps()
         ϕ_tf[i,:] = @. (Cell.Neg.L*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_n*pt/Cell.Neg.L)) - pt*ν_n*sinh(ν_n) + Cell.Neg.L*(cosh(ν_n)-cosh(ν_n*(Cell.Neg.L-pt)/Cell.Neg.L)))/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_n)*ν_n) #Lee. Eqn. 4.22
         zero_tf = @. -(pt^2)/(2*Cell.Const.CC_A*κ_eff_Neg*Cell.Neg.L)
         D[i]  =  @. (Cell.Neg.L*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_n_∞*pt/Cell.Neg.L)) - pt*ν_n_∞*sinh(ν_n_∞) + Cell.Neg.L*(cosh(ν_n_∞)-cosh(ν_n_∞*(Cell.Neg.L-pt)/Cell.Neg.L)))/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_n_∞)*ν_n_∞) #Lee. Eqn. 4.22 @ ∞
         ϕ_tf[i,findall(s.==0)] .= zero_tf
         
      elseif pt <= Cell.Neg.L + Cell.Sep.L + eps()
         ϕ_tf[i,:] = @. (Cell.Neg.L - pt)/(Cell.Const.CC_A*κ_eff_Sep) + (Cell.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_n/2)-ν_n))/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_n) #Lee. Eqn. 4.23
         zero_tf = @. (2*κ_eff_Neg*Cell.Neg.L-κ_eff_Sep*Cell.Neg.L-2*κ_eff_Neg*pt)/(2*Cell.Const.CC_A*κ_eff_Neg*κ_eff_Sep)
         D[i] = @. (Cell.Neg.L - pt)/(Cell.Const.CC_A*κ_eff_Sep) + (Cell.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_n_∞/2)-ν_n_∞))/(Cell.Const.CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_n_∞) #Lee. Eqn. 4.23 @ ∞
         ϕ_tf[i,findall(s.==0)] .= zero_tf
         
      else
         ϕ_tf[i,:] = @. -Cell.Sep.L/(Cell.Const.CC_A*κ_eff_Sep) + Cell.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_n/2)-ν_n)/(Cell.Const.CC_A*(σ_eff_Neg+ κ_eff_Neg)*ν_n) + (Cell.Pos.L*(-σ_eff_Pos*cosh(ν_p) + σ_eff_Pos*cosh((Cell.Const.Ltot-pt)*ν_p/Cell.Pos.L) +  κ_eff_Pos*(cosh((pt-Cell.Neg.L-Cell.Sep.L)*ν_p/Cell.Pos.L)-1)) - (pt-Cell.Neg.L-Cell.Sep.L)*κ_eff_Pos*sinh(ν_p)*ν_p)/(Cell.Const.CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_p*sinh(ν_p))
         zero_tf = @. - Cell.Neg.L/(2*Cell.Const.CC_A*κ_eff_Neg) - Cell.Sep.L/(Cell.Const.CC_A*κ_eff_Sep) - (Cell.Pos.L-(Cell.Const.Ltot-pt)^2/Cell.Pos.L)/(2*Cell.Const.CC_A*κ_eff_Pos) 
         D[i] = @. -Cell.Sep.L/(Cell.Const.CC_A*κ_eff_Sep) + Cell.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_n_∞/2)-ν_n_∞)/(Cell.Const.CC_A*(σ_eff_Neg+ κ_eff_Neg)*ν_n_∞) + (Cell.Pos.L*(-σ_eff_Pos*cosh(ν_p_∞) + σ_eff_Pos*cosh((Cell.Const.Ltot-pt)*ν_p_∞/Cell.Pos.L) +  κ_eff_Pos*(cosh((pt-Cell.Neg.L-Cell.Sep.L)*ν_p_∞/Cell.Pos.L)-1)) - (pt-Cell.Neg.L-Cell.Sep.L)*κ_eff_Pos*sinh(ν_p_∞)*ν_p_∞)/(Cell.Const.CC_A*κ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_p_∞*sinh(ν_p_∞))
         ϕ_tf[i,findall(s.==0)] .= zero_tf
         end
      i+=1
   end

   res0 .= zeros(length(z))

end
