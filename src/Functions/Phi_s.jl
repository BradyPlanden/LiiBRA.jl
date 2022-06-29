function Phi_s(Cell,s,z,Def,ϕ_tf,D,res0)
   """ 
   Solid Potential Transfer Function

   Phi_s(Cell,s,z,Def)

   """

   if Def == "Pos"
      Electrode = Cell.Pos #Electrode Length
   else
      Electrode = Cell.Neg #Electrode Length
   end

   κ_eff = Cell.Const.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
   σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 
   comb_cond_eff = κ_eff+σ_eff #combining into single variable

   #Defining SOC
   θ = Cell.Const.SOC * (Electrode.θ_100-Electrode.θ_0) + Electrode.θ_0

   #Prepare for j0
   cs0 = Electrode.cs_max * θ

   #Current Flux Density
   if Cell.Const.CellTyp == "Doyle_94"
      κ = Electrode.k_norm/Electrode.cs_max/Cell.Const.ce0^(1-Electrode.α)
      j0 = κ*(Cell.Const.ce0*(Electrode.cs_max-cs0))^(1-Electrode.α)*cs0^Electrode.α
   else
      j0 = Electrode.k_norm*(Cell.Const.ce0*cs0*(Electrode.cs_max-cs0))^(1-Electrode.α)
   end

   #Resistance
   Rtot = R*Cell.Const.T/(j0*F^2) + Electrode.RFilm

   #∂Uocp_Def
   ∂Uocp_elc = Cell.Const.∂Uocp(Def,θ)/Electrode.cs_max

   ν = @. Electrode.L*sqrt((Electrode.as/σ_eff+Electrode.as/κ_eff)/(Rtot+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(Electrode.β)/(tanh(Electrode.β)-Electrode.β)))) #Condensing Variable - eq. 4.13
   ν_∞ = @. Electrode.L*sqrt((Electrode.as/κ_eff+Electrode.as/σ_eff)/(Rtot))

   ϕ_tf .= @. (-Electrode.L*(κ_eff*(cosh(ν)-cosh((z-1)*ν)))-Electrode.L*(σ_eff*(1-cosh(z*ν)+z*ν*sinh(ν))))/(Cell.Const.CC_A*σ_eff*(comb_cond_eff)*ν*sinh(ν)) #Transfer Function - eq. 4.19
   D .= @. (-Electrode.L*(κ_eff*(cosh(ν_∞)-cosh((z-1)*ν_∞)))-Electrode.L*(σ_eff*(1-cosh(z*ν_∞)+z*ν_∞*sinh(ν_∞))))/(Cell.Const.CC_A*σ_eff*(comb_cond_eff)*ν_∞*sinh(ν_∞)) # Contribution to D as G->∞
   zero_tf = @. Electrode.L*(z-2)*z/(2*Cell.Const.CC_A*σ_eff)
   ϕ_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]
   res0 .= zeros(length(z))

   if Def == "Pos"
      ϕ_tf .= -ϕ_tf
      D .= -D
   end

end
