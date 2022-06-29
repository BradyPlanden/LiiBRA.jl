function Phi_se(Cell,s,z,Def,ϕ_tf,D,res0)
   """ 
   Solid-Electrolyte Potential Transfer Function

   Phi_se(Cell,s,z,Def)

   """

   if Def == "Pos"   
      Electrode = Cell.Pos #Electrode Length
   else
      Electrode = Cell.Neg #Electrode Length
   end

   κ_eff = Cell.Const.κ*Electrode.ϵ_e^Electrode.κ_brug #Effective Electrolyte Conductivity 
   σ_eff = Electrode.σ*Electrode.ϵ_s^Electrode.σ_brug #Effective Electrode Conductivity 

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
   ∂Uocp_elc = Cell.Const.∂Uocp(Def,θ)/Electrode.cs_max #Open Circuit Potential Partial
   res0 .= @. -3*∂Uocp_elc/(Electrode.as*F*Electrode.L*Cell.Const.CC_A*Electrode.Rs) # residual for pole removal

   ν = @. Electrode.L*sqrt((Electrode.as/σ_eff+Electrode.as/κ_eff)/(Rtot+∂Uocp_elc*(Electrode.Rs/(F*Electrode.Ds))*(tanh(Electrode.β)/(tanh(Electrode.β)-Electrode.β)))) #Condensing Variable - eq. 4.13
   ν_∞ = @. Electrode.L*sqrt(Electrode.as*((1/κ_eff)+(1/σ_eff))/(Rtot))

   ϕ_tf .= @. Electrode.L/(Cell.Const.CC_A*ν*sinh(ν))*((1/κ_eff)*cosh(ν*z)+(1/σ_eff)*cosh(ν*(z-1)))-res0/s #Transfer Function - eq. 4.14

   zero_tf = @. (6*(5*Electrode.Ds*F*Rtot-∂Uocp_elc*Electrode.Rs)*σ_eff)/(30*Cell.Const.CC_A*Electrode.as*Electrode.Ds*F*σ_eff*Electrode.L) + (5*Electrode.as*Electrode.Ds*F*Electrode.L^2*(σ_eff*(-1+3*z^2)+κ_eff*(2-6*z+3*z^2)))/(30*Cell.Const.CC_A*Electrode.as*Electrode.Ds*F*σ_eff*κ_eff*Electrode.L)
   D .= @. Electrode.L/(Cell.Const.CC_A*ν_∞*sinh(ν_∞))*((1/κ_eff)*cosh(ν_∞*z)+(1/σ_eff)*cosh(ν_∞*(z-1))) # Contribution to D as G->∞
   ϕ_tf[:,findall(s.==0)] .= zero_tf[:,findall(s.==0)]
   res0 .= zeros(length(z))

   if Def == "Pos" #Double check this implementation
      ϕ_tf .= -ϕ_tf
      D .= -D
   end

end
