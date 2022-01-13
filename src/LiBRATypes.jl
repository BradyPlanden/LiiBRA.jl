# @with_kw struct D_Linear
#     #C_e.jl
#     D_ce::Function = z -> zeros(length(z))

#     #Phie.jl
#     D_ϕe_neg::Function = (Lneg,σ_eff_Neg,κ_eff_Neg,ν_neg,ν_pos,pt,CC_A) -> @. (Lneg*(σ_eff_Neg/κ_eff_Neg)*(1-cosh(ν_neg*pt/Lneg)) - pt*ν_neg*sinh(ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg) + (Lneg*(cosh(ν_neg)-cosh(ν_neg*(Lneg-pt)/Lneg)/(CC_A*κ_eff_Neg*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg))) #Lee. Eqn. 4.22
#     D_ϕe_sep::Function = (Lneg,σ_eff_Neg,κ_eff_Neg,ν_neg,ν_pos,pt,CC_A) -> @. (CellData.Neg.L - pt)/(CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) #Lee. Eqn. 4.23
#     D_ϕe_pos::Function = (Lneg,σ_eff_Neg,κ_eff_Neg,ν_neg,ν_pos,pt,CC_A) -> @. -CellData.Sep.L/(CC_A*κ_eff_Sep) + (CellData.Neg.L*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg))/(CC_A*(κ_eff_Neg+σ_eff_Neg)*ν_neg) - (CellData.Pos.L*(1+(σ_eff_Pos/κ_eff_Pos)*cosh(ν_pos)) + CellData.Pos.L*cosh((CellData.Neg.L+CellData.Sep.L-pt)*ν_pos/CellData.Pos.L) + CellData.Pos.L*((σ_eff_Pos/κ_eff_Pos)*cosh((Ltot-pt)*ν_pos/CellData.Pos.L)))/(CC_A*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos)*ν_pos) + (CellData.Neg.L + CellData.Sep.L - pt)/(CC_A*((κ_eff_Pos+σ_eff_Pos)))  #Lee. Eqn. 4.24

#     #Flux.jl
#     D_j::Function = (Lneg,σ_eff_Neg,κ_eff_Neg,ν_neg,ν_pos,pt,CC_A) -> @. ν*(σ_eff*cosh(ν*z)+κ_eff*cosh(ν*(z-1)))/(as*F*Electrode.L*CC_A*(κ_eff+σ_eff)*sinh(ν))

#     #C_se.jl
#     D_cse::Function = z -> zeros(length(z))

#     #Phi_s.jl
#     D_ϕs::Function = (Lneg,σ_eff_Neg,κ_eff_Neg,ν_neg,ν_pos,pt,CC_A) -> @. -Electrode.L*(κ_eff*(cosh(ν)-cosh(z-1)*ν))/(CellData.Const.CC_A*σ_eff*(comb_cond_eff)*ν*sinh(ν))-Electrode.L*(σ_eff*(1-cosh(z*ν)+z*ν*sinh(ν)))/(CellData.Const.CC_A*σ_eff*(comb_cond_eff)*ν*sinh(ν)) #Transfer Function - eq. 4.19

#     #Phi_se.jl
#     D_ϕse::Function = (Lneg,σ_eff_Neg,κ_eff_Neg,ν_neg,ν_pos,pt,CC_A) -> @. Electrode.L/(CC_A*ν*sinh(ν))*((1/κ_eff)*cosh(ν*z)+(1/σ_eff)*cosh(ν*(z-1))) #Transfer Function - eq. 4.14

# end
