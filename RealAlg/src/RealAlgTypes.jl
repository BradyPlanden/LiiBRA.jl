using Parameters

@with_kw mutable struct TransferFun
tfs =   [[C_e, Phi_e, C_se, Phi_s, Phi_se, j, C_se, Phi_s, j, Phi_se] ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] [Number[0, 128e-6, 204e-6, 394e-6], Number[128e-6, 204e-6, 394e-6], Number[0,1], Number[1],Number[0,1],Number[0,1],Number[0,1],Number[1],Number[0,1],Number[0,1]]]
end