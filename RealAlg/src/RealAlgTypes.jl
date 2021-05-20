using Parameters

@with_kw mutable struct TransferFun
tfs =   [[C_e, Phi_e, C_se, Phi_s, Phi_se, j, C_se, Phi_s, j, Phi_se] ["Na", "Na", "Pos", "Pos", "Pos", "Pos", "Neg", "Neg", "Neg", "Neg"] [Number[0, 4.26E-05, 8.52E-05, 9.72E-05, 1.35E-04, 1.73E-04], Number[4.26E-05, 8.52E-05, 9.72E-05, 1.35E-04, 1.73E-04], Number[0,1], Number[1],Number[0,1],Number[0,1],Number[0,1],Number[1],Number[0,1],Number[0,1]]]
end