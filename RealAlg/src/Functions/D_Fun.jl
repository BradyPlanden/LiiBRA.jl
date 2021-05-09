@inline function (CellData, ν_neg, ν_pos, σ_eff_Neg, σ_eff_Pos, κ_eff_Neg, κ_eff_Sep, κ_eff_Pos)

    x5 = (0.000128*(σ_eff_Neg + κ_eff_Neg*cosh(ν_neg) - κ_eff_Neg*cosh(0*ν_neg) - σ_eff_Neg*cosh(1*ν_neg))- 0.000128*κ_eff_Neg*sinh(ν_neg)*ν_neg)/(1*κ_eff_Neg*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg)*ν_neg); 
    x6 = (-7.6e-05)/(1*κ_eff_Sep) +0.000128*((1-σ_eff_Neg/κ_eff_Neg)*tanh(ν_neg/2)-ν_neg)/(1*(κ_eff_Neg +σ_eff_Neg)*ν_neg);
    x7 = -7.6e-05/κ_eff_Sep + 0.000128*((1-σ_eff_Neg/κ_eff_Neg)* tanh(ν_neg/2)-ν_neg)/ (1*(σ_eff_Neg+κ_eff_Neg)*ν_neg) + (0.00019*(-σ_eff_Pos*cosh(ν_pos) + σ_eff_Pos*cosh(2.85316e-16*ν_pos) +κ_eff_Pos*(cosh(1*ν_pos)-1)) - 0.00019*κ_eff_Pos*sinh(ν_pos)*ν_pos)/ (1*κ_eff_Pos*(κ_eff_Pos + σ_eff_Pos)*ν_pos*sinh(ν_pos));
    x10 = -0.000128*(κ_eff_Neg*(cosh(ν_neg)-cosh((1-1)*ν_neg))+σ_eff_Neg*(1-cosh(1*ν_neg)+1*ν_neg*sinh(ν_neg)))/(1*σ_eff_Neg*(κ_eff_Neg+σ_eff_Neg)*ν_neg*sinh(ν_neg));
    x11 = -0.00019/(1*ν_pos*sinh(ν_pos))*(1/κ_eff_Pos*cosh(ν_pos*0)+1/σ_eff_Pos*cosh(ν_pos*(0-1)));
    x12 = -0.00019/(1*ν_pos*sinh(ν_pos))*(1/κ_eff_Pos*cosh(ν_pos*1)+1/σ_eff_Pos*cosh(ν_pos*(1-1)));
    x13 = -ν_pos/(1.92165e+06*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos))*(σ_eff_Pos*cosh(ν_pos*0)+ κ_eff_Pos*cosh(ν_pos*(0-1)));
    x14 = -ν_pos/(1.92165e+06*(κ_eff_Pos+σ_eff_Pos)*sinh(ν_pos))*(σ_eff_Pos*cosh(ν_pos*1)+ κ_eff_Pos*cosh(ν_pos*(1-1)));
    x17 = 0.00019*(κ_eff_Pos*(cosh(ν_pos)-cosh((1-1)*ν_pos))+σ_eff_Pos*(1-cosh(1*ν_pos)+1*ν_pos*sinh(ν_pos)))/(1*σ_eff_Pos*(κ_eff_Pos+σ_eff_Pos)*ν_pos*sinh(ν_pos));
    x18 = ν_neg/(1.39606e+06*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg))*(σ_eff_Neg*cosh(ν_neg*0)+ κ_eff_Neg*cosh(ν_neg*(0-1)));
    x19 = ν_neg/(1.39606e+06*(κ_eff_Neg+σ_eff_Neg)*sinh(ν_neg))*(σ_eff_Neg*cosh(ν_neg*1)+ κ_eff_Neg*cosh(ν_neg*(1-1)));
    x20 = 0.000128/(1*ν_neg*sinh(ν_neg))*(1/κ_eff_Neg*cosh(ν_neg*0)+1/σ_eff_Neg*cosh(ν_neg*(0-1)));
    x21 = 0.000128/(1*ν_neg*sinh(ν_neg))*(1/κ_eff_Neg*cosh(ν_neg*1)+1/σ_eff_Neg*cosh(ν_neg*(1-1)));

D = [0;0;0;0;x5;x6;x7;0;0;x10;x11;x12;x13;x14;0;0;x17;x18;x19;x20;x21];
