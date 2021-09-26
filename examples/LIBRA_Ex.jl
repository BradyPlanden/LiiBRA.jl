using LIBRA, Plots
plotly()

#---------- Cell Definition -----------------#
cell("LG_M50") #Alternative "Doyle_94"

#---------- DRA Loop -----------------#
@inline function DRA_loop(CellData)
    A = B = C = D = tuple()
    for Temp in 25.0:25.0:25.0
        CellData.Const.T = 273.15+Temp
        Arr_Factor = (1/CellData.Const.T_ref-1/CellData.Const.T)/R
        CellData.Const.κ = CellData.Const.κf(CellData.Const.ce0)*exp(CellData.Const.Ea_κ*Arr_Factor)
            for i in 2500:500:2500
                CellData.RA.H1 = 0:i
                CellData.RA.H2 = 0:i
                Nfft = 2^(ceil(log2(CellData.RA.Fs*CellData.RA.Tlen)))
                f = 0:Nfft-1
                s = transpose(((2im.*CellData.RA.Fs)*tan.(pi.*f./Nfft)))
                for SOC in 0.77:0.77:0.77
                    CellData.Const.SOC = SOC
                    A_DRA, B_DRA, C_DRA, D_DRA = DRA(CellData,s,f)
                    A = flatten(A,A_DRA)
                    B = flatten(B,B_DRA)
                    C = flatten(C,C_DRA)
                    D = flatten(D,D_DRA)
                end
            end
    end

return A, B, C, D
end

#---------- Generate Model -----------------#
A, B, C, D = DRA_loop(CellData)


#---------- Simulate Model -----------------#
Tk = ones(401)*298.15 #Cell Temperature
Iapp = [ones(1)*0.; ones(40)*4.8181; ones(160)*0.; ones(40)*-3.613; ones(160)*0.] #1C HPPC Experiment Current Profile
time = 0:0.25:99.75
CellV, Ce, jNeg, jPos, RtotNeg, RtotPos, η0, ηL, η_neg, η_pos, ϕ_ẽ1, ϕ_ẽ2, Uocp_Neg, Uocp_Pos, ϕ_e = Sim_Model(CellData,Iapp,Tk,A,B,C,D)


#----------- Plotting ---------------------------#
display(plot(time,CellV[1:end-1], legend=:topright,color=:blue,label="Voltage",bottom_margin=5Plots.mm, left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Voltage [V]"))
display(plot!(twinx(),time,Iapp[1:end-1],legend=:bottomright,color=:red,label="Current", left_margin = 5Plots.mm, right_margin = 15Plots.mm, ylabel = "Current [A]"))
display(plot!(title = "LG M50 - HPPC @ 80% SOC", xlabel = "Time [s]"))

display(plot(Ce))
display(plot(ϕ_e))