using FEM

E = 200000.e0
nu = 0.3e0
l = 1.e-2
Hg = 4.e7
Hl = 10_000
m = 2.0
sy = 0.0
tstar = 1000.0
angles = [45.0, 105.0]
nslip = 2

mat = FEM.gradmekhprimalmod().GradMekh(E, nu, l,
               Hg, Hl, m, sy, tstar,
               angles, nslip)

matstat = FEM.GradMekhPrimalMod.GradMekhPrimalMS(zeros(2), zeros(2),zeros(2), zeros(6), zeros(6))
temp_matstat = copy(matstat)

#τ = [46.15384615384616,-66.37819326614503]
τ = [0.0, 0.0]
γ = [-0.0,-0.0]

matstat = FEM.GradMekhPrimalMod.GradMekhPrimalMS([0.0,0.0],[0.0,0.0],[0.1,0.1],[0.0,0.0,0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0,0.0,0.0])
temp_matstat = FEM.GradMekhPrimalMod.GradMekhPrimalMS([-0.0,-0.0],[63062.62890470481,63062.82377457908],[-0.006410256410256407,0.003205128205128202],[0.0,0.0,0.0,0.0,0.0,0.0],[-0.005128205128205127,0.005128205128205127,0.0,0.005128205128205128,0.0,0.0])



FEM.stress(mat, matstat, temp_matstat,τ, γ)
NSLIP = 2

∆t = 2.000000000000000e-003 #TODO: Fix
α = 1
H = 1.0
∆γ = γ[α] - matstat.n_γ[α]
tstar = mat.tstar
m = mat.m
κ = zeros(NSLIP)

R, J = FEM.GradMekhPrimalMod.resid(mat.σy, ∆γ, ∆t, τ[α], κ[α], mat.m, tstar)
#κ[α] += H
#R_f = FEM.GradMekhPrimalMod.resid(mat.σy, ∆γ, ∆t, τ[α], κ[α], mat.m, tstar)
#κ[α] -= H
#J = (R_f - R) / H
κ[α] = κ[α] - R / J
