using FEM
import FEM.write_data
import FEM.mod_value
import FEM.AddFun
import FEM.create_sparse_structure
import FEM.get_colptrs
import FEM.assembleK!
# Nodes

# t > 1.8
# breakpoint in magnus code
nodes = [
 GeoNode2(1,[4.99999388655151, 10.0000000000000]),
 GeoNode2(2,[0.00000000000000, 10.0000000000000]),
 GeoNode2(3,[3.33333333333333, 4.99999270067665]),
 GeoNode2(4,[2.49999694327576, 10.0000000000000]),
 GeoNode2(5,[1.66666055321818, 7.49999635033832]),
 GeoNode2(6,[4.16665749649394, 7.49999635033832])]


geomesh = GeoMesh()

for node in nodes
    push!(geomesh, node)
end

elems =
[GeoQTrig(1, FEM.Vertex6(1,2,3,4,5,6))]

for elem in elems
    push!(geomesh, elem)
end
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))



E = 200000.e0
nu = 0.3e0
n = 1.e0
l = 1.e-2
Hg = 4.e7
Hl = 10_000.0
m = 2.0
sy = 1000.0
tstar = 1000.0
angles = [45.0, 105.0]
nslip = 2

mat = FEM.gradmekhmodjlsmall().GradMekh(E, nu, n, l,
               Hg, Hl, m, sy, tstar,
               angles, nslip)

mat_section = MaterialSection(mat)
push!(mat_section, geomesh.element_sets["all"])

ele_section = ElementSection(FEM.gradtrigmod().GradTrig)
push!(ele_section, geomesh.element_sets["all"])


# Create the fe problem
fp = FEM.create_feproblem_grad("grad", geomesh, [ele_section], [mat_section])

elem = fp.sections[1].elements[1]
mat = fp.sections[1].material
ulem =
 [
 0.224750000000000
 0.000000000000000e+000
 0.224750000000000
 0.000000000000000e+000
 0.116070360910861
 -6.535293225632476e-005
 0.224750000000000
 0.000000000000000e+000
 0.180637013102801
 1.165487436149526e-002
 0.178834621542997
 1.670855365259811e-003
 -3.705262145828408e-004
 -7.275567962847590e-023
 4.859544264794582e-004
 -1.371021480168575e-019
 -2.323140690434515e-020
 -5.970224756865314e-020
 7.384073653782159e-021
 -2.108672096650753e-021
 1.058646882120961e-003
 -5.238405128859649e-006
 -2.258113026569606e-003
 4.373293207049114e-004 ]


mod_value(fp.dof_vals, fp.nodes[elem.vertices[1]].dofs[1], AddFun(), ulem[1])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[1]].dofs[2], AddFun(), ulem[2])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[1]].dofs[3], AddFun(), ulem[13])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[1]].dofs[4], AddFun(), ulem[14])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[1]].dofs[5], AddFun(), ulem[15])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[1]].dofs[6], AddFun(), ulem[16])

mod_value(fp.dof_vals, fp.nodes[elem.vertices[2]].dofs[1], AddFun(), ulem[3])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[2]].dofs[2], AddFun(), ulem[4])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[2]].dofs[3], AddFun(), ulem[17])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[2]].dofs[4], AddFun(), ulem[18])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[2]].dofs[5], AddFun(), ulem[19])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[2]].dofs[6], AddFun(), ulem[20])

mod_value(fp.dof_vals, fp.nodes[elem.vertices[3]].dofs[1], AddFun(), ulem[5])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[3]].dofs[2], AddFun(), ulem[6])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[3]].dofs[3], AddFun(), ulem[21])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[3]].dofs[4], AddFun(), ulem[22])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[3]].dofs[5], AddFun(), ulem[23])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[3]].dofs[6], AddFun(), ulem[24])

mod_value(fp.dof_vals, fp.nodes[elem.vertices[4]].dofs[1], AddFun(), ulem[7])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[4]].dofs[2], AddFun(), ulem[8])

mod_value(fp.dof_vals, fp.nodes[elem.vertices[5]].dofs[1], AddFun(), ulem[9])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[5]].dofs[2], AddFun(), ulem[10])

mod_value(fp.dof_vals, fp.nodes[elem.vertices[6]].dofs[1], AddFun(), ulem[11])
mod_value(fp.dof_vals, fp.nodes[elem.vertices[6]].dofs[2], AddFun(), ulem[12])


I = Float64[1,1,1,0,0,0,0,0,0]

#=
inv_nFpstate_gp_3 =
[ 8.935443305939650e-004
-8.935443305939650e-004
0.000000000000000e+000
2.394244780597564e-004
0.000000000000000e+000
0.000000000000000e+000
0.000000000000000e+000
2.394244780597564e-004
0.000000000000000e+000]
=#

inv_nFpstate_gp_3 =
[ 8.935443305939650e-004
-8.935443305939650e-004
0.000000000000000e+000
2.394244780597564e-004
0.000000000000000e+000
0.000000000000000e+000
0.000000000000000e+000
-3.334752868526071e-003
0.000000000000000e+000]

n_k_alpha_gp_3 = [
0.000000000000000e+000
-3.564628423099781e-003 ]

n_lambda_alpha_gp_3 = [
0.000000000000000e+000
2.147755004754846e-0059 ]

FEM.fill_from_start!(elem.matstats[1].n_ε_p, zeros(9))
FEM.fill_from_start!(elem.matstats[2].n_ε_p, zeros(9))
FEM.fill_from_start!(elem.matstats[3].n_ε_p, inv_nFpstate_gp_3)
FEM.fill_from_start!(elem.matstats[3].n_k, n_k_alpha_gp_3)
FEM.fill_from_start!(elem.matstats[3].n_∆λ, n_lambda_alpha_gp_3)


fe = FEM.intf(elem, mat, fp.nodes, fp.dof_vals)
#println(fe)

    K = create_sparse_structure(fp::FEProblem)
    colptrs = get_colptrs(K, fp::FEProblem)

assembleK!(K, fp, colptrs, fp.dof_vals)
#println(full(K))


#println(fe)
