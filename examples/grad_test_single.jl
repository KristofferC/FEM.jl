using FEM
import FEM.write_data
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
kinf = 1e+010
lambda_0 = 4.e-2
Hg = 4.e7
Hl = 10_000.0
m = 2.0
faktor = 0.5
sy = 1000.0
tstar = 1000.0
c_dam = 1.0
angles = [45.0, 105.0]
nslip = 2

mat = GradMekh(E, nu, n, l, kinf, lambda_0,
               Hg, Hl, m, faktor, sy, tstar,
               c_dam, angles, nslip)

mat_section = MaterialSection(mat)
push!(mat_section, geomesh.element_sets["all"])

ele_section = ElementSection(GradTrig)
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


fp.nodes[elem.vertices[1]].dofs[1].value = ulem[1]
fp.nodes[elem.vertices[1]].dofs[2].value = ulem[2]
fp.nodes[elem.vertices[1]].dofs[3].value = ulem[13]
fp.nodes[elem.vertices[1]].dofs[4].value = ulem[14]
fp.nodes[elem.vertices[1]].dofs[5].value = ulem[15]
fp.nodes[elem.vertices[1]].dofs[6].value = ulem[16]

fp.nodes[elem.vertices[2]].dofs[1].value = ulem[3]
fp.nodes[elem.vertices[2]].dofs[2].value = ulem[4]
fp.nodes[elem.vertices[2]].dofs[3].value = ulem[17]
fp.nodes[elem.vertices[2]].dofs[4].value = ulem[18]
fp.nodes[elem.vertices[2]].dofs[5].value = ulem[19]
fp.nodes[elem.vertices[2]].dofs[6].value = ulem[20]

fp.nodes[elem.vertices[3]].dofs[1].value = ulem[5]
fp.nodes[elem.vertices[3]].dofs[2].value = ulem[6]
fp.nodes[elem.vertices[3]].dofs[3].value = ulem[21]
fp.nodes[elem.vertices[3]].dofs[4].value = ulem[22]
fp.nodes[elem.vertices[3]].dofs[5].value = ulem[23]
fp.nodes[elem.vertices[3]].dofs[6].value = ulem[24]

fp.nodes[elem.vertices[4]].dofs[1].value = ulem[7]
fp.nodes[elem.vertices[4]].dofs[2].value = ulem[8]

fp.nodes[elem.vertices[5]].dofs[1].value = ulem[9]
fp.nodes[elem.vertices[5]].dofs[2].value = ulem[10]

fp.nodes[elem.vertices[6]].dofs[1].value = ulem[11]
fp.nodes[elem.vertices[6]].dofs[2].value = ulem[12]


state_gp_1 = [
  1.000000000000000
  1.000000000000000
  1.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000]


state_gp_2 = [
  1.000000000000000
  1.000000000000000
  1.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000
  0.000000000000000 ]

 state_gp_3 =
 [ 1+8.935443305939650e-004
  1+-8.935443305939650e-004
  1+0.000000000000000e+000
  2.394244780597564e-004
  0.000000000000000e+000
  0.000000000000000e+000
  0.000000000000000e+000
 -3.334752868526071e-003
  0.000000000000000e+000
  0.000000000000000e+000
  -3.564628423099781e-003
  0.000000000000000e+000
  2.147755004754846e-005
  7.910158601032613e-003 ]

FEM.fill_from_start!(elem.matstats[1].state, state_gp_1)
FEM.fill_from_start!(elem.matstats[2].state , state_gp_2)
FEM.fill_from_start!(elem.matstats[3].state , state_gp_3)

fe = FEM.intf(elem, mat, fp.nodes)
println(fe)

Ke = FEM.assembleK(fp)
#println(full(Ke))


#println(fe)