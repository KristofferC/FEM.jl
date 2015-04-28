using FEM
import FEM.write_data
# Nodes

# Generate geomesh and node / elementsets

# Generates a mesh of the shape known as the "Cook membrane"
# Possible mesh elements are GeoTrig for 3 node triangles,
# GeoQTrig for 6 node triangles and GeoQuad for 4 node quadraterials

nodes = [
 GeoNode2(1,[-0.664376,0.112257]),
 GeoNode2(2,[-1.07331,0.112257]),
 GeoNode2(3,[-0.800687,-0.230241]),
 GeoNode2(4,[-0.868843,0.112257]),
 GeoNode2(5,[-0.936999,-0.058992]),
 GeoNode2(6,[-0.732532,-0.058992]),
 GeoNode2(7,[-1.07331,-0.230241]),
 GeoNode2(8,[-1.07331,-0.572738]),
 GeoNode2(9,[-1.07331,-0.401489]),
 GeoNode2(10,[-0.936999,-0.401489]),
 GeoNode2(11,[-0.936999,-0.230241]),
 GeoNode2(12,[-0.255441,0.112257]),
 GeoNode2(13,[-0.255441,-0.230241]),
 GeoNode2(14,[-0.459908,0.112257]),
 GeoNode2(15,[-0.459908,-0.058992]),
 GeoNode2(16,[-0.255441,-0.058992]),
 GeoNode2(17,[-0.255441,-0.572738]),
 GeoNode2(18,[-0.664376,-0.572738]),
 GeoNode2(19,[-0.255441,-0.401489]),
 GeoNode2(20,[-0.459908,-0.401489]),
 GeoNode2(21,[-0.459908,-0.572738]),
 GeoNode2(22,[-0.528064,-0.230241]),
 GeoNode2(23,[-0.732532,-0.401489]),
 GeoNode2(24,[-0.868843,-0.572738]),
 GeoNode2(25,[-1.07331,-0.058992])]


geomesh = GeoMesh()

for node in nodes
    push!(geomesh, node)
end

elems =
[GeoQTrig(1, FEM.Vertex6(1,2,3,4,5,6)),
 GeoQTrig(2, FEM.Vertex6(7,8,3,9,10,11))  ,
 GeoQTrig(3, FEM.Vertex6(12,1,13,14,15,16)) ,
 GeoQTrig(4, FEM.Vertex6(17,13,18,19,20,21)),
 GeoQTrig(5, FEM.Vertex6(13,1,3,15,6,22)),
 GeoQTrig(6, FEM.Vertex6(18,13,3,20,22,23)),
 GeoQTrig(7, FEM.Vertex6(8,18,3,24,23,10)),
 GeoQTrig(8, FEM.Vertex6(2,7,3,25,11,5))]

for elem in elems
    push!(geomesh, elem)
end
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))


push!(geomesh, gennodeset(n->n.coords[2]>0.1122, "top", geomesh.nodes))
push!(geomesh, gennodeset(n->n.coords[2]<-0.572, "bottom", geomesh.nodes))
push!(geomesh, ElementSet("all", collect(1:length(geomesh.elements))))

# Boundary conditions
bcs = [DirichletBC(0.0, [FEM.Du, FEM.Dv], geomesh.node_sets["bottom"]),
       DirichletBC(0.0, [FEM.Dv], geomesh.node_sets["top"]),
        DirichletBC(2.806425e-5/(0.002), [FEM.Du], geomesh.node_sets["top"])]



E = 200000.e0
nu = 0.3e0
n = 1.e0
l = 1.e-2
kinf = 1e+010
lambda_0 = 4.e-2
Hg = 4.e7
Hl = 10_000.0
m = 2.0
faktor = 1
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
fp = FEM.create_feproblem_grad("grad", geomesh, [ele_section], [mat_section], bcs)

vtkexp = VTKExporter()
set_binary!(vtkexp, false)
push!(vtkexp, Stress)
push!(vtkexp, Strain)
push!(vtkexp, VonMises)
push!(vtkexp, InvFp)
push!(vtkexp, KAlpha)


solver = NRSolver(0.01, 20)

solve(solver, fp, vtkexp)
