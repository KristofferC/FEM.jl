facts("FEM.Mesh") do

mesh = GeoMesh()

context("FEM.Mesh.Nodes") do

node_1 = GeoNode2(1, [0, 0])
node_2 = GeoNode2(2, [1, 1])
node_3 = GeoNode2(3, [1, 2])
node_4 = GeoNode2(4, [0, 1])
push!(mesh, node_1)
push!(mesh, node_2)
push!(mesh, node_3)
push!(mesh, node_4)

@fact is(node_2, mesh.nodes[2]) => true
@fact length(mesh.nodes) => 4

end # context

context("FEM.GeoMesh.Elements") do

# Elements
element_1 = GeoTrig(1, [1, 2, 3])
element_2 = GeoTrig(2, [1, 2, 4])
push!(mesh, element_1)
push!(mesh, element_2)

@fact length(mesh.elements) => 2
@fact is(element_2, mesh.elements[2]) => true

end # context

context("FEM.GeoMesh.NodeSet") do

# Sets
bottom_set = NodeSet("y0", [1])
top_set = NodeSet("x0", [2, 3])

push!(mesh, bottom_set)
push!(mesh, top_set)

@fact length(mesh.node_sets) => 2
@fact is(bottom_set, mesh.node_sets["y0"]) => true

node_s = gennodeset(x->x[1] > 0.5, "test", mesh.nodes)
@fact node_s.name =>"test"
@fact node_s.node_ids => [2, 3]

end # context

context("FEM.GeoMesh.ElementSet") do

# Element set
element_set = ElementSet("all", [1, 2])
push!(mesh, element_set)

@fact length(mesh.node_sets) => 2
@fact is(element_set, mesh.element_sets["all"]) => true

end # context

end # facts


