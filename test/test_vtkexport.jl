import FEM.write_vtk_file

facts("FEM.vtkexport") do

context("FEM.vtkexport") do

mesh = GeoMesh()

node_1 = GeoNode2([0, 0], 1)
node_2 = GeoNode2([1, 1], 2)
node_3 = GeoNode2([1, 2], 3)
node_4 = GeoNode2([0, 1], 4)
addnode!(mesh, node_1)
addnode!(mesh, node_2)
addnode!(mesh, node_3)
addnode!(mesh, node_4)

# Elements
element_1 = GeoTrig([1, 2, 3], 1)
element_2 = GeoTrig([1, 2, 4], 2)
addelem!(mesh, element_1)
addelem!(mesh, element_2)


write_vtk_file(mesh, "test.vtp", true)

end # context

end # facts