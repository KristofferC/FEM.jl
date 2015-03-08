import FEM.stiffness

facts("FEM.Element") do


context("FEM.Element.LinTrig") do


    nodes = [Node([0, 0], 1), Node([1, 1], 2), Node([1, 2], 3)]

    elem = LinTrig([1, 2, 3], 1)

    mat = LinearIsotropic(200e9, 0.3)

    stiffness(elem, nodes, mat)


end # context

end # facts
