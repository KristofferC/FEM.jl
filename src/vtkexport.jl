function create_vtk_object(mesh::Mesh)

    element_vtk_table = Dict{Element, Any}
    element_vtk_table[LinTrig]

    n_nodes = length(mesh.nodes)

    # Add points to vtkPoints
    points = pycall(vtk.vtkPoints, PyAny)
    for node in mesh.nodes
        d_coords = node.coords
        push!(d_coords, 0)
        points[:InsertNextPoint](d_coords[1], d_coords[2], d_coords[3])
    end

    # Add elements to vtkCellArray
    elems_vtk = pycall(vtk.vtkCellArray, PyAny)
    for element in mesh.elements
        elem_vtk = pycall(vtk.vtkTriangle, PyAny)
        for (i, vertex) in enumerate(element.vertices)
            ids = elem_vtk[:GetPointIds]()
            ids[:SetId](i-1, vertex - 1)
        end
        elems_vtk[:InsertNextCell](elem_vtk)
    end

   # Add displacements
    disp_array = pycall(vtk.vtkDoubleArray, PyAny)
    disp_array[:SetNumberOfComponents](3)
    disp_array[:SetName]("Displacement")
    disp_array[:SetNumberOfTuples](n_nodes)
    for node in mesh.nodes:
        # TODO: Remove hardcoding, should make like a
        # get displacements in node class
        disp = get_displacements(node)
        disp_array[:SetTuple3](node.n, disp[1], disp[2], disp[3])
    end



    # Create a polydata to store everything in
    polydata = pycall(vtk.vtkPolyData, PyAny)

    point_data = polydata[:GetPointData]()
    # Add the points to the polydata object
    point_data[:SetPoints](points)

    # Add the elements to the polydata object
    point_data[:SetPolys](elems_vtk)

    # Add displacements
    point_data[:AddArray](disp_array)

    return polydata
end


function write_vtk_file(mesh::FEMesh, name::ASCIIString, write_ascii::Bool = false)
    polydata = create_vtk_object(mesh)

    writer = pycall(vtk.vtkXMLPolyDataWriter, PyAny)
    writer[:SetFileName](name)
    if vtk.VTK_MAJOR_VERSION <= 5
        writer[:SetInput](polydata)
    else
        writer[:SetInputData](polydata)
    end

    if write_ascii
        writer[:SetDataModeToAscii]()
    end
    writer[:Write]()
end
