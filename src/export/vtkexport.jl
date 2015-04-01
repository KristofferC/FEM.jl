function create_vtk_object(fp::FEProblem)

    n_nodes = length(fp.nodes)

    # Add points to vtkPoints
    points = pycall(vtk.vtkPoints, PyAny)
    for node in fp.nodes
        d_coords = node.coords
        points[:InsertNextPoint](d_coords[1], d_coords[2], 0)
    end

    # Add elements to vtkCellArray
    elems_vtk = pycall(vtk.vtkCellArray, PyAny)
    for section in fp.sections
        for element in section.elements
            elem_vtk = pycall(vtk.vtkTriangle, PyAny)
            for (i, vertex) in enumerate(element.vertices)
                ids = elem_vtk[:GetPointIds]()
                ids[:SetId](i-1, vertex - 1)
            end
            elems_vtk[:InsertNextCell](elem_vtk)
        end
    end

   # Add displacements
    disp_array = pycall(vtk.vtkDoubleArray, PyAny)
    disp_array[:SetNumberOfComponents](3)
    disp_array[:SetName]("Displacement")
    disp_array[:SetNumberOfTuples](n_nodes)
    for node in fp.nodes
        disp = get_displacement(node)
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


function write_vtk_file(fp::FEProblem, name::ASCIIString, write_ascii::Bool = false)
    polydata = create_vtk_object(fp)

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
