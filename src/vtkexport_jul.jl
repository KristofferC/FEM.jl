# Modified from Christian Denglers post at
# https://groups.google.com/forum/#!topic/julia-users/d6CEi2VFV94

##### exportVTK
function exportVTK(fp::FEProblem, filename)

    nr_of_elements = 0
    for section in fp.sections
        nr_of_elements += length(section.elements)
    end

    fid = open(filename, "w");

    #ASCII file header
    println(fid, "# vtk DataFile Version 3.0");
    println(fid, "FEM.jl export");
    println(fid, "ASCII\n");
    println(fid, "DATASET POLYDATA");

    println(fid, "POINTS $(length(fp.nodes)) double");
    for node in fp.nodes
        print(fid, node.coords.x, " ");
        print(fid, node.coords.y, " ");
        print(fid, 0.0, "\n");
    end

    print(fid,"\nCELLS ",nr_of_elements," " ,nr_of_elements * (3+1),"\n");

    for section in fp.sections
        for element in values(section.elements)
            print(fid, 3, " ");
            print(fid, element.vertices.v1-1, " ");
            print(fid, element.vertices.v2-1, " ");
            print(fid, element.vertices.v3-1, "\n");
        end
    end

    print(fid,"\nCELL_TYPES ",nr_of_elements,"\n");
    for section in fp.sections
        for element in values(section.elements)
            print(fid, 5, "\n")
        end
    end

    print(fid, "\nPOINT_DATA $(length(fp.nodes))\n");

    println(fid, "\nVECTORS Displacement double");
    for node in fp.nodes
        u = get_displacement(node)
        print(fid, u.x, " ")
        print(fid, u.y, " ")
        print(fid, u.z, "\n")
    end

    print(fid,"CELL_DATA ",length(scalars),"\n");
    print(fid,"TENSORS ", strain," double\n");
    for section in fp.sections
        for element in values(section.elements)
            element.ms.get_strain()

        end
    end

    close(fid);
end
