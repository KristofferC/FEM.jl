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
    println(fid, "DATASET UNSTRUCTURED_GRID");

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

    print(fid, "\nCELL_DATA $(nr_of_elements)\n");

    println(fid, "\nTENSORS Strain double");
    strain_buf = zeros(3,3)
     for section in fp.sections
        mat = section.material
        for element in values(section.elements)
            mat_stats = mat.matstats[element.n]
            strain = mat.matstats[element.n][1].strain
            strain_buf[1,1] = strain[1]
            strain_buf[2,2] = strain[2]
            strain_buf[3,3] = strain[3]
            strain_buf[1,3] = strain[3,1] = strain[4]

            print(fid, "$(strain_buf)"[2:end-1]);
            print(fid, "\n");
        end
    end


    println(fid, "\nTENSORS Stress double");
    stress_buf = zeros(3,3)
     for section in fp.sections
        mat = section.material
        for element in values(section.elements)
            mat_stats = mat.matstats[element.n]
            stress = mat.matstats[element.n][1].stress
            stress_buf[1,1] = stress[1]
            stress_buf[2,2] = stress[2]
            stress_buf[3,3] = stress[3]
            stress_buf[1,3] = stress[3,1] = stress[4]

            print(fid, "$(stress_buf)"[2:end-1]);
            print(fid, "\n");
        end
    end


    #print(fid, "\nCELL_DATA $(nr_of_elements)\n");

    #print(fid,"TENSORS ", strain," double\n");
    #for section in fp.sections
    #    matstats = section.material.matstats
    #    end
    #end

    close(fid);
end
