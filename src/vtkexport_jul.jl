# Modified from Christian Denglers post at
# https://groups.google.com/forum/#!topic/julia-users/d6CEi2VFV94

#TODO, ASCII export is currently slow
function exportVTK(fp::FEProblem, filename, ascii::Bool=false)

    nr_of_elements = 0
    for section in fp.sections
        nr_of_elements += length(section.elements)
    end

    fid = open(filename, "w");

    #ASCII file header
    println(fid, "# vtk DataFile Version 3.0");
    println(fid, "FEM.jl export");
    if ascii
        println(fid, "ASCII\n");
    else
        println(fid, "BINARY\n");
    end
    println(fid, "DATASET UNSTRUCTURED_GRID");

    println(fid, "POINTS $(length(fp.nodes)) double");
    for node in fp.nodes
        c = get_coord(node)
        if ascii
            print(fid, c.x, " ");
            print(fid, c.y, " ");
            print(fid, c.z, "\n");
        else
           write(fid, bswap(c.x))
           write(fid, bswap(c.y))
           write(fid, bswap(c.z))
       end
    end


    println(fid,"CELLS ",nr_of_elements," " ,nr_of_elements * (3+1));
    for section in fp.sections
        for element in values(section.elements)
            if ascii
                print(fid, 3, " ");
                print(fid, element.vertices.v1-1, " ");
                print(fid, element.vertices.v2-1, " ");
                print(fid, element.vertices.v3-1, "\n");
            else
                write(fid, bswap(Int32(3)))
                write(fid, bswap(Int32(element.vertices.v1-1)))
                write(fid, bswap(Int32(element.vertices.v2-1)))
                write(fid, bswap(Int32(element.vertices.v3-1)))
            end
        end
    end

    println(fid,"\nCELL_TYPES ",nr_of_elements);
    for section in fp.sections
        for element in values(section.elements)
            if ascii
                print(fid, 5, "\n")
            else
                write(fid, bswap(Int32(5)))
            end
        end
    end


    print(fid, "\nPOINT_DATA $(length(fp.nodes))\n");

    println(fid, "\nVECTORS Displacement double");
    for node in fp.nodes
        u = get_displacement(node)
        if ascii
            print(fid, u.x, " ")
            print(fid, u.y, " ")
            print(fid, u.z, "\n")
        else
            write(fid, bswap(u.x))
            write(fid, bswap(u.y))
            write(fid, bswap(u.z))
        end
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
            if ascii
                print(fid, "$(strain_buf)"[2:end-1]);
                print(fid, "\n");
            else
                for strain in strain_buf
                    write(fid, bswap(strain))
                end
            end
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

            if ascii
                print(fid, "$(stress_buf)"[2:end-1]);
                print(fid, "\n");
            else
                for stress in stress_buf
                    write(fid, bswap(stress))
                end
            end
        end
    end

    close(fid);
end
