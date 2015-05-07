# Functionality for reading a mesh from COMSOL
# Currently 2D, only 1 object, only domain elements, currently only triangles

# 21 lines before the coordinates start
# this is not robust
const LINES_BEFORE_CORDS = 22
const LINES_BETWEEN_CORDS_ELEMENT = 10

function read_mphtxt(filename)
    mesh = GeoMesh()
    open(filename, "r") do f
        for i in 1:LINES_BEFORE_CORDS
            readline(f)
        end


        read_coordinates(f, mesh)

        for i in 1:LINES_BETWEEN_CORDS_ELEMENT
            readline(f)
        end

        read_elements(f, mesh)
    end
    return mesh
end

function read_coordinates(f, mesh)
    n_nodes = 1
    while true
        line = strip(readline(f))
        # The list of coordinates is broken by an empty line
        line == "" && break
        coords = split(line, " ")
        coords = [parse(Float64, c) for c in coords]
        push!(mesh, GeoNode2(n_nodes, Point2(coords[1], coords[2])))
        n_nodes += 1
    end
end

function read_elements(f, mesh)
    n_eles = 1
    while true
        line = strip(readline(f))
        # The list of elements is broken by an empty line
        line == "" && break
        vs = split(line, " ")
        vs = [parse(Int, v)+1 for v in vs]
        push!(mesh, GeoTrig(n_eles, Vertex3(vs[1], vs[2], vs[3])))
        n_eles += 1
    end
end