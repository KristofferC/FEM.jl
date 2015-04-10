# Meshes a quadraterial area given by four corners with nx * xy * 2 triangles.

# TODO: Add quadraterial option
# TODO: Add test
# TODO: Check if the corners does not represent a convex set.

   # Node number for column i and row j


function meshquad(nx::Int, ny::Int, corners::Matrix{Float64},
                 ele_type::Union(Type{GeoQuad}, Type{GeoTrig}, Type{GeoQTrig}))

    if size(corners) != (2,4)
        error("corner argument needs to be a 2x4 matrix")
    end

    if (nx <= 0 || ny <= 0)
        error("Need at least 1x1 elements")
    end

    mesh = GeoMesh()

    #... you what Lint?
    @lintpragma("Ignore use of undeclared variable :")
    # The four corners
    LL = corners[:, 1]
    UL = corners[:, 2]
    UR = corners[:, 3]
    LR = corners[:, 4]




    if ele_type == GeoQTrig
        np_y = 2*ny
        np_x = 2*nx
    else
        np_y = ny
        np_x = ny
    end

     @inline node_nr(i::Int, j::Int) = (np_x + 1) * (i-1) + j


    # Add the nodes
    for i in 0:np_y
        ratio_bounds = i / np_y

        x0 = LL[1] * (1 - ratio_bounds) + ratio_bounds * UL[1]
        x1 = LR[1] * (1 - ratio_bounds) + ratio_bounds * UR[1]

        y0 = LL[2] * (1 - ratio_bounds) + ratio_bounds * UL[2]
        y1 = LR[2] * (1 - ratio_bounds) + ratio_bounds * UR[2]

        for j in 0:np_x
            ratio = j / np_x
            x = x0 * (1 - ratio) + ratio * x1
            y = y0 * (1 - ratio) + ratio * y1
            push!(mesh, GeoNode2(node_nr(i+1, j+1), [x,y]))
        end
    end

    # Add the elements
    n_elem = 0
       if ele_type == GeoQTrig
        jump = 2
    else
        jump = 1
    end

    for i in 1:jump:jump*ny
        for j in 1:jump:jump*nx
            if ele_type == GeoTrig
                n_elem += 1
                elem1 = GeoTrig(n_elem, Vertex3(node_nr(i, j), node_nr(i+1, j), node_nr(i+1, j+1)))
                push!(mesh, elem1)
                n_elem += 1
                elem2 = GeoTrig(n_elem, Vertex3(node_nr(i, j), node_nr(i+1, j+1), node_nr(i, j+1)))
                push!(mesh, elem2)
            elseif ele_type == GeoQTrig
                n_elem += 1
                elem1 = GeoQTrig(n_elem, Vertex6(node_nr(i, j), node_nr(i+2, j), node_nr(i+2, j+2),
                                                 node_nr(i+1, j), node_nr(i+2, j+1), node_nr(i+1, j+1)))
                push!(mesh, elem1)
                n_elem += 1
                elem1 = GeoQTrig(n_elem, Vertex6(node_nr(i, j), node_nr(i, j+2), node_nr(i+2, j+2),
                                                 node_nr(i, j+1), node_nr(i+1, j+2), node_nr(i+1, j+1)))
                push!(mesh, elem1)
            else
                n_elem += 1
                elem = GeoQuad(n_elem, Vertex4(node_nr(i,j), node_nr(i,j+1), node_nr(i+1,j+1), node_nr(i+1,j)))
                push!(mesh, elem)
            end
        end
    end
    return mesh
end

# Generates what is known as the "Cook membrane"
function gencook(nx, ny, ele_type::Union(Type{GeoQuad}, Type{GeoTrig}, Type{GeoQTrig}), scale = 1.0)
    m = scale * [[0.0; 0.0] [0.0; 44.0] [48.0; 60.0] [48.0; 44.0]]
    meshquad(nx, ny, m, ele_type)
end
