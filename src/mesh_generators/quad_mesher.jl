import FEM.LinTrigStorage

# Meshes a quadraterial area given by four corners with nx * xy * 2 triangles.

# TODO: Add quadraterial option
# TODO: Add test
# TODO: Check if the corners does not represent a convex set.
function meshquad(nx::Int, ny::Int, corners::Matrix{Float64})

    if size(corners) != (2,4)
        error("corner argument needs to be a 2x4 matrix")
    end

    if (nx <= 0 || ny <= 0)
        error("Need at least 1x1 elements")
    end

    mesh = GeoMesh()

    # The four corners
    LL = corners[:,1]
    UL = corners[:,2]
    UR = corners[:,3]
    LR = corners[:,4]

    # Node number for column i and row j
    node_nr(i::Int, j::Int) =  (nx + 1) * (i-1) + j

    # Add the nodes
    for i in 0:ny
        ratio_bounds = i / ny

        x0 = LL[1] * (1 - ratio_bounds) + ratio_bounds * UL[1]
        x1 = LR[1] * (1 - ratio_bounds) + ratio_bounds * UR[1]

        y0 = LL[2] * (1 - ratio_bounds) + ratio_bounds * UL[2]
        y1 = LR[2] * (1 - ratio_bounds) + ratio_bounds * UR[2]

        for j in 0:nx
            ratio = j / nx
            x = x0 * (1 - ratio) + ratio * x1
            y = y0 * (1 - ratio) + ratio * y1
            push!(mesh, GeoNode2(node_nr(i+1, j+1), [x,y]))
        end
    end

    # Add the elements
    n_elem = 0
    for i in 1:ny
        for j in 1:nx
            n_elem += 1
            elem1 = GeoTrig(n_elem, Vertex3(node_nr(i, j), node_nr(i+1, j), node_nr(i+1, j+1)))
            push!(mesh, elem1)
            n_elem += 1
            elem2 = GeoTrig(n_elem, Vertex3(node_nr(i, j), node_nr(i+1, j+1), node_nr(i, j+1)))
            push!(mesh, elem2)
        end
    end
    return mesh
end

# Generates what is known as the "Cook membrane"
function gencook(nx, ny, scale = 1.0)
    m = scale * [[0.0; 0.0] [0.0; 44.0] [48.0; 60.0] [48.0; 44.0]]
    meshquad(nx, ny, m)
end
