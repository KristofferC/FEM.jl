import FEM.LinTrigStorage

# Meshes a quadraterial area given by four corners with nx * xy * 2 triangles.

# TODO: Add quadraterial option
# TODO: Add test
# TODO: Check if the corners does not represent a convex set.
function meshquad(nx::Int, ny::Int, corners::Matrix{Float64})

    if size(corners) != (2,4)
        error("corner argument needs to be a 2x4 matrix")
    end

    if (nx <= 1 || ny <= 1)
        error("Need at least 2x2 elements")
    end

    mesh = Mesh()

    # The four corners
    LL = corners[:,1]
    UL = corners[:,2]
    UR = corners[:,3]
    LR = corners[:,4]

    # Node number for column i and row j
    node_nr(i, j) =  (nx + 1) * (i-1) + j

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
            addnode!(mesh, Node2([x,y], node_nr(i+1, j+1)))
        end
    end

    interp = LinTrigInterp()
    lts = LinTrigStorage()
    gps = [GaussPoint2(Point2(1/3, 1/3), 0.5)]
    # Add the elements
    n_elem = 0
    for i in 1:ny
        for j in 1:nx
            n_elem += 1
            elem1 = LinTrig([node_nr(i, j), node_nr(i+1, j), node_nr(i+1, j+1)], n_elem, interp, lts, gps)
            addelem!(mesh, elem1)
            n_elem += 1
            elem2 = LinTrig([node_nr(i, j), node_nr(i+1, j+1), node_nr(i, j+1)], n_elem, interp, lts, gps)
            addelem!(mesh, elem2)
        end
    end
    return mesh
end

# Generates what is known as the "Cook membrane"
function gencook(nx, ny)
    m = [[0.0; 0.0] [0.0; 0.044] [0.048; 0.06] [0.048; 0.044]]
    meshquad(nx, ny, m)
end
