using FEM
function gen_quadmesh(n_x::Int, n_y::Int, corners::Matrix)

    # Corners in lower left, upper left, upper right, lower right order

    mesh = Mesh()
    mesh_coords_x = zeros(n_y+1, n_x+1)
    mesh_coords_y = zeros(n_y+1, n_x+1)

    LL = corners[1,:]
    UL = corners[2,:]
    UR = corners[3,:]
    LR = corners[4,:]

    println(LL)
    println(UL)
    println(UR)
    println(LR)

    for i in 0:n_y
        ratio_y = i / n_y

        y0 = LL[2]*(1-ratio_y) + ratio_y * UL[2]
        y1 = LR[2]*(1-ratio_y) + ratio_y * UR[2]

        println(y1)

        for j in 0:n_x
            ratio_x = j / n_x
            x0 = LL[1]*(1-ratio_x) + ratio_x * LR[1]
            x1 = UL[1]*(1-ratio_x) + ratio_x * UR[1]
            x = x1 * ratio_x
            y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            mesh_coords_x[i+1, j+1] = x
            mesh_coords_y[i+1, j+1] = y
        end
    end

    node_nr(i, j) =  (n_x + 1) * (i-1) + j

    # Add nodes to mesh
    for i in 1:n_y+1
        for j in 1:n_x+1
            node_n = n_x * (i-1) + j
            x = mesh_coords_x[i, j]
            y = mesh_coords_y[i, j]
            addnode!(mesh, Node([x,y], node_nr(i, j)))
        end
    end

    # Add elements to mesh
    n_elem = 0
    for i in 1:n_y
        for j in 1:n_y
            n_elem += 1
            elem1 = LinTrig([node_nr(i, j), node_nr(i+1, j), node_nr(i+1, j+1)], n_elem)
            addelem!(mesh, elem1)
            n_elem += 1
            elem2 = LinTrig([node_nr(i, j), node_nr(i+1, j+1), node_nr(i, j+1)], n_elem)
            addelem!(mesh, elem2)
        end
    end
    return mesh
end

function gen_cook(nx, ny)
    m = [[0.0 0.0];
         [0.0 0.044];
         [0.048 0.06];
         [0.048 0.044]]
    gen_quadmesh(nx, ny,m)
end
