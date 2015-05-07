function lintrig2quadtrig(mesh::GeoMesh)
    geomesh = GeoMesh()
    n_nodes = length(mesh.nodes)

    nn_dict = Dict{(Int, Int), Int}()

    for node in mesh.nodes
        push!(geomesh, node)
    end
    combs = [[1 2]; [2 3]; [3 1]]

    n_eles = 0
    for lintrig in mesh.elements
        vs = lintrig.vertices
        vs_quad = [vs[1],vs[2],vs[3]]
        for i in 1:size(combs,1)
            c = combs[i,:]
            n1 = mesh.nodes[vs[c[1]]]
            n2 = mesh.nodes[vs[c[2]]]
            node_numbs = sort([n1.n, n2.n])

            if !haskey(nn_dict, (node_numbs[1], node_numbs[2]))
                n_nodes += 1
                nn_dict[(node_numbs[1], node_numbs[2])] = n_nodes
                c1 = get_coord(n1)
                c2 = get_coord(n2)

                node = GeoNode2(n_nodes, Point2((c1.x+c2.x)/2, (c1.y+c2.y)/2))
                push!(geomesh, node)
            end
            node_n =  nn_dict[(node_numbs[1], node_numbs[2])]
            push!(vs_quad, node_n)
        end
        n_eles += 1
        push!(geomesh, GeoQTrig(n_eles, vs_quad))
    end
    return geomesh
end