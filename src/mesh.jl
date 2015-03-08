immutable Mesh
    nodes::Vector{Node}
    elements::Vector{Element}
    element_sets::Dict{ASCIIString, ElementSet}
    node_sets::Dict{ASCIIString, NodeSet}
end


function Mesh()
    nodes = Array(Node, 0)
    elements = Array(Element, 0)
    element_sets = Dict{ASCIIString, ElementSet}()
    node_sets =  Dict{ASCIIString, NodeSet}()

    Mesh(nodes, elements, element_sets, node_sets)
end

addnode!(mesh::Mesh, node::Node) = push!(mesh.nodes, node)
function addnodes!(mesh::Mesh, nodes::Vector{Node})
    for n in nodes
        addnode!(mesh, n)
    end
end

addelem!(mesh::Mesh, elem::Element) = push!(mesh.elements, elem)

function addelemset!(mesh::Mesh, elem_set::ElementSet)
    mesh.element_sets[elem_set.name] = elem_set
end

function addnodeset!(mesh::Mesh, node_set::NodeSet)
     mesh.node_sets[node_set.name] = node_set
end



