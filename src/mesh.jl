immutable Mesh
    nodes::Vector{Node2}
    elements::Vector{LinTrig}
    element_sets::Dict{ASCIIString, ElementSet}
    node_sets::Dict{ASCIIString, NodeSet}
end


function Mesh()
    nodes = Array(Node2, 0)
    elements = Array(LinTrig, 0)
    element_sets = Dict{ASCIIString, ElementSet}()
    node_sets =  Dict{ASCIIString, NodeSet}()

    Mesh(nodes, elements, element_sets, node_sets)
end

addnode!(mesh::Mesh, node::AbstractNode) = push!(mesh.nodes, node)
function addnodes!(mesh::Mesh, nodes::Vector{AbstractNode})
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



immutable Section
    material::Material
    elements::Set{Int} # The elements in the section
end

function Section(mat::Material)
    elements = Set(Int[])
    Section(mat, elements)
end

function addelemset!(section::Section, elemset::ElementSet)
    for i in elemset.element_ids
        push!(section.elements, i)
    end
end
