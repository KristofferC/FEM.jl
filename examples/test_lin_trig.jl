using FEM

# Nodes
node_1 = Node(1, [0, 0])
node_2 = Node(2, [1, 1])
node_3 = Node(3, [1, 2])
node_4 = Node(4, [0, 1])


# Elements
element_1 = TrigPlaneStrain(1, [1, 2, 3])
element_2 = TrigPlaneStrain(2, [1, 2, 4])

# Sets
bottom_set = NodeSet("y0", [1])
top_set = NodeSet("x0", [2, 3])
node_sets = [bottom_set, top_set]

# Element set
element_set = ElementSet("all", [1, 2])

# Create mesh
mesh = Mesh()
mesh.add_node(node_1)
mesh.add_node(node_2)
mesh.add_node(node_3)
mesh.add_node(node_4)
mesh.add_element(element_1)
mesh.add_element(element_2)
mesh.add_element_set(element_set)
mesh.add_node_set(bottom_set)
mesh.add_node_set(top_set)
