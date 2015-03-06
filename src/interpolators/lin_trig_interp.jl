immutable LinTrigInterp <: Interpolator
end


# Shape functions in local coordinates
function setN!(interp::LinTrigInterp, loc_coords::Vector{Float64}, N::Vector{Float64})
    ksi = local_coords[0]
    eta = local_coords[1]

    N[0] = ksi
    N[1] = eta
    N[2] = 1.0 - ksi - eta
end

function area(interp::LinTrigInterp,, vertices, mesh):
        [x1, x2, x3, y1, y2, y3] = self._give_cords(vertices, mesh)
        return 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
end

function eval_dNdx(self, local_coords, vertices, mesh):

        dN = setN(local_coords)
        J = self.give_J(local_coords, vertices, mesh)
        J_inv = inv_2x2(J)
        dNdx = dN.dot(J_inv.transpose())

        return dNdx
end


function derivatives!(interp::LinTrigInterp, loc_coords::Vector{Float64}, dN::Vector{Float64}):

        # Derivative w.r.t ksi
        dN[0, 0] = 1.0
        dN[1, 0] = 0.0
        dN[2, 0] = -1.0

        # Derivative w.r.t eta
        dN[0, 1] = 0.0
        dN[1, 1] = 1.0
        dN[2, 1] = -1.0
end

    def give_det_J(self, local_coords, vertices, mesh):
                J = self.give_J(local_coords, vertices, mesh)
       ireturn determinant_2x2(J)

def give_J(self, local_coords, vertices, mesh):

    J = np.zeros((2, 2), dtype=np.float64)
    dN = self. give_derivatives(local_coords)

    for row in xrange(0, 3):
        x = mesh.nodes[vertices[row]].coordinates[0]
        y = mesh.nodes[vertices[row]].coordinates[1]

        J[0, 0] += dN[row, 0] * x
        J[0, 1] += dN[row, 0] * y
        J[1, 0] += dN[row, 1] * x
        J[1, 1] += dN[row, 1] * y

    return J





    def _give_cords(self, vertices, mesh):
        return([mesh.nodes[vertices[0]].coordinates[0],
                mesh.nodes[vertices[1]].coordinates[0],
                mesh.nodes[vertices[2]].coordinates[0],
                mesh.nodes[vertices[0]].coordinates[1],
                mesh.nodes[vertices[1]].coordinates[1],
                mesh.nodes[vertices[2]].coordinates[1]])