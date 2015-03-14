import FEM.Nvec
import FEM.dNmatrix
import FEM.Jmatrix
import FEM.dNdxmatrix

facts("FEM.Interpolator") do

context("FEM.Interpolator.LinTrigInterp") do
      interp = LinTrigInterp()

      node_1 = Node2([1.3, 0.5], 1)
      node_2 = Node2([1.6, 1.0], 2)
      node_3 = Node2([0.2, 0.9], 3)

      loc_coords = [0.5, 0.3, 0.0]

      N = Nvec(interp, loc_coords)

      @fact norm(N - [0.5, 0.3, 1.0 - 0.5 - 0.3]) => roughly(0.0)

end # context

end # facts
