import FEM.Nvec
import FEM.dNmatrix
import FEM.Jmatrix
import FEM.dNdxmatrix
import FEM.FENode2
import FEM.Point2

facts("FEM.Interpolator") do

context("FEM.Interpolator.LinTrigInterp") do
      interp = LinTrigInterp()

      node_1 = FENode2(1, [1.3, 0.5])
      node_2 = FENode2(2, [1.6, 1.0])
      node_3 = FENode2(3, [0.2, 0.9])

      loc_coords = Point2(0.5, 0.3)

      N = Nvec(interp, loc_coords)

      @fact norm(N - [0.5, 0.3, 1.0 - 0.5 - 0.3]) => roughly(0.0)

end # context

end # facts
