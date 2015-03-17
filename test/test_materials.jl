import FEM.stiffness
import FEM.GaussPoint2
import FEM.Point2

facts("FEM.AbstractMaterial") do

context("FEM.Material.LinearIsotropic") do
    mat = LinearIsotropic(200e9, 0.3)
    gp = GaussPoint2(Point2(0.0, 0.0), 0.0)
    context("FEM.Material.LinearIsotropic.stiffness") do

        # Stiffness
        #########################################################################################
        stiff = stiffness(mat, gp)
        # Below calculated using Calfem with command hooke(4, 200e9, 0.3)
        stiff_calfem = [
        [269230769230.769    115384615384.615    115384615384.615                  0];
        [115384615384.615    269230769230.769    115384615384.615                  0];
        [115384615384.615    115384615384.615    269230769230.769                  0];
        [               0                   0                   0   76923076923.0769]]

        @fact norm(stiff - stiff_calfem) / norm(stiff) => roughly(0.0, atol=10.0^(-14))
    end

end # context

end # facts


