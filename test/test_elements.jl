import FEM.stiffness
import FEM.Bmatrix
import FEM.get_gps
import FEM.get_interp
import FEM.get_storage
import FEM.Vertex3
import FEM.Vertex4
facts("FEM.Element") do


context("FEM.Element.LinTrig") do


    nodes = [FENode2(1, [0, 0]), FENode2(2, [1, 1]), FENode2(3, [1, 2])]

    # Create the gps, interpolator, storage...

    gps = get_gps(LinTrig)
    interp = get_interp(LinTrig)
    storage = get_storage(LinTrig)
    elem = LinTrig(Vertex3(1, 2, 3), gps, 1, interp, storage)
    mat = LinearIsotropic(200e9, 0.3)

    Ke = stiffness(elem, nodes, mat)

    # Calculated with plante([0.0, 1.0, 1.0], [0.0, 1.0, 2.0], [2, 1], hooke(2, 200e9, 0.3))
    Ke_calfem =  1.0e+11 * [
    [ 1.346153846153846                   0  -2.692307692307692   0.576923076923077   1.346153846153846  -0.576923076923077];
    [                 0   0.384615384615385   0.384615384615385  -0.769230769230769  -0.384615384615385   0.384615384615385];
    [-2.692307692307692   0.384615384615385   5.769230769230769  -1.923076923076923  -3.076923076923077   1.538461538461538];
    [ 0.576923076923077  -0.769230769230769  -1.923076923076923   2.884615384615385   1.346153846153846  -2.115384615384615];
    [ 1.346153846153846  -0.384615384615385  -3.076923076923077   1.346153846153846   1.730769230769231  -0.961538461538461];
    [-0.576923076923077   0.384615384615385   1.538461538461538  -2.115384615384615  -0.961538461538461   1.730769230769231]]

    @fact norm(Ke - Ke_calfem) / norm(Ke) => roughly(0.0)

end # context


context("FEM.Element.LinQuad") do

    nodes = [FENode2(1, [0.0, 0.0]), FENode2(2, [1.0, 0.0]), FENode2(3, [1.0, 2.5]), FENode2(4, [0.0, 1.5])]

    gps = get_gps(LinQuad)
    interp = get_interp(LinQuad)
    storage = get_storage(LinQuad)
    elem = LinQuad(Vertex4(1, 2, 3, 4), gps, 1, interp, storage)
    mat = LinearIsotropic(200e9, 0.3)



    elem = LinQuad(Vertex4(1, 2, 3, 4), gps, 1, interp, storage)

    Ke = stiffness(elem, nodes, mat)

    # Calculated with plani4e([0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 2.5, 1.5], [2, 1, 2], hooke(2, 200e9, 0.3))
    Ke_calfem = 1.0e+11 * [
    [ 1.665302782324059   0.296644844517185  -1.653027823240589   0.014320785597381  -0.366202945990180  -0.398936170212766   0.353927986906710   0.087970540098200];
    [ 0.296644844517185   0.949263502454992  -0.177986906710311  -0.261865793780688  -0.398936170212766  -0.315057283142389   0.280278232405892  -0.372340425531915];
    [-1.653027823240589  -0.177986906710311   2.491816693944353  -0.623977086743044   0.873567921440262   0.239361702127660  -1.712356792144026   0.562602291325696];
    [ 0.014320785597381  -0.261865793780687  -0.623977086743044   1.080196399345336   0.047054009819967  -0.118657937806874   0.562602291325696  -0.699672667757774];
    [-0.366202945990180  -0.398936170212766   0.873567921440262   0.047054009819967   1.145662847790507   0.337561374795417  -1.653027823240589   0.014320785597381];
    [-0.398936170212766  -0.315057283142389   0.239361702127660  -0.118657937806874   0.337561374795417   0.695581014729951  -0.177986906710311  -0.261865793780688];
    [ 0.353927986906710   0.280278232405892  -1.712356792144026   0.562602291325696  -1.653027823240589  -0.177986906710311   3.011456628477905  -0.664893617021277];
    [ 0.087970540098200  -0.372340425531915   0.562602291325696  -0.699672667757774   0.014320785597381  -0.261865793780687  -0.664893617021277   1.333878887070377]]

    @fact norm(Ke - Ke_calfem) / norm(Ke) => roughly(0.0, 10.0^(-15))



end # context


end # facts
