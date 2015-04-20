const NSLIP = 2

function get_u_dof_idxs(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, nunodes * nuvars)
    count = 1
    for i = 1:nunodes
        for j = 1:nuvars
            if count < nuvars * ngradnodes + 1
                idxs[count] =  (i - 1) * (ngradvars + nuvars) + j
            else
                idxs[count] = count + ngradvars * ngradnodes
            end
            count += 1
        end
    end
    return idxs
end

function get_grad_dof_idxs(nunodes::Int, nuvars::Int, ngradnodes::Int, ngradvars::Int)
    idxs = Array(Int, nunodes * nuvars)
    count = 1
    for i = 1:ngradnodes
        for j = 1:ngradvars
            idxs[count] = i * nuvars + (i - 1) * ngradvars + j
            count += 1
        end
    end
    return idxs
end

function get_grad_idxs_plane(ngradnodes::Int, alpha::Int)
    const NDIM = 2
    idxs = Array(Int, NDIM * ngradnodes)
    count = 1

    for i=1:ngradnodes
        for j=1:NDIM
            idxs[count] = (alpha - 1) * NDIM + (i - 1) * (NDIM * NSLIP) + j;
            count += 1
        end
    end
    return idxs
end


type GradTrigStorage <: ElemStorage
    B::Matrix{Float64}
    Bdiv::Matrix{Float64}
    DeBe::Matrix{Float64}
    Ke::Matrix{Float64}
    ɛ::Vector{Float64}
    f::Vector{Float64}
    f_u::Vector{Float64}
    f_grad::Vector{Float64}
end

function GradTrigStorage()
   B = zeros(4, 12)
   Bdiv = zeros(6)
   DeBe = zeros(4,12)
   Ke = zeros(12 + 2*3*NSLIP, 12 + 2*3*NSLIP)
   ɛ = zeros(4)
   f_u = zeros(12)
   f_grad = zeros(2 * 3 * NSLIP)
   f = zeros(12 + 2*3*NSLIP)
   GradTrigStorage(B, Bdiv, DeBe, Ke, ɛ, f, f_u, f_grad)
end


type GradTrig{T <: AbstractMaterialStatus} <: AbstractFElement{T}
    vertices::Vertex6
    gps::Vector{GaussPoint2}
    n::Int
    interp_u::QuadTrigInterp
    interp_grad::LinTrigInterp
    storage::GradTrigStorage
    matstats::Vector{T}
    temp_matstats::Vector{T}
end
gausspoints(elem::GradTrig) = elem.gps

# Constructor
function GradTrig{T <: AbstractMaterialStatus}(vertices::Vertex6, n, interp_u::QuadTrigInterp,
                interp_grad::LinTrigInterp, storage::GradTrigStorage, gps::Vector{GaussPoint2}, matstat::T)
    matstats = T[]
    temp_matstats = T[]
    for i in 1:length(gps)
        push!(matstats, copy(matstat))
        push!(temp_matstats, copy(matstat))
    end
    GradTrig(vertices, gps, n, interp, storage, matstats, temp_matstats)
end

get_ndofs(::GradTrig) = 12 + 2*3*NSLIP
get_geoelem(ele::GradTrig) = GeoQTrig(ele.n, ele.vertices)
get_geotype(::GradTrig) = GeoQTrig

createstorage(::Type{GradTrig}) = GradTrigStorage()
createinterp(::Type{GradTrig}) = GradTrigInterp()

#=
function creategps(::Type{GradTrig})
    p1 = 1/3
    p2 = 0.2
    p3 = 0.6
    [GaussPoint2(Point2(p1, p1), -0.281250000000000);
     GaussPoint2(Point2(p2, p3), 0.260416666666667);
     GaussPoint2(Point2(p2, p2), 0.260416666666667)
     GaussPoint2(Point2(p3, p2), 0.260416666666667)]
end
=#

function creategps(::Type{GradTrig})
    p1 = 2/3
    p2 = 1/6
    w = 1/3
    [GaussPoint2(Point2(p1, p2), w);
     GaussPoint2(Point2(p2, p1), w);
     GaussPoint2(Point2(p2, p2), w)]
end

@inline function get_ndofs(::GradTrig)
    return 12 + 2*3*NSLIP
end

function doftypes(::GradTrig, v::Int)
    if v <= 3
        return [Du, Dv, Gu1, Gv1, Gu2, Gv2] # TODO: Fix from arbitrary NSLIP
    else
        return [Du, Dv]
    end
end

#=
computeKappa(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)



    FloatArray g, para = mekhmat->givePara();
    double l = para.at(4);
    double Hg = para.at(7);
    double factor = para.at(10);



    this->getOverrideInterp()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(elem) );

    computeSlipDegreesOfFreedom(g, stepN);

    int nDim = 3;

    FloatMatrix hessian(nDim, nDim);
    for ( int i = 1; i <= this->getNslipPlanes(); i++ ) {
        hessian.zero();

        // Take out the slip dofse for the current slip plane
        FloatArray g_alpha;
        IntArray locK_alpha;
        setSlipDofsForSlipPlaneArray(locK_alpha, i, stepN);
        g_alpha.beSubArrayOf(g, locK_alpha);

        // Calculate hessian
        // Rewrite this to matrix multiplication?
        for ( int j = 1; j <= nDim; j++ ) {
            for ( int k = 1; k <= nDim; k++ ) {
                for ( int node = 1; node <= this->getNslipNodes(); node++ ) {
                    hessian.at(j, k) += dnx.at(node, k) * g_alpha.at( ( node - 1 ) * nDim + j );
                }
            }
        }

        // Isotropic contribution
        answer.at(i) = Hg * l * l * factor * hessian.giveTrace();
    }
=#

function stiffness{P <: AbstractMaterial}(elem::GradTrig,
                                          nodes::Vector{FENode2},
                                          material::P)
    const H = 10e-7
    Ke = elem.storage.Ke
    fill!(Ke, 0.0)
    f = intf(elem, mat, nodes)
    f = copy(f)
    col = 1
    for node in nodes
        for dof in node.dofs
            if !dof.active
                continue
            end
            dof.value += H
            f_pert = intf(elem, mat, nodes)
            @devec Ke[:, col] = (f_pert .- f) ./ H
            dof.value -= H
            col += 1
        end
    end
    return Ke
end

function intf_u{P <: AbstractMaterial}(elem::GradTrig, mat::P, nodes::Vector{FENode2})
    u = get_field(elem, nodes)
    uf = u[get_u_dof_idxs(6, 2, 3, 4)]
    ɛ = elem.storage.ɛ
    fill!(elem.storage.f_u, 0.0)
    for (i, gp) in enumerate(elem.gps)
        B = Bmatrix(elem, gp, nodes)
        A_mul_B!(ɛ, B, u)
        fill_from_start!(elem.temp_matstats[i].strain, ɛ)

        F_grad_mekh = [ɛ[1] + 1.0, ɛ[2] + 1.0,  1.0,
                       ɛ[4], 0.0, 0.0, 0.0, ɛ[4], 0.0]

        σ = stress(mat, F_grad_mekh, matstat[i])
        @debug("σ = $σ")
        fill_from_start!(elem.temp_matstats[i].stress, σ)

        dV = weight(elem, gp, nodes)

        # f_u += B' * σ * dV
        BLAS.gemv!('T', dV, B, σ, 1.0, elem.storage.f_u)
    end
    return elem.storage.f_u
end


function intf_grad{P <: AbstractMaterial}(elem::GradTrig, mat::P, nodes::Vector{FENode2})
    fill!(f_grad, 0.0)
    u = get_field(elem, nodes)
    M = mass_matrix(elem.interp_grad, nodes)

    for slip in 1:NSLIP
        idxs = get_grad_idxs_plane(6, 2, 3, 4, slip)
        fe[idxs] += M * u[idxs]
    end

    # Dummy gp, constant
    B = Bdiv(elem, elem.gps[1], nodes)
    dN = dNmatrix(elem.interp_u, gp.local_coords)
    J = Jmatrix(elem.interp_u, gp.local_coords, elem.vertices, nodes, dN)

    for i in 1:NSLIP
        k_alpha_tot = 0.0
        for i in length(gausspoints)
            k_alpha_tot += elem.matstats[i].k_alpha
        end
        k_alpha /= length(gausspoints)

        fe += k_alpha * B * A0
    end
end


function intf{P <: AbstractMaterial}(elem::GradTrig, mat::P, nodes::Vector{FENode2})
 #   f_u = intf_u(elem, mat, nodes)
  #  f_grad = intf_grad(elem, matm, nodes)
   # f[] = f_u
    #f[] = f_grad
    return
end



function Bdiv(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp_grad, gp.local_coords, elem.vertices, nodes)
    Bdiv = elem.storage.Bdiv
    Bdiv[1] = dNdx[1, 1]
    Bdiv[2] = dNdx[1, 2]
    Bdiv[3] = dNdx[2, 1]
    Bdiv[4] = dNdx[2, 2]
    Bdiv[5] = dNdx[3, 1]
    Bdiv[6] = dNdx[3, 2]
    return Bdiv
end


function Bmatrix(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dNdx = dNdxmatrix(elem.interp_u, gp.local_coords, elem.vertices, nodes)
    B = elem.storage.B
    for i in 1:6
        B[1, 2*i - 1] = dNdx[i, 1]
        B[2, 2*i]     = dNdx[i, 2]
        B[3, 2*i - 1] = dNdx[i, 2]
        B[3, 2*i]     = dNdx[i, 1]
    end
    return B
end


function weight(elem::GradTrig, gp::GaussPoint2, nodes::Vector{FENode2})
    dN = dNmatrix(elem.interp, gp.local_coords)
    J = Jmatrix(elem.interp, gp.local_coords, elem.vertices, nodes, dN)
    return abs(det2x2(J)) * gp.weight
end

# Get the stress/strain in gausspoint i
get_field(elem::GradTrig, ::Type{Stress}, i::Int) = elem.matstats[i].stress
get_field(elem::GradTrig, ::Type{Strain}, i::Int) = elem.matstats[i].strain


