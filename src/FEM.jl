module FEM

import Base.push!
import Base.show
import Base.copy

using Base.LinAlg

using Compat
using FixedSizeArrays # 0.5 seconds
using Devectorize # 0.8 seconds
using FastAnonymous # 0.04 seconds
using Requires # 0.08 seconds
using Logging # 0.20 seconds

macro lintpragma(s) end


export NodeSet, gennodeset, ElementSet, DirichletBC, NodeLoad, Dof
export FENode2, FENode3
export GeoMesh, GeoTrig, GeoQTrig, GeoQuad, GeoTetra, GeoNode3, GeoNode2

# Elements, interpolators, materials, materialstatuses
export AbstractFElement, LinTrig, LinQuad, QuadTrig, GradTrig, activedofs
export AbstractInterpolator, LinTrigInterp, LinQuadInterp, QuadTrigInterp
export AbstractMaterial, LinearIsotropic, GradMekh
export AbstractElemStorage

export FESection, MaterialSection, ElementSection
export FEProblem
export Stress, Strain, InvFp, KAlpha, VonMises, Kappa, τ
export Solver, NRSolver, solve
export meshquad, gencook, read_mphtxt
export create_feproblem
export AbstractDataExporter, AbstractField


@Logging.configure(level=DEBUG, output = open("FEMlog.log", "w"))


t = @elapsed include("geomesh.jl")
@debug("Time to include geomesh.jl: $t s")

t = @elapsed include("core.jl")
@debug("Time to include core.jl: $t s")

t = @elapsed include("fields.jl")
@debug("Time to include fields.jl: $t s")

t = @elapsed include("materials/material.jl")
@debug("Time to include material.jl: $t s")

t = @elapsed include("interpolators/interpolator.jl")
@debug("Time to include interpolator.jl: $t s")

t = @elapsed include("elements/element.jl")
@debug("Time to include element.jl: $t s")

t = @elapsed include("sections.jl")
@debug("Time to include sections.jl: $t s")

t = @elapsed include("fe_problem.jl")
@debug("Time to include fe_problem.jl: $t s")


include("sparse_tools.jl")

include("export/export.jl")
include("solver.jl")
include("mesh/quad_mesher.jl")
include("mesh/read_mphtxt.jl")
include("mesh/mesh_conversions.jl")

end
