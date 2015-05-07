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
export Stress, Strain, InvFp, KAlpha, VonMises
export Solver, NRSolver, solve
export meshquad, gencook, read_mphtxt
export create_feproblem
export AbstractDataExporter, AbstractField
#export write_data, VTKExporter, set_binary!, set_compress!

@Logging.configure(level=CRITICAL, filename="log.log")

include("geomesh.jl")
include("core.jl")

include("fields.jl")
include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")
include("sections.jl")
include("fe_problem.jl")

include("sparse_tools.jl")

include("export/export.jl")
include("solver.jl")
include("mesh/quad_mesher.jl")
include("mesh/read_mphtxt.jl")
include("mesh/mesh_conversions.jl")



#include("vtkexport.jl")


end
