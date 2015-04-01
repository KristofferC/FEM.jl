module FEM

import Base.push!

import Base.show
import Base.copy

using Base.LinAlg


using Compat
using FixedSizeArrays
using GeometryTypes
#using Devectorize

using Zlib

macro lintpragma(s) end

include("geomesh.jl")
include("core.jl")


export NodeSet, gennodeset, ElementSet, DirichletBC, NodeLoad, Dof
export FENode2, FENode3
export GeoMesh, GeoTrig, GeoQuad, GeoNode3, GeoNode2

# Elements, interpolators, materials, materialstatuses
export AbstractFElement, LinTrig, LinQuad
export AbstractInterpolator, LinTrigInterp, LinQuadInterp
export AbstractMaterial, LinearIsotropic
export AbstractMaterialStatus, LinearIsotropicMS

export FESection, MaterialSection, ElementSection
export FEProblem
export Solver, NRSolver, solve
export meshquad, gencook
export create_feproblem
export exportVTK

#using Logging

#@Logging.configure(level=DEBUG, filename="log.log")



include("sparse.jl")
include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")
include("sections.jl")
include("fe_problem.jl")
include("solver.jl")
include("mesh_generators/quad_mesher.jl")
include("export/vtkexport_jul.jl")
#include("vtkexport.jl")


end
