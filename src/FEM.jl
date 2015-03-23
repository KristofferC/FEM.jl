module FEM

import Base.push!

using Compat
using FixedSizeArrays
#using PyCall
#@pyimport vtk

import Base.show
macro lintpragma(s) end

include("geomesh.jl")


include("core.jl")

export NodeSet, gennodeset, ElementSet, DirichletBC, NodeLoad, Du, Dv, Dw, Dof
export FEMesh, FENode2, FENode3
export GeoMesh, GeoTrig, GeoQuad, GeoNode3, GeoNode2
export AbstractFElement, LinTrig #, BilinQuad
export AbstractInterpolator, LinTrigInterp
export AbstractMaterial, LinearIsotropic
export AbstractMaterialStatus, LinearIsotropicMS
export FESection, MaterialSection, ElementSection
export FEProblem
export Solver, NRSolver, solve
export meshquad, gencook
export create_feproblem


include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")

include("femesh.jl")

include("fe_problem.jl")
include("solver.jl")
include("mesh_generators/quad_mesher.jl")
#include("vtkexport.jl")


end
