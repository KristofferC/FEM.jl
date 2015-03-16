module FEM

import Base.push!

using Compat
using FixedSizeArrays
#using PyCall
#@pyimport vtk

import Base.show
macro lintpragma(s) end




include("geomesh.jl")
export GeoMesh, GeoTrig, GeoQuad, GeoTetra, GeoNode2, GeoNode3

include("core.jl")

export Node2, Node3, NodeSet, gennodeset, ElementSet, DirichletBC, NodeLoad, Du, Dv, Dw, Dof
export Element, LinTrig, BilinQuad
export Interpolator, LinTrigInterp
export Material, LinearIsotropic
export  addnode!, addelem!, addelemset!, addnodeset!, addnodes!, Section
export FEProblem
export Solver, NRSolver, solve
export meshquad, gencook


include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")


include("femesh.jl")
export FEMesh, FENode2, FENode3

#include("fe_problem.jl")
#include("solver.jl")
include("mesh_generators/quad_mesher.jl")
#include("vtkexport.jl")


end
