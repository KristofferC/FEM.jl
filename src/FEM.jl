module FEM

using Compat
using FixedSizeArrays
#using PyCall
#@pyimport vtk

macro lintpragma(s) end


export Node, NodeSet, gennodeset, ElementSet, GaussPoint, DirichletBC, NodeLoad, Du, Dv, Dw, Dof
export Element, LinTrig, BilinQuad
export Interpolator, LinTrigInterp
export Material, LinearIsotropic
export Mesh, addnode!, addelem!, addelemset!, addnodeset!, addnodes!, Section
export FEProblem
export Solver, NRSolver, solve
export meshquad, gencook

include("core.jl")
include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")
include("mesh.jl")
include("fe_problem.jl")
include("solver.jl")
include("mesh_generators/quad_mesher.jl")
#include("vtkexport.jl")


end
