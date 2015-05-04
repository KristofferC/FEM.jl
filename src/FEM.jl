module FEM

import Base.push!

import Base.show
import Base.copy

using Base.LinAlg

using Compat
using FixedSizeArrays

using Devectorize
using FastAnonymous

macro lintpragma(s) end




export NodeSet, gennodeset, ElementSet, DirichletBC, NodeLoad, Dof
export FENode2, FENode3
export GeoMesh, GeoTrig, GeoQTrig, GeoQuad, GeoNode3, GeoNode2

# Elements, interpolators, materials, materialstatuses
export AbstractFElement, LinTrig, LinQuad, QuadTrig, GradTrig, activedofs
export AbstractInterpolator, LinTrigInterp, LinQuadInterp, QuadTrigInterp
export AbstractMaterial, LinearIsotropic, GradMekh

export FESection, MaterialSection, ElementSection
export FEProblem
export Stress, Strain, InvFp, KAlpha, VonMises
export Solver, NRSolver, solve
export meshquad, gencook
export create_feproblem
export write_data, VTKExporter, set_binary!, set_compress!

using Logging

@Logging.configure(level=CRITICAL, filename="log.log")

include("geomesh.jl")
include("core.jl")

include("fields.jl")
include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")
include("sections.jl")
include("fe_problem.jl")
include("vtkexport_xml.jl")
include("sparse_tools.jl")
include("solver.jl")
include("mesh_generators/quad_mesher.jl")

#include("export/vtkexport_jul.jl")

#include("vtkexport.jl")


end
