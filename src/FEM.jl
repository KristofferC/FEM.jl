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
export GeoMesh, GeoTrig, GeoQTrig, GeoQuad, GeoNode3, GeoNode2

# Elements, interpolators, materials, materialstatuses
export AbstractFElement, LinTrig, LinQuad, QuadTrig
export AbstractInterpolator, LinTrigInterp, LinQuadInterp, QuadTrigInterp
export AbstractMaterial, LinearIsotropic
export AbstractMaterialStatus, LinearIsotropicMS

export FESection, MaterialSection, ElementSection
export FEProblem
export Stress, Strain
export Solver, NRSolver, solve
export meshquad, gencook
export create_feproblem
export write_VTKXML, VTKExporter, set_binary!, set_compress!

using Logging

@Logging.configure(level=DEBUG, filename="log.log")

include("fields.jl")
include("materials/material.jl")
include("interpolators/interpolator.jl")
include("elements/element.jl")
include("sections.jl")
include("fe_problem.jl")
include("vtkexport_xml.jl")
include("solver.jl")
include("mesh_generators/quad_mesher.jl")
#include("export/vtkexport_jul.jl")

#include("vtkexport.jl")


end
