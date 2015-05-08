abstract AbstractMaterial

abstract AbstractMaterialStatus


stiffness() = error("Not implemented")
create_matstat() = error("Not implemented")
stress() = error("Not implemented")

@lazymod LinearIsotropicMod "linear_isotropic.jl"
@lazymod GradMekhMod "grad_mekh.jl"