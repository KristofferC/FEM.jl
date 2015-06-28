abstract AbstractMaterial

abstract AbstractMaterialStatus


stiffness() = error("Not implemented")
create_matstat() = error("Not implemented")
stress() = error("Not implemented")
get_kalpha() = error("Not implemented")
get_kalphas() = error("Not implemented")

@lazymod LinearIsotropicMod "linear_isotropic.jl"