abstract AbstractMaterial

abstract AbstractMaterialStatus


stiffness() = error("Not implemented")
create_matstat() = error("Not implemented")
stress() = error("Not implemented")
get_kalpha() = error("Not implemented")
get_kalphas() = error("Not implemented")

@lazymod LinearIsotropicMod "linear_isotropic.jl"
@lazymod GradMekhMod "grad_mekh.jl"
@lazymod GradMekhModJl "grad_mekh_jul.jl"
@lazymod GradMekhModJlSmall "grad_mekh_jul_small.jl"
