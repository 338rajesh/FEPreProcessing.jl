get_dielectric_tensor(mat::OrthotropicMaterial)::Matrix{Float64} = [
        mat.eps_ij[:eps_11] 0.0 0.0
        0.0 mat.eps_ij[:eps_22] 0.0
        0.0 0.0 mat.eps_ij[:eps_33]
    ]

get_magnetic_permeability_tensor(mat::OrthotropicMaterial)::Matrix{Float64} = [
    mat.μ_ij[:μ11] 0.0 0.0
    0.0 mat.μ_ij[:μ22] 0.0
    0.0 0.0 mat.μ_ij[:μ33]
]

get_piezoelectric_coupling_tensor(mat::OrthotropicMaterial)::Matrix{Float64} = [
        mat.e_ij[:e11] 0.0 0.0
        mat.e_ij[:e21] 0.0 0.0
        mat.e_ij[:e21] 0.0 0.0
        0.0 0.0 0.0
        0.0 0.0 mat.e_ij[:e15]
        0.0 mat.e_ij[:e15] 0.0
    ]

get_piezomagnetic_coupling_tensor(mat::OrthotropicMaterial)::Matrix{Float64} = [
    mat.q_ij[:q11] 0.0 0.0
    mat.q_ij[:q21] 0.0 0.0
    mat.q_ij[:q21] 0.0 0.0
    0.0 0.0 0.0
    0.0 0.0 mat.q_ij[:q15]
    0.0 mat.q_ij[:q15] 0.0
]


# ==============================================
#           METHODS ON IsotropicMaterial
# ===============================================

function get_elastic_tensor(
    analysis_type::DataType,
    mat::IsotropicMaterial,
)::Matrix{Float64}
    E = mat.E
    nu = mat.nu
    @assert E>0 "Youngs moduli must be > 0."
    kk1::Float64 = E/((1.0+nu)*(1.0-(2.0*nu)))
    kk2::Float64 = 0.5*(1.0 - (2.0 * nu))
    if analysis_type in (Elastic_3DFEA, Thermo_Elastic_3DFEA)
        return kk1*[
            1-nu nu nu 0.0 0.0 0.0
            nu 1-nu nu 0.0 0.0 0.0
            nu nu 1-nu 0.0 0.0 0.0
            0.0 0.0 0.0 kk2 0.0 0.0 
            0.0 0.0 0.0 0.0 kk2 0.0 
            0.0 0.0 0.0 0.0 0.0 kk2
        ]
    elseif analysis_type in (PlaneStrain_2DFEA, Thermo_Elastic_PlaneStrain_2DFEA)
        return kk1 * [
            1-nu nu 0.0
            nu 1-nu 0.0
            0.0 0.0 kk2
        ]
    elseif analysis_type===PlaneStress_2DFEA
        return (E/(1-(nu*nu))) * [
            1 nu 0.0
            nu 1 0.0
            0.0 0.0 0.5*(1-nu)
        ]
    else
        @error "Encountered an unknown analysis type $analysis_type while finding elastic tensor"
    end
end


get_thermal_expansion_vector(mat::IsotropicMaterial) = reshape([mat.alpha mat.alpha mat.alpha 0.0 0.0 0.0], 6, 1)

thermal_expansion_vector_2D(mat::IsotropicMaterial) = reshape([mat.alpha mat.alpha 0.0], 3, 1)

thermal_expansion_vector_plane_strain(mat::IsotropicMaterial) = (1.0 + mat.nu) * reshape([mat.alpha mat.alpha 0.0], 3, 1)

get_thermal_conductivity_tensor(mat::IsotropicMaterial) = [mat.K 0.0 0.0;0.0 mat.K 0.0; 0.0 0.0 mat.K]

# =====================================================
#           METHODS ON TransverselyIsotropicMaterial
# =====================================================

get_thermal_expansion_vector(mat::TransverselyIsotropicMaterial) = reshape([mat.alpha_ij[:alpha_11] mat.alpha_ij[:alpha_22] mat.alpha_ij[:alpha_22] 0.0 0.0 0.0], 6, 1)


# =====================================================
#           GENERAL methods
# =====================================================

function get_material_tensors(
    analysis_type::DataType,
    mat::Material,
)
    if analysis_type in (PlaneStrain_2DFEA, Elastic_3DFEA)
        return get_elastic_tensor(analysis_type, mat)
    elseif analysis_type == Thermo_Elastic_PlaneStrain_2DFEA
        return (
            get_elastic_tensor(analysis_type, mat),
            # thermal_expansion_vector_2D(mat),
            thermal_expansion_vector_plane_strain(mat),
            mat.cv,
        )
    elseif analysis_type==Thermo_Elastic_3DFEA
        return (
            get_elastic_tensor(analysis_type, mat),
            get_thermal_expansion_vector(mat),
            mat.cv,
        )
    elseif analysis_type==Thermal_Conductivity_3DFEA
        return get_thermal_conductivity_tensor(mat)
    elseif analysis_type==Piezo_Electric_3DFEA
        Cmat = get_elastic_tensor(analysis_type, mat)
        e_mat = get_piezoelectric_coupling_tensor(mat)
        eps_mat = get_dielectric_tensor(mat)
        return [
            Cmat -e_mat;
            -transpose(e_mat) -eps_mat;
        ]
    elseif analysis_type==Magneto_Electro_Elastic_3DFEA
        Cmat = get_elastic_tensor(analysis_type, mat)
        e_mat = get_piezoelectric_coupling_tensor(mat)
        q_mat = get_piezomagnetic_coupling_tensor(mat)
        eps_mat = get_dielectric_tensor(mat)
        μ_mat = get_magnetic_permeability_tensor(mat)
        return [
            Cmat -e_mat -q_mat;
            -transpose(e_mat) -eps_mat -zeros(Float64, 3, 3);
            -transpose(q_mat) -zeros(Float64, 3, 3) -μ_mat;
        ]
    else
        @error "Encountered an unknown analysis type while writing material tensors..!" 
    end
end
