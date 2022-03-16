"""
for a given type of analysis
    Get element stiffness matrix formulation
        Get B matrix


"""

function get_strain_disp_matrix_3D(
    shfg_xyz::Union{SMatrix,Matrix},
)::Matrix{Float64}
    num_shape_func::Int64 = size(shfg_xyz)[2]
    Bmatrix::Matrix{Float64} = zeros(Float64, 6, 3 * num_shape_func)
    for i in 1:num_shape_func
        dNdX, dNdY, dNdZ = shfg_xyz[:, i]
        Bmatrix[:, ((3*i)-2):(3*i)] = [
            dNdZ 0.0 0.0
            0.0 dNdX 0.0
            0.0 0.0 dNdY
            0.0 dNdY dNdX
            dNdY 0.0 dNdZ
            dNdX dNdZ 0.0
        ]
    end
    return Bmatrix
end


function get_piezo_Bmat(
    shfg_xyz::Union{SMatrix,Matrix},
)::Matrix{Float64}
    num_shape_func::Int64 = size(shfg_xyz)[2]
    Bmatrix::Matrix{Float64} = zeros(Float64, 9, 4 * num_shape_func)
    for i in 1:num_shape_func
        dNdX, dNdY, dNdZ = shfg_xyz[:, i]
        Bmatrix[:, ((4*i)-3):(4*i)] = [
            dNdZ 0.0 0.0 0.0
            0.0 dNdX 0.0 0.0
            0.0 0.0 dNdY 0.0
            0.0 dNdY dNdX 0.0
            dNdY 0.0 dNdZ 0.0
            dNdX dNdZ 0.0 0.0
            0.0 0.0 0.0 -dNdZ
            0.0 0.0 0.0 -dNdX
            0.0 0.0 0.0 -dNdY
        ]
    end
    return Bmatrix
end


function get_MEE_Bmat(
    shfg_xyz::Union{SMatrix,Matrix},
)::Matrix{Float64}
    num_shape_func::Int64 = size(shfg_xyz)[2]
    Bmatrix::Matrix{Float64} = zeros(Float64, 12, 5 * num_shape_func)
    for i in 1:num_shape_func
        dNdX, dNdY, dNdZ = shfg_xyz[:, i]
        Bmatrix[:, ((5*i)-4):(5*i)] = [
            dNdZ 0.0 0.0 0.0 0.0
            0.0 dNdX 0.0 0.0 0.0
            0.0 0.0 dNdY 0.0 0.0
            0.0 dNdY dNdX 0.0 0.0
            dNdY 0.0 dNdZ 0.0 0.0
            dNdX dNdZ 0.0 0.0 0.0
            0.0 0.0 0.0 -dNdZ 0.0
            0.0 0.0 0.0 -dNdX 0.0
            0.0 0.0 0.0 -dNdY 0.0
            0.0 0.0 0.0 0.0 -dNdZ
            0.0 0.0 0.0 0.0 -dNdX
            0.0 0.0 0.0 0.0 -dNdY
        ]
    end
    return Bmatrix
end


function get_strain_disp_matrix_2D(
    shfg_xyz::Union{SMatrix,Matrix}
)::Matrix{Float64}
    num_shape_func::Int64 = size(shfg_xyz)[2]
    Bmatrix::Matrix{Float64} = zeros(Float64, 3, 2 * num_shape_func)
    for i in 1:num_shape_func
        dNdX, dNdY = shfg_xyz[:, i]
        Bmatrix[:, ((2*i)-1):(2*i)] = [
            dNdX 0.0
            0.0 dNdY
            dNdY dNdX
        ]
    end
    return Bmatrix
end


function get_thermal_conductivity_Bmat(
    shfg_xyz::Union{SMatrix,Matrix}
)::Matrix{Float64}
    num_shape_func::Int64 = size(shfg_xyz)[2]
    Bmatrix::Matrix{Float64} = zeros(Float64, 3, num_shape_func)
    for i in 1:num_shape_func
        dNdX, dNdY, dNdZ = shfg_xyz[:, i]
        Bmatrix[:, i:i] = [dNdZ, dNdX, dNdY ]
    end
    return Bmatrix
end
 

function get_B_matrix(
    shf_grad_xyz::Union{SMatrix,Matrix},
    analysis_type::DataType
)::Matrix{Float64}
    if analysis_type in (PlaneStress_2DFEA, PlaneStrain_2DFEA)
        return get_strain_disp_matrix_2D(shf_grad_xyz)
    elseif analysis_type in (Elastic_3DFEA, Thermo_Elastic_3DFEA, )
        return get_strain_disp_matrix_3D(shf_grad_xyz)
    elseif analysis_type === Piezo_Electric_3DFEA
        return get_piezo_Bmat(shf_grad_xyz)
    elseif analysis_type in (Thermal_Conductivity_3DFEA, )
        return get_thermal_conductivity_Bmat(shf_grad_xyz)
    elseif analysis_type in (Magneto_Electro_Elastic_3DFEA, )
        return get_MEE_Bmat(shf_grad_xyz)
    else
        @error "Encountered unknown analysis type :: $analysis_type, in B matrix evaluation."
    end
end