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


"""
    preapre_constraints_matrix_csc() -> SparseMatrixCSC

# Args
- **node_tags** :: Vector{T}
- **parent_child_node_tag_pairs** :: Vector{Tuple{T, T}}
- **dof_per_node** :: T

"""
function preapre_constraints_matrix_csc(
    analysis_type::DataType,
    node_tags::Vector{T},
    parent_child_node_tag_pairs::Vector{Tuple{T, T}};
    hdbc::Dict{T, Vector{Symbol}}=Dict{T, Vector{Symbol}}(),
)::SparseMatrixCSC{Float64, T} where {T}    
    # ==============================
    # term: meaning
    # -------------
    # rnt: retaining node tag
    # pnt: parent node tag
    # cnt: child node tag
    # hdbc: homogenous Dirichlet boundary conditions
    # ==============================
    println("Preparing CSC type matrix for eliminating child node tags..!")
    dof_quantities = get_nodal_dof(analysis_type)
    dof_per_node = length(dof_quantities)
    child_node_tags::Vector{T} = [i[2] for i in parent_child_node_tag_pairs]
    retaining_node_tags::Vector{T} = [i for i in node_tags if !(i in child_node_tags)]
    total_dof::T = dof_per_node*length(node_tags)
    retaining_dof::T = dof_per_node * length(retaining_node_tags)
    #
    row_indices::Vector{Int64} = T[]
    col_indices::Vector{Int64} = T[]
    for ant in node_tags
        append!(row_indices, (1+((ant-1) * dof_per_node)):(ant*dof_per_node))
        if !(ant in child_node_tags)
            a_rnt_index = findfirst(isequal(ant), retaining_node_tags)
        else
            a_cnt_index = findfirst(isequal(ant), child_node_tags)
            a_pnt = parent_child_node_tag_pairs[a_cnt_index][1]
            a_rnt_index = findfirst(isequal(a_pnt), retaining_node_tags)
        end
        append!(col_indices, (1+((a_rnt_index-1) * dof_per_node)):(a_rnt_index*dof_per_node))
    end
    if length(hdbc) > 0
        #
        # prepare columns where hdbc needs to be imposed.
        rm_column_indices::Vector{T} = T[]
        for (a_hdbc_nt, spec_dof) in hdbc
            a_hdbc_nt_index = findfirst(isequal(a_hdbc_nt), retaining_node_tags)
            if length(spec_dof)==1 && spec_dof[1]==:ALL
                dof_indices = collect(1:dof_per_node)
            else
                dof_indices = [kount for (kount, i) in enumerate(dof_quantities) if i in spec_dof]
            end
            append!(
                rm_column_indices,
                ((1+((a_hdbc_nt_index-1)*dof_per_node)):(a_hdbc_nt_index*dof_per_node))[dof_indices]
            )
        end
        # remove columns
        #
        row_indices_1::Vector{T} = T[]
        col_indices_1::Vector{T} = T[]
        for (ii, jj) in zip(row_indices, col_indices)
            if !(jj in rm_column_indices)
                # reduce jj by jj1 where jj1 is the number of elements in rm_column_indices which are < jj
                jj1 = length(findall(rm_column_indices .< jj))
                push!(row_indices_1, ii)
                push!(col_indices_1, jj - jj1)
            end
        end
        #
        retaining_dof -= length(rm_column_indices)
        #
        return sparse(row_indices_1, col_indices_1, 1.0, total_dof, retaining_dof)
    else
        return sparse(row_indices, col_indices, 1.0, total_dof, retaining_dof)
    end
end


