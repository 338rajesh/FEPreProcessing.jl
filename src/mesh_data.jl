# ==================================================================================#
#
#                       Prepare Finite Element Sets
#
# ==================================================================================#

"""
    get_felement_sets() -> Vector{AbstractFElementSet}
    
 Args:
 -----
 - **element_connectivity**::Dict{Int64, Matrix{Int64}},
     Element connectivity gmsh ids are used to identify element sets which will eventually be converted to ABAQUS convention of element naming.
 - **nodal_data**::Dict{Int64, Tuple{Vararg{Float64}}},    
 - **elset_tag**::Union{Int64, String};
 - **num_ip**::Int64=-1,
 - **gq_order**::Int64=-1,
 - **print_debug_info**::Bool=false,

"""


function make_finite_element_set(
    element_type::Int,
    element_connectivity::Matrix{Int},
    nodal_data::Dict{Int,Vector{Float64}};
    elset_tag::Union{String,Int}="Element-Set",
    material::Union{Nothing,Material}=nothing,
    num_ip::Int64=-1,
    gq_order::Int64=-1
)::FiniteElementSet
    #
    elements::Vector{AbstractFElement} = AbstractFElement[]
    ivols::Vector{Float64} = Float64[]
    shfg_xyz::Vector{SMatrix} = SMatrix[]
    num_node_per_fele::Int = size(element_connectivity, 1) - 1  # first row is for ele_tag => size-1
    # num_fele::Int = size(element_connectivity, 2)
    num_shape_func::Int = num_node_per_fele
    dim, fele_data_type, integration_points = get_felement_properties(
        element_type, gq_order, num_ip
    )
    parent_csys_shf_data = get_shape_functions_info(
        fele_data_type,
        integration_points
    )
    #
    # Running over each element of the connectivity table.
    for a_ele_info in eachcol(element_connectivity)
        empty!(ivols)
        empty!(shfg_xyz)
        ele_tag = a_ele_info[1]
        ele_nodes = a_ele_info[2:end]
        nodal_coordinates = [nodal_data[i][1:dim] for i in ele_nodes]
        for (a_ip, shfg_parent) in parent_csys_shf_data
            Jacobian = zeros(Float64, dim, dim)
            # run over pairs of coordinates and shape function data
            for (k2, a_node_coor) in enumerate(nodal_coordinates)
                Jacobian +=
                    (shfg_parent[2:end, k2:k2] * reshape(collect(a_node_coor), 1, :))
            end
            Jacobian_det = det(Jacobian)
            @assert Jacobian_det > 0.0 "Determinant of a Jacobian must be >= 0.0 but the element $(ele_tag) has $Jacobian_det."
            push!(ivols, Jacobian_det * a_ip.wt)
            aip_shf_grad_xyz =
                [shfg_parent[1:1, :]; inv(Jacobian) * shfg_parent[2:end, :]]
            push!(shfg_xyz, SMatrix{1 + dim,num_shape_func,Float64}(aip_shf_grad_xyz))
        end
        push!(
            elements,
            fele_data_type(
                ele_tag,
                SVector{num_node_per_fele,Int}(ele_nodes),
                tuple(ivols...),
                tuple(shfg_xyz...),
            )
        )  # end of updating elements vector.
    end
    return FiniteElementSet(
        elset_tag,
        elements,
        material,
    )
end


function make_finite_element_sets(
    element_connectivities::Dict{Int64,Matrix{Int64}},
    nodal_data::Dict{Int64,Vector{Float64}},
    ele_material::Material,
    elset_tag::Union{Int64,String};
    num_ip::Int64=-1,
    gq_order::Int64=-1
)::Vector{FiniteElementSet}
    #
    return [make_finite_element_set(
        a_ele_type,
        a_ele_type_ele_connectivity,
        nodal_data;
        elset_tag=elset_tag,
        material=ele_material,
        num_ip=num_ip,
        gq_order=gq_order
    ) for (a_ele_type, a_ele_type_ele_connectivity) in element_connectivities]
end



# =================================================================================
#                       Prepare Node Sets
# =================================================================================

"""
It returns a dictionary with the following key-value pairs

key:
---- 
    Node tag
Value:
-----
    A tuple of nodal coordinates. If not specified, it will have x, y, z coordinates but one can choose to have either (x, y) or (y, z) or (z, x) or (x, ) or (y, ) or (z, ) using ** option

"""
function make_fe_node_set(
    node_tags::Union{Vector{Int64},Vector{Int32}},
    nodal_coor::Matrix{Float64};
    nodal_dof::Tuple{Vararg{String}}=("x", "y", "z")
)::Dict{Int64,Vector{Float64}}  # TODO make node sets object oriented
    #
    @assert size(nodal_coor, 1) >= 2 "Number of rows must be 2 or 3 being in the order of x, y, z"
    @assert length(node_tags) == size(nodal_coor, 2) "Number of node tags and number of columns in nodal coordinates matrix doesn't match; $(length(node_tags)) ≠ $(size(nodal_coor)[1]))"
    #
    req_dofs::Vector{Int64} = []
    for a_dof in nodal_dof
        if uppercase(a_dof) == "X"
            push!(req_dofs, 1)
        elseif uppercase(a_dof) == "Y"
            push!(req_dofs, 2)
        elseif uppercase(a_dof) == "Z"
            push!(req_dofs, 3)
        end
    end
    #
    nodal_coor = nodal_coor[req_dofs, :]  # extracting dof, as required.
    #
    node_tags = convert(Vector{Int64}, node_tags)
    #
    return Dict((i, j) for (i, j) in zip(node_tags, eachcol(nodal_coor)))
end

# =================================================================================
#                       Misc.
# =================================================================================

function get_total_volume_of_elements(
    finite_element_sets::Vector{FiniteElementSet};
    print_info::Bool=false
)
    ele_set_volumes::Vector{Float64} = Float64[]
    element_counter::Int = 0
    a_ele_set_vol::Float64 = 0.0
    if print_info
        println("\n EVALUATING VOLUMES OF ELEMENT SETS...!")
    end
    for a_finite_element_set in finite_element_sets
        a_ele_set_vol = 0.0
        for a_finte_element in a_finite_element_set.elements
            a_ele_set_vol += sum(a_finte_element.ip_volumes)
            element_counter += 1
        end
        if print_info
            println(
                "Volume of element set with ID: ", a_finite_element_set.tag, " is ", a_ele_set_vol
            )
        end
        push!(ele_set_volumes, a_ele_set_vol)
    end
    total_volume::Float64 = sum(ele_set_volumes)
    if print_info
        println("Sum of volumes of ", element_counter, " elements is  ", total_volume)
    end
    return total_volume
end



