function unit_cell_FEModel(
    ntags::Vector{Int},
    ncoor::Matrix{Float64},
    phases_data::Dict;
    ϵ::Float64 = 1e-06,
    numIP::Int = 4,
    gqOrder::Int = 2,
)::Dict{String, Any}
    node_pairs = make_rve_node_pairs(ntags, ncoor; small_para = ϵ)
    parent_child_node_tags = [(i.p_nodetag, i.c_nodetag) for i in node_pairs]
    println("Number of parent-child node pairs ", length(parent_child_node_tags))
    constraining_node_tag = parent_child_node_tags[1][1]  # a_parent_node_tag for constraining DOF 
    all_node_set = make_fe_node_set(ntags, ncoor; nodal_dof = ("X", "Y", "Z"))
    #
    element_sets = Dict(
        (
            get_felement_sets(
                phase_ele_conn,
                all_node_set,
                phase_ID;
                num_ip = numIP,
                gq_order = gqOrder,
                print_debug_info = false,
            ),
            phase_material,
        ) for (phase_ID, (phase_ele_conn, phase_material)) in phases_data
    )
    #
    return Dict{String,Any}(
        "all_node_set" => all_node_set,
        "element_sets" => element_sets,
        "pc_node_tags" => parent_child_node_tags,
        "constraing_ntag" => constraining_node_tag,
    )
end
