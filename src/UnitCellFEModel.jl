"""
    make_unit_cell_FEA_model() -> Dict

This function does, the following in order,

 - Node pairing (from node tags and their coordinate data)
 - 

"""
function make_unit_cell_FEA_model(
    ntags::Vector{Int},
    ncoor::Matrix{Float64},
    material_phases_info::Vector{PhaseFiniteElementConnectivity};
    ϵ::Float64 = 1e-06,
    numIP::Int = 4,
    gqOrder::Int = 2,
)::Dict{String, Any}
    #
    # =============================================
    #          MAKING FINITE ELEMENT NODE PAIRS
    # =============================================
    node_pairs = make_rve_node_pairs(ntags, ncoor; small_para = ϵ)
    parent_child_node_tags = [(i.p_nodetag, i.c_nodetag) for i in node_pairs]
    println("Number of parent-child node pairs ", length(parent_child_node_tags))
    constraining_node_tag = parent_child_node_tags[1][1]  # a_parent_node_tag for constraining DOF 
    #
    # =============================================
    #          MAKING FINITE ELEMENT NODE SET
    # =============================================
    all_node_set::Dict{Int, Vector{Float64}} = make_fe_node_set(
        ntags,
        ncoor;
        nodal_dof = ("X", "Y", "Z")
    )
    #
    # =============================================
    #          MAKING FINITE ELEMENT SET
    # =============================================

    element_sets::Vector{FiniteElementSet} = FiniteElementSet[]

    for a_phase_info in material_phases_info
        a_phase_element_sets = make_finite_element_sets(
                a_phase_info.ele_conn,
                all_node_set,
                a_phase_info.material,
                a_phase_info.tag;
                num_ip = numIP,
                gq_order = gqOrder,
            )
        append!(element_sets, a_phase_element_sets)
    end

    return Dict{String,Any}(
        "all_node_set" => all_node_set,
        "element_sets" => element_sets,
        "pc_node_tags" => parent_child_node_tags,
        "constraing_ntag" => constraining_node_tag,
    )
end
