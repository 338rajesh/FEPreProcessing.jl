""" make_unit_cell_FEA_model() -> Dict """

function make_unit_cell_FEA_model(
    ntags::Vector{Int},
    ncoor::Matrix{Float64},
    material_phases_info::Vector{PhaseFiniteElementConnectivity};
    ϵ::Float64 = 1e-06,
    numIP::Int = 4,
    gqOrder::Int = 2,
    verbose::Int = 1,
)::Dict{String, Any}
    verbose > 0 ? printstyled("FINITE ELEMENT PRE-PROCESSING\n"; color=:yellow, bold=true) : nothing
    verbose > 0 ? println(repeat("=", 50)) : nothing
    # =============================================
    #          MAKING FINITE ELEMENT NODE PAIRS
    # =============================================
    uc_node_set::NodeSet = NodeSet(tags=ntags, coor=ncoor, tag="uc_all_nodes")
    verbose > 0 ? println("> Making RVE node pairs,") : nothing
    node_pairs::Vector{NodePair} = make_node_pairs(uc_node_set; small_par=ϵ)
    parent_child_node_tags = [(i.parent_tag, i.child_tag) for i in node_pairs]
    constraining_node_tag = parent_child_node_tags[1][1]  # a_parent_node_tag for constraining DOF 
    verbose > 0 ? println("> Number of parent-child node pairs ", length(parent_child_node_tags), ",") : nothing
    verbose > 0 ? println("> Constraining node tag ", constraining_node_tag, ",") : nothing
    # =============================================
    #          MAKING FINITE ELEMENT NODE SET
    # =============================================
    verbose > 0 ? println("> Making finite element node sets,") : nothing
    all_node_set::Dict{Int, Vector{Float64}} = make_fe_node_set(ntags, ncoor; nodal_dof = ("X", "Y", "Z"))
    # =============================================
    #          MAKING FINITE ELEMENT SET
    # =============================================
    element_sets::Vector{FiniteElementSet} = FiniteElementSet[]
    for a_phase_info in material_phases_info
        verbose > 0 ? println("> Making finite element set with ID: '", a_phase_info.tag, "'") : nothing
        a_phase_element_sets = make_finite_element_sets(
                a_phase_info.ele_conn, all_node_set, a_phase_info.material, a_phase_info.tag;
                num_ip = numIP,
                gq_order = gqOrder,
            )
        append!(element_sets, a_phase_element_sets)
    end
    verbose > 0 ? println(repeat("=", 50)) : nothing
    return Dict{String,Any}(
        "all_node_set" => all_node_set,
        "element_sets" => element_sets,
        "pc_node_tags" => parent_child_node_tags,
        "constraing_ntag" => constraining_node_tag,
    )
end
