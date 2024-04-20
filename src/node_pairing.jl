
function get_boundary_node_sets(
    node_set::NodeSet,
    bbox::NTuple{6,Float64}; 
    small_par::Float64=1e-06, 
)::Vector{NodeSet}
    n_tags = node_set.tags
    n_coor = node_set.coor
    x_min, y_min, z_min, x_max, y_max, z_max = bbox
    x::Vector{Float64} = node_set.coor[1, :]
    y::Vector{Float64} = node_set.coor[2, :]
    z::Vector{Float64} = node_set.coor[3, :]
    xn_mask::BitVector = x .<= x_min + small_par
    yn_mask::BitVector = y .<= y_min + small_par
    zn_mask::BitVector = z .<= z_min + small_par
    xp_mask::BitVector = x .>= x_max - small_par
    yp_mask::BitVector = y .>= y_max - small_par
    zp_mask::BitVector = z .>= z_max - small_par

    # Nodes exclusively on the vertices
    xn_yn_zn::BitVector = xn_mask .& yn_mask .& zn_mask
    xp_yn_zn::BitVector = xp_mask .& yn_mask .& zn_mask
    xp_yp_zn::BitVector = xp_mask .& yp_mask .& zn_mask
    xn_yp_zn::BitVector = xn_mask .& yp_mask .& zn_mask
    xn_yn_zp::BitVector = xn_mask .& yn_mask .& zp_mask
    xp_yn_zp::BitVector = xp_mask .& yn_mask .& zp_mask
    xp_yp_zp::BitVector = xp_mask .& yp_mask .& zp_mask
    xn_yp_zp::BitVector = xn_mask .& yp_mask .& zp_mask

    # Nodes exclusively on the edges
    xn_yn::BitVector = xn_mask .& yn_mask .& .!xn_yn_zn .& .!xn_yn_zp
    xp_yn::BitVector = xp_mask .& yn_mask .& .!xp_yn_zn .& .!xp_yn_zp
    xp_yp::BitVector = xp_mask .& yp_mask .& .!xp_yp_zn .& .!xp_yp_zp
    xn_yp::BitVector = xn_mask .& yp_mask .& .!xn_yp_zn .& .!xn_yp_zp
    # 
    xn_zn::BitVector = xn_mask .& zn_mask .& .!xn_yn_zn .& .!xn_yp_zn
    xp_zn::BitVector = xp_mask .& zn_mask .& .!xp_yn_zn .& .!xp_yp_zn
    xp_zp::BitVector = xp_mask .& zp_mask .& .!xp_yn_zp .& .!xp_yp_zp
    xn_zp::BitVector = xn_mask .& zp_mask .& .!xn_yn_zp .& .!xn_yp_zp
    # 
    yn_zn::BitVector = yn_mask .& zn_mask .& .!xn_yn_zn .& .!xp_yn_zn
    yp_zn::BitVector = yp_mask .& zn_mask .& .!xp_yp_zn .& .!xn_yp_zn
    yp_zp::BitVector = yp_mask .& zp_mask .& .!xp_yp_zp .& .!xn_yp_zp
    yn_zp::BitVector = yn_mask .& zp_mask .& .!xn_yn_zp .& .!xp_yn_zp

    # Nodes exclusively on the faces
    xn_face_mask::BitVector = xn_mask .& .!yn_mask .& .!zn_mask .& .!yp_mask .& .!zp_mask
    yn_face_mask::BitVector = .!xn_mask .& yn_mask .& .!zn_mask .& .!xp_mask .& .!zp_mask
    zn_face_mask::BitVector = .!xn_mask .& .!yn_mask .& zn_mask .& .!xp_mask .& .!yp_mask
    xp_face_mask::BitVector = .!yn_mask .& .!zn_mask .& xp_mask .& .!yp_mask .& .!zp_mask
    yp_face_mask::BitVector = .!xn_mask .& .!zn_mask .& .!xp_mask .& yp_mask .& .!zp_mask
    zp_face_mask::BitVector = .!xn_mask .& .!yn_mask .& .!xp_mask .& .!yp_mask .& zp_mask

    make_node_set(bv::BitVector, t::String)::NodeSet = NodeSet(n_tags[bv], n_coor[:, bv], t)

    return [
        make_node_set(xn_yn_zn, "xn_yn_zn"),  # Vertices node sets
        make_node_set(xp_yn_zn, "xp_yn_zn"),
        make_node_set(xp_yp_zn, "xp_yp_zn"),
        make_node_set(xn_yp_zn, "xn_yp_zn"),
        make_node_set(xn_yn_zp, "xn_yn_zp"),
        make_node_set(xp_yn_zp, "xp_yn_zp"),
        make_node_set(xp_yp_zp, "xp_yp_zp"),
        make_node_set(xn_yp_zp, "xn_yp_zp"),

        make_node_set(xn_yn, "xn_yn"),  # Edges node sets
        make_node_set(xp_yn, "xp_yn"),
        make_node_set(xp_yp, "xp_yp"),
        make_node_set(xn_yp, "xn_yp"),
        make_node_set(xn_zn, "xn_zn"),
        make_node_set(xp_zn, "xp_zn"),
        make_node_set(xp_zp, "xp_zp"),
        make_node_set(xn_zp, "xn_zp"),
        make_node_set(yn_zn, "yn_zn"),
        make_node_set(yp_zn, "yp_zn"),
        make_node_set(yp_zp, "yp_zp"),
        make_node_set(yn_zp, "yn_zp"),
        # 
        make_node_set(xn_face_mask, "xn"),  # Faces node sets
        make_node_set(yn_face_mask, "yn"),
        make_node_set(zn_face_mask, "zn"),
        make_node_set(xp_face_mask, "xp"),
        make_node_set(yp_face_mask, "yp"),
        make_node_set(zp_face_mask, "zp"), 
    ]
end

function parent_child_map(lx::Float64, ly::Float64, lz::Float64)::Dict{String, Vector{String}}
    pc_map::Dict{String, Vector{String}} = Dict(
        "xn" => ["xp",],
        "yn" => ["yp",],
        "zn" => ["zp",],
        "xn_yn" => ["xn_yp", "xp_yn", "xp_yp",],
        "yn_zn" => ["yn_zp", "yp_zn", "yp_zp",],
        "zn_xn" => ["zn_xp", "zp_xn", "zp_xp",],
        "xn_yn_zn" => ["xn_yn_zp", "xn_yp_zn", "xn_yp_zp", "xp_yn_zn", "xp_yn_zp", "xp_yp_zn", "xp_yp_zp",],
    )
    if lx > 0.0 && ly > 0.0 && lz > 0.0  # XYZ
        return pc_map
    elseif lx > 0.0 && ly > 0.0 && lz == 0.0  # XY
        return Dict(k => v for (k, v) in pc_map if !(occursin("z", k)))
    elseif lx == 0.0 && ly > 0.0 && lz > 0.0  # YZ
        return Dict(k => v for (k, v) in pc_map if !(occursin("x", k)))
    elseif lx > 0.0 && ly == 0.0 && lz > 0.0  # ZX
        return Dict(k => v for (k, v) in pc_map if !(occursin("y", k)))
    else
        error("The bounding box is not valid. lx = $lx, ly = $ly, lz = $lz")
    end
end


function get_pc_pairs(
    bns::Vector{NodeSet},
    lx::Float64, ly::Float64, lz::Float64,
)::Dict{NodeSet, Vector{NodeSet}}
    pc_map::Dict{String, Vector{String}} = parent_child_map(lx, ly, lz)
    pc_pairs::Dict{NodeSet, Vector{NodeSet}} = Dict(i=> NodeSet[] for i in bns if haskey(pc_map, i.tag))
    for a_pns in keys(pc_pairs)
        for a_child_tag in pc_map[a_pns.tag]  # a_pns is a parent node set
            for a_ns in bns
                if a_ns.tag == a_child_tag  # a_ns is a child node set
                    push!(pc_pairs[a_pns], a_ns)
                end
            end
        end
    end
    return pc_pairs
end


function assert_node_set_sizes(pns::NodeSet, cns::NodeSet)
    @assert length(pns.tags) > 0 "Empty set of parent node tags; $(pns.tags)"
    @assert length(cns.tags) > 0 "Empty set of child node tags; $(cns.tags)"
    # 
    @assert length(pns.tags) == length(cns.tags) """
    Number of parent node tags and child node tags must be same;
    $(length(pns.tags)) tags on parent $(pns.tag) ≠ $(length(cns.tags)) tags on child $(cns.tag)
    """
    @assert size(pns.coor) == size(cns.coor) "Parent and child node sets must have same shape of coordinates"
    @assert size(pns.tags, 1) == size(pns.coor, 2) "Number of parent tags and parent coordinates must be same,
    Number of parent tags = $(length(pns.tags)) and number of parent coordinates = $(size(pns.coor, 2))"
    @assert size(cns.tags, 1) == size(cns.coor, 2) "Number of child tags and child coordinates must be same,
    Number of child tags = $(length(cns.tags)) and number of child coordinates = $(size(cns.coor, 2))"
    
end


function get_node_pairs(
    parent_node_set::NodeSet,
    child_node_set::NodeSet,
    lx::Float64,
    ly::Float64,
    lz::Float64;
    ϵ::Float64=1.0e-6,
    verbose::Int = 1,
)::Vector{NodePair}
    # d_vec = [0.0, 0.0, 0.0]
    x_shift = occursin("xp", child_node_set.tag) ? lx : 0.0
    y_shift = occursin("yp", child_node_set.tag) ? ly : 0.0
    z_shift = occursin("zp", child_node_set.tag) ? lz : 0.0
    # 
    assert_node_set_sizes(parent_node_set, child_node_set)
    num_nodes::Int = length(parent_node_set.tags)
    pns_tags = parent_node_set.tags
    cns_tags = child_node_set.tags
    pns_coor = parent_node_set.coor
    cns_coor = child_node_set.coor
    node_pairs::Vector{NodePair} = NodePair[]
    # 
    # Sorting nodes on the parent and child node sets
    parent_sorted_data = sortslices(vcat(pns_coor, pns_tags'); dims=2, by=x -> (x[1], x[2], x[3]), rev=false)
    child_sorted_data = sortslices(vcat(cns_coor, cns_tags'); dims=2, by=x -> (x[1], x[2], x[3]), rev=false)
    # It is expected that shifting the parent node by x_shift, y_shift and z_shift
    # should take it to the child node, with tolerance sphere of radius ϵ.
    shifted_parents_coor::Matrix{Float64} = parent_sorted_data[1:3, :] .+ [x_shift; y_shift; z_shift]
    proximities = mapreduce(x -> x^2, +, shifted_parents_coor - child_sorted_data[1:3, :]; dims=1) .< ϵ^2
    if all(proximities)
        if verbose > 2
            @info "Parent-child nodes are matched for all $(num_nodes) nodes 
            of parent $(parent_node_set.tag) and child $(child_node_set.tag)"
        end
        for i = 1:num_nodes
            ap_nx, ap_ny, ap_nz, ap_nt = parent_sorted_data[:, i]
            ac_nx, ac_ny, ac_nz, ac_nt = child_sorted_data[:, i]
            push!(node_pairs,
            NodePair(
                ap_nt,
                ac_nt,
                [ap_nx, ap_ny, ap_nz,],
                [ac_nx, ac_ny, ac_nz,],
                parent_node_set.tag * "<=>" * child_node_set.tag
            )
            )
        end
    else
        println("Parent-child node matching failed for
        $(sum(proximities)) nodes of parent $(parent_node_set.tag) and child $(child_node_set.tag)")
        mask = findall(proximities)
        failed_p_nodes = parent_sorted_data[:, mask]
        failed_c_nodes = child_sorted_data[:, mask]
        println("Failed nodes information:")
        for idx in sum(proximities)
            println("Parent Tag:$(failed_p_nodes[4, idx]), Child Tag:$(failed_c_nodes[4, idx])")
            println("Parent Coor:$(failed_p_nodes[1:3, idx]), Child Coor:$(failed_c_nodes[1:3, idx])")
        end
    end
    return node_pairs
end


function get_node_pairs(
    pc_bns_pairs::Dict{NodeSet, Vector{NodeSet}}, lx::Float64, ly::Float64, lz::Float64;
    small_par::Float64=1.0e-6,
    verbose::Int = 1,
)::Vector{NodePair}
    node_pairs::Vector{NodePair} = NodePair[]
    for (parent_bns, children_bns) in pc_bns_pairs
        for child_bns in children_bns
            append!(node_pairs, get_node_pairs(parent_bns, child_bns, lx, ly, lz; ϵ=small_par, verbose=verbose))
        end
    end
    return node_pairs
end


function make_node_pairs(
    node_set::NodeSet;
    small_par::Float64=1e-06,
    verbose::Int=1
)::Vector{NodePair}
    x_min, y_min, z_min, x_max, y_max, z_max = bbox = bounding_box(node_set.coor)
    lx, ly, lz = x_max - x_min, y_max - y_min, z_max - z_min
    boundary_node_sets::Vector{NodeSet} = get_boundary_node_sets(node_set, bbox; small_par)
    verbose > 2 ? print_node_sets_summary(boundary_node_sets) : nothing
    pc_bns_pairs::Dict{NodeSet, Vector{NodeSet}} = get_pc_pairs(boundary_node_sets, lx, ly, lz)
    return get_node_pairs(pc_bns_pairs, lx, ly, lz; small_par=small_par, verbose=verbose)
end
