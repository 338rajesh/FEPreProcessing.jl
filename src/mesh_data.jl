
function get_rve_dim_tag(bbox::NTuple{6,Float64})
    x_min, y_min, z_min, x_max, y_max, z_max = bbox
    if (x_min != x_max) && (y_min != y_max) && (z_min != z_max)
        return "XYZ"
    elseif (x_min == x_max && y_min != y_max && z_min != z_max)
        return "YZ"
    elseif (x_min != x_max && y_min == y_max && z_min != z_max)
        return "ZX"
    elseif (x_min != x_max && y_min != y_max && z_min == z_max)
        return "XY"
    end
end


"""
    get_rve_node_groups() -> Dict{String,Matrix{IF64}}

Check periodicity and get the groups of nodes on faces, edges and vertices.

Tasks:
    1. Get the groups of nodes on faces, edges and vertices
    2. Verify, if the number of nodes on opposite faces, edges match
    3. Verify, if the vertices groups has single node
    4. Verify, if a node on one face is pairing with more than one node on opposite faces
    5. 

Returns:
--------

    3D Node groups as a dictionary with the following keys\n
    `outer_nodes`,\n
    faces: `xp`, `yp`, `zp`, `xn`, `yn`, `zn` \n
    edges: `zpxp`, `zpxn`, `znxp`, `znxn`, `xpyp`, `xpyn`, `xnyp`, `xnyn`, `ypzp`, `ypzn`, `ynzp`, `ynzn`,\n
    vertices: `xpypzp`, `xpypzn`, `xpynzp`, `xnypzp`, `xnynzp`, `xnypzn`, `xnypzp`, `xnynzn`.

    2D Node groups for zx-plane as a dictionary with the following keys \n
    `outer_nodes`,\n
    edges: `xp`, `xn`, `yp`, `yn`, \n
    vertices: `xpyp`, `xpyn`, `xnyp`, `xnyn`.

"""
function get_rve_node_groups(
    node_tags::Vector{Int},
    node_coordinates::Matrix{Float64},
    rve_bbox::NTuple{6,Float64};
    small_par::Float64=1e-06
)::Dict{String,Matrix{IF64}}
    x_min, y_min, z_min, x_max, y_max, z_max = rve_bbox
    xlb_node::Bool, ylb_node::Bool, zlb_node::Bool = false, false, false
    xub_node::Bool, yub_node::Bool, zub_node::Bool = false, false, false
    #
    a_node_coordinates::NTuple{3,Float64} = (0.0, 0.0, 0.0)
    xn_ntags::Vector{Int} = Int[]
    xn_ncoor::Vector{Float64} = Float64[]
    yn_ntags::Vector{Int} = Int[]
    yn_ncoor::Vector{Float64} = Float64[]
    zn_ntags::Vector{Int} = Int[]
    zn_ncoor::Vector{Float64} = Float64[]
    xp_ntags::Vector{Int} = Int[]
    xp_ncoor::Vector{Float64} = Float64[]
    yp_ntags::Vector{Int} = Int[]
    yp_ncoor::Vector{Float64} = Float64[]
    zp_ntags::Vector{Int} = Int[]
    zp_ncoor::Vector{Float64} = Float64[]
    #
    # edges about x-axis
    ynzn_ntags::Vector{Int} = Int[]
    ynzn_ncoor::Vector{Float64} = Float64[]
    ypzn_ntags::Vector{Int} = Int[]
    ypzn_ncoor::Vector{Float64} = Float64[]
    ypzp_ntags::Vector{Int} = Int[]
    ypzp_ncoor::Vector{Float64} = Float64[]
    ynzp_ntags::Vector{Int} = Int[]
    ynzp_ncoor::Vector{Float64} = Float64[]
    # edges about y-axis
    znxn_ntags::Vector{Int} = Int[]
    znxn_ncoor::Vector{Float64} = Float64[]
    zpxn_ntags::Vector{Int} = Int[]
    zpxn_ncoor::Vector{Float64} = Float64[]
    zpxp_ntags::Vector{Int} = Int[]
    zpxp_ncoor::Vector{Float64} = Float64[]
    znxp_ntags::Vector{Int} = Int[]
    znxp_ncoor::Vector{Float64} = Float64[]
    # edges about z-axis
    xnyn_ntags::Vector{Int} = Int[]
    xnyn_ncoor::Vector{Float64} = Float64[]
    xpyn_ntags::Vector{Int} = Int[]
    xpyn_ncoor::Vector{Float64} = Float64[]
    xpyp_ntags::Vector{Int} = Int[]
    xpyp_ncoor::Vector{Float64} = Float64[]
    xnyp_ntags::Vector{Int} = Int[]
    xnyp_ncoor::Vector{Float64} = Float64[]
    #
    xpypzp_ntags::Vector{Int} = Int[]
    xpypzp_ncoor::Vector{Float64} = Float64[]
    xpypzn_ntags::Vector{Int} = Int[]
    xpypzn_ncoor::Vector{Float64} = Float64[]
    xpynzn_ntags::Vector{Int} = Int[]
    xpynzn_ncoor::Vector{Float64} = Float64[]
    xpynzp_ntags::Vector{Int} = Int[]
    xpynzp_ncoor::Vector{Float64} = Float64[]
    xnypzp_ntags::Vector{Int} = Int[]
    xnypzp_ncoor::Vector{Float64} = Float64[]
    xnypzn_ntags::Vector{Int} = Int[]
    xnypzn_ncoor::Vector{Float64} = Float64[]
    xnynzn_ntags::Vector{Int} = Int[]
    xnynzn_ncoor::Vector{Float64} = Float64[]
    xnynzp_ntags::Vector{Int} = Int[]
    xnynzp_ncoor::Vector{Float64} = Float64[]
    #
    outer_ntags::Vector{Int} = Int[]
    outer_ncoor::Vector{Float64} = Float64[]
    #
    # Making the groups of nodes
    for (k, anode_tag) in enumerate(node_tags)
        nx, ny, nz = node_coordinates[:, k]
        xlb_node = x_min - small_par < nx < x_min + small_par
        ylb_node = y_min - small_par < ny < y_min + small_par
        zlb_node = z_min - small_par < nz < z_min + small_par
        xub_node = x_max - small_par < nx < x_max + small_par
        yub_node = y_max - small_par < ny < y_max + small_par
        zub_node = z_max - small_par < nz < z_max + small_par
        # a_node_coordinates = (nx, ny, nz)
        if (xlb_node || xub_node || ylb_node || yub_node || zlb_node || zub_node)
            #
            push!(outer_ntags, anode_tag)
            push!(outer_ncoor, nx, ny, nz)
            # only faces
            if (xlb_node && !ylb_node && !zlb_node && !yub_node && !zub_node)
                push!(xn_ntags, anode_tag)
                push!(xn_ncoor, nx, ny, nz)  # exclusively xn
            elseif (xub_node && !ylb_node && !zlb_node && !yub_node && !zub_node)
                push!(xp_ntags, anode_tag)
                push!(xp_ncoor, nx, ny, nz)  # exclusively xp
            elseif (ylb_node && !xlb_node && !zlb_node && !xub_node && !zub_node)
                push!(yn_ntags, anode_tag)
                push!(yn_ncoor, nx, ny, nz)  # exclusively yn
            elseif (yub_node && !xlb_node && !zlb_node && !xub_node && !zub_node)
                push!(yp_ntags, anode_tag)
                push!(yp_ncoor, nx, ny, nz)  # exclusively yp
            elseif (zlb_node && !xlb_node && !ylb_node && !xub_node && !yub_node)
                push!(zn_ntags, anode_tag)
                push!(zn_ncoor, nx, ny, nz)  # exclusively zn
            elseif (zub_node && !xlb_node && !ylb_node && !xub_node && !yub_node)
                push!(zp_ntags, anode_tag)
                push!(zp_ncoor, nx, ny, nz)  # exclusively zp
            #
            # only edges
            elseif (ylb_node && zlb_node && !xlb_node && !xub_node)
                push!(ynzn_ntags, anode_tag)
                push!(ynzn_ncoor, nx, ny, nz)  # exclusively ynzn
            elseif (yub_node && zlb_node && !xlb_node && !xub_node)
                push!(ypzn_ntags, anode_tag)
                push!(ypzn_ncoor, nx, ny, nz)  # exclusively ypzn
            elseif (yub_node && zub_node && !xlb_node && !xub_node)
                push!(ypzp_ntags, anode_tag)
                push!(ypzp_ncoor, nx, ny, nz)  # exclusively ypzp
            elseif (zub_node && ylb_node && !xlb_node && !xub_node)
                push!(ynzp_ntags, anode_tag)
                push!(ynzp_ncoor, nx, ny, nz)  # exclusively ynzp
            #
            elseif (xlb_node && zlb_node && !ylb_node && !yub_node)
                push!(znxn_ntags, anode_tag)
                push!(znxn_ncoor, nx, ny, nz)  # exclusively znxn
            elseif (xlb_node && zub_node && !ylb_node && !yub_node)
                push!(zpxn_ntags, anode_tag)
                push!(zpxn_ncoor, nx, ny, nz)  # exclusively zpxn
            elseif (xub_node && zub_node && !ylb_node && !yub_node)
                push!(zpxp_ntags, anode_tag)
                push!(zpxp_ncoor, nx, ny, nz)  # exclusively zpxp
            elseif (xub_node && zlb_node && !ylb_node && !yub_node)
                push!(znxp_ntags, anode_tag)
                push!(znxp_ncoor, nx, ny, nz)  # exclusively znxp
            #
            elseif (xlb_node && ylb_node && !zlb_node && !zub_node)
                push!(xnyn_ntags, anode_tag)
                push!(xnyn_ncoor, nx, ny, nz)  # exclusively xnyn
            elseif (xub_node && ylb_node && !zlb_node && !zub_node)
                push!(xpyn_ntags, anode_tag)
                push!(xpyn_ncoor, nx, ny, nz)  # exclusively xpyn
            elseif (xub_node && yub_node && !zlb_node && !zub_node)
                push!(xpyp_ntags, anode_tag)
                push!(xpyp_ncoor, nx, ny, nz)  # exclusively xpyp
            elseif (xlb_node && yub_node && !zlb_node && !zub_node)
                push!(xnyp_ntags, anode_tag)
                push!(xnyp_ncoor, nx, ny, nz)  # exclusively xnyp
            #
            # only vertices
            elseif (xub_node && yub_node && zub_node)
                push!(xpypzp_ntags, anode_tag)
                push!(xpypzp_ncoor, nx, ny, nz)  # exclusively xpypzp
            elseif (xub_node && yub_node && zlb_node)
                push!(xpypzn_ntags, anode_tag)
                push!(xpypzn_ncoor, nx, ny, nz)  # exclusively xpypzn
            elseif (xub_node && ylb_node && zlb_node)
                push!(xpynzn_ntags, anode_tag)
                push!(xpynzn_ncoor, nx, ny, nz)  # exclusively xpynzn
            elseif (xub_node && ylb_node && zub_node)
                push!(xpynzp_ntags, anode_tag)
                push!(xpynzp_ncoor, nx, ny, nz)  # exclusively xpynzp
            elseif (xlb_node && yub_node && zub_node)
                push!(xnypzp_ntags, anode_tag)
                push!(xnypzp_ncoor, nx, ny, nz)  # exclusively xnypzp
            elseif (xlb_node && yub_node && zlb_node)
                push!(xnypzn_ntags, anode_tag)
                push!(xnypzn_ncoor, nx, ny, nz)  # exclusively xnypzn
            elseif (xlb_node && ylb_node && zlb_node)
                push!(xnynzn_ntags, anode_tag)
                push!(xnynzn_ncoor, nx, ny, nz)  # exclusively xnynzn
            elseif (xlb_node && ylb_node && zub_node)
                push!(xnynzp_ntags, anode_tag)
                push!(xnynzp_ncoor, nx, ny, nz)  # exclusively xnynzp
            end
        end
    end
    #
    return Dict(
        "xp" => IF64[reshape(xp_ntags, 1, :); reshape(xp_ncoor, 3, :)],
        "xn" => IF64[reshape(xn_ntags, 1, :); reshape(xn_ncoor, 3, :)],
        "yp" => IF64[reshape(yp_ntags, 1, :); reshape(yp_ncoor, 3, :)],
        "yn" => IF64[reshape(yn_ntags, 1, :); reshape(yn_ncoor, 3, :)],
        "zp" => IF64[reshape(zp_ntags, 1, :); reshape(zp_ncoor, 3, :)],
        "zn" => IF64[reshape(zn_ntags, 1, :); reshape(zn_ncoor, 3, :)],
        "ypzn" => IF64[reshape(ypzn_ntags, 1, :); reshape(ypzn_ncoor, 3, :)],
        "ynzp" => IF64[reshape(ynzp_ntags, 1, :); reshape(ynzp_ncoor, 3, :)],
        "ypzp" => IF64[reshape(ypzp_ntags, 1, :); reshape(ypzp_ncoor, 3, :)],
        "ynzn" => IF64[reshape(ynzn_ntags, 1, :); reshape(ynzn_ncoor, 3, :)],
        "znxp" => IF64[reshape(znxp_ntags, 1, :); reshape(znxp_ncoor, 3, :)],
        "zpxn" => IF64[reshape(zpxn_ntags, 1, :); reshape(zpxn_ncoor, 3, :)],
        "zpxp" => IF64[reshape(zpxp_ntags, 1, :); reshape(zpxp_ncoor, 3, :)],
        "znxn" => IF64[reshape(znxn_ntags, 1, :); reshape(znxn_ncoor, 3, :)],
        "xnyp" => IF64[reshape(xnyp_ntags, 1, :); reshape(xnyp_ncoor, 3, :)],
        "xpyn" => IF64[reshape(xpyn_ntags, 1, :); reshape(xpyn_ncoor, 3, :)],
        "xpyp" => IF64[reshape(xpyp_ntags, 1, :); reshape(xpyp_ncoor, 3, :)],
        "xnyn" => IF64[reshape(xnyn_ntags, 1, :); reshape(xnyn_ncoor, 3, :)],
        "xpynzn" => IF64[reshape(xpynzn_ntags, 1, :); reshape(xpynzn_ncoor, 3, :)],
        "xnypzn" => IF64[reshape(xnypzn_ntags, 1, :); reshape(xnypzn_ncoor, 3, :)],
        "xnynzp" => IF64[reshape(xnynzp_ntags, 1, :); reshape(xnynzp_ncoor, 3, :)],
        "xnypzp" => IF64[reshape(xnypzp_ntags, 1, :); reshape(xnypzp_ncoor, 3, :)],
        "xpynzp" => IF64[reshape(xpynzp_ntags, 1, :); reshape(xpynzp_ncoor, 3, :)],
        "xpypzn" => IF64[reshape(xpypzn_ntags, 1, :); reshape(xpypzn_ncoor, 3, :)],
        "xpypzp" => IF64[reshape(xpypzp_ntags, 1, :); reshape(xpypzp_ncoor, 3, :)],
        "xnynzn" => IF64[reshape(xnynzn_ntags, 1, :); reshape(xnynzn_ncoor, 3, :)],
        "outer_nodes" => IF64[reshape(outer_ntags, 1, :); reshape(outer_ncoor, 3, :)],
    )
end

#
function get_rve_node_groups(
    node_tags::Vector{Int},
    node_coordinates::Matrix{Float64},
    rve_bbox::NTuple{4,Float64};
    small_par::Float64=1e-06
)::Dict{String,Matrix{IF64}}
    x_min, y_min, x_max, y_max = rve_bbox
    xlb_node::Bool, ylb_node::Bool = false, false
    xub_node::Bool, yub_node::Bool = false, false
    #
    a_node_coordinates::NTuple{3,Float64} = (0.0, 0.0, 0.0)
    xn_ntags::Vector{Int} = Int[]
    xn_ncoor::Vector{Float64} = Float64[]
    yn_ntags::Vector{Int} = Int[]
    yn_ncoor::Vector{Float64} = Float64[]
    xp_ntags::Vector{Int} = Int[]
    xp_ncoor::Vector{Float64} = Float64[]
    yp_ntags::Vector{Int} = Int[]
    yp_ncoor::Vector{Float64} = Float64[]
    #
    # Vertices
    xnyn_ntags::Vector{Int} = Int[]
    xnyn_ncoor::Vector{Float64} = Float64[]
    xpyn_ntags::Vector{Int} = Int[]
    xpyn_ncoor::Vector{Float64} = Float64[]
    xpyp_ntags::Vector{Int} = Int[]
    xpyp_ncoor::Vector{Float64} = Float64[]
    xnyp_ntags::Vector{Int} = Int[]
    xnyp_ncoor::Vector{Float64} = Float64[]
    #
    outer_ntags::Vector{Int} = Int[]
    outer_ncoor::Vector{Float64} = Float64[]
    #
    # Making the groups of nodes
    for (k, anode_tag) in enumerate(node_tags)
        nx, ny = node_coordinates[:, k]
        xlb_node = x_min - small_par < nx < x_min + small_par
        ylb_node = y_min - small_par < ny < y_min + small_par
        xub_node = x_max - small_par < nx < x_max + small_par
        yub_node = y_max - small_par < ny < y_max + small_par
        # a_node_coordinates = (nx, ny, nz)
        if (xlb_node || xub_node || ylb_node || yub_node)
            #
            push!(outer_ntags, anode_tag)
            push!(outer_ncoor, nx, ny,)
            # only edges
            if (xlb_node && !ylb_node && !yub_node)   # exclusively xn
                push!(xn_ntags, anode_tag)
                push!(xn_ncoor, nx, ny,)
            elseif (xub_node && !ylb_node && !yub_node)  # exclusively xp
                push!(xp_ntags, anode_tag)
                push!(xp_ncoor, nx, ny,)
            elseif (ylb_node && !xlb_node && !xub_node)  # exclusively yn
                push!(yn_ntags, anode_tag)
                push!(yn_ncoor, nx, ny,)
            elseif (yub_node && !xlb_node && !xub_node)  # exclusively yp
                push!(yp_ntags, anode_tag)
                push!(yp_ncoor, nx, ny,)
            #
            # only vertices
            elseif (xlb_node && ylb_node)  # exclusively xnyn
                push!(xnyn_ntags, anode_tag)
                push!(xnyn_ncoor, nx, ny,)
            elseif (xub_node && ylb_node)  # exclusively xpyn
                push!(xpyn_ntags, anode_tag)
                push!(xpyn_ncoor, nx, ny,)
            elseif (xub_node && yub_node)  # exclusively xpyp
                push!(xpyp_ntags, anode_tag)
                push!(xpyp_ncoor, nx, ny,)
            elseif (xlb_node && yub_node)  # exclusively xnyp
                push!(xnyp_ntags, anode_tag)
                push!(xnyp_ncoor, nx, ny,)
            end
        end
    end
    #
    return Dict(
        "xp" => IF64[reshape(xp_ntags, 1, :); reshape(xp_ncoor, 2, :)],
        "xn" => IF64[reshape(xn_ntags, 1, :); reshape(xn_ncoor, 2, :)],
        "yp" => IF64[reshape(yp_ntags, 1, :); reshape(yp_ncoor, 2, :)],
        "yn" => IF64[reshape(yn_ntags, 1, :); reshape(yn_ncoor, 2, :)],
        "xnyp" => IF64[reshape(xnyp_ntags, 1, :); reshape(xnyp_ncoor, 2, :)],
        "xpyn" => IF64[reshape(xpyn_ntags, 1, :); reshape(xpyn_ncoor, 2, :)],
        "xpyp" => IF64[reshape(xpyp_ntags, 1, :); reshape(xpyp_ncoor, 2, :)],
        "xnyn" => IF64[reshape(xnyn_ntags, 1, :); reshape(xnyn_ncoor, 2, :)],
        "outer_nodes" => IF64[reshape(outer_ntags, 1, :); reshape(outer_ncoor, 2, :)],
    )
end


function _check_periodicity_error(
    g1::Tuple{SubString{String}, Matrix{IF64}},
    g2::Tuple{SubString{String}, Matrix{IF64}},
)
    g1_id, g2_id = g1[1], g2[1]
    g1_ninfo, g2_ninfo = g1[2], g2[2]
    num_g1 = size(g1_ninfo, 2)
    num_g2 = size(g2_ninfo, 2)

    if num_g1 != num_g2
        msg = "Unequal number of nodes on opposite groups $(num_g1)@$(g1_id) != $(num_g2)@$(g2_id)."
        high_node_group_id = num_g1 > num_g2 ? g1_id : g2_id
        high_node_group = num_g1 > num_g2 ? g1_ninfo : g2_ninfo
        msg *= "\n$high_node_group_id contains $(abs(num_g1-num_g2)) extra nodes."
        throw(FEPrepError(msg))
    end
end


""" Creates Vector of Node Pairs for node groups in 3D """
function node_pairs(
    png::Matrix{IF64},  # png: parent node group
    cng::Matrix{IF64},  # cng: child node group
    shift::NTuple{3,Float64},
    pair_tag::String;
    ϵ::Float64=1e-06
)::Vector{FENodePair}
    δx, δy, δz = shift
    node_pairs::Vector{FENodePair} = []
    #
    parent_id, child_id = split(pair_tag, "_")
    # rounding to a place higher than tolerance order. This is to avoid the influence of machine
    # precision on the sorting process.
    png[2:end, :] = round.(png[2:end, :]; digits=(1 + ceil(Int, log10(1 / ϵ))))  # rounding to a place higher than tolerance order
    cng[2:end, :] = round.(cng[2:end, :]; digits=(1 + ceil(Int, log10(1 / ϵ))))
    #
    sorted_png = sortslices(png, dims=2, by=x -> (x[2], x[3], x[4]), rev=false)
    sorted_cng = sortslices(cng, dims=2, by=x -> (x[2], x[3], x[4]), rev=false)
    _check_periodicity_error((parent_id, sorted_png), (child_id, sorted_cng))
    #
    num_nodes = size(png, 2)
    for i = 1:num_nodes
        apg_nt, apg_nx, apg_ny, apg_nz = sorted_png[:, i]
        acg_nt, acg_nx, acg_ny, acg_nz = sorted_cng[:, i]  # TODO check if vectorization increases performance
        if (
            abs(abs(apg_nx - acg_nx) - δx) < ϵ &&
            abs(abs(apg_ny - acg_ny) - δy) < ϵ &&
            abs(abs(apg_nz - acg_nz) - δz) < ϵ
        )
            push!(node_pairs, FENodePair("$parent_id-$apg_nt-$child_id-$acg_nt", apg_nt, acg_nt))
        else
            throw(MeshPeriodicityError(apg_nt, "$parent_id-$apg_nt-$child_id-$acg_nt"))
        end
    end
    return node_pairs
end

""" Created Vector of Node Pairs for node groups in 2D """
function node_pairs(
    png::Matrix{IF64},  # png: parent node group
    cng::Matrix{IF64},  # cng: child node group
    shift::NTuple{2,Float64},
    pair_tag::String;
    ϵ::Float64=1e-06
)::Vector{FENodePair}
    δx, δy = shift
    node_pairs::Vector{FENodePair} = []
    #
    parent_id, child_id = split(pair_tag, "_")
    # rounding to a place higher than tolerance order. This is to avoid the influence of machine
    # precision on the sorting process.
    png[2:end, :] = round.(png[2:end, :]; digits=(1 + ceil(Int, log10(1 / ϵ))))  # rounding to a place higher than tolerance order
    cng[2:end, :] = round.(cng[2:end, :]; digits=(1 + ceil(Int, log10(1 / ϵ))))
    #
    sorted_png = sortslices(png, dims=2, by=x -> (x[2], x[3]), rev=false)
    sorted_cng = sortslices(cng, dims=2, by=x -> (x[2], x[3]), rev=false)
    if size(sorted_png, 2) != size(sorted_cng, 2)
        throw(FEPrepError("Unequal number of nodes on opposite groups"))
    end
    num_nodes = size(png, 2)
    for i = 1:num_nodes
        apg_nt, apg_nx, apg_ny = sorted_png[:, i]
        acg_nt, acg_nx, acg_ny = sorted_cng[:, i]  # TODO check if vectorization increases performance
        if (
            abs(abs(apg_nx - acg_nx) - δx) < ϵ &&
            abs(abs(apg_ny - acg_ny) - δy) < ϵ
        )
            push!(node_pairs, FENodePair("$parent_id-$apg_nt-$child_id-$acg_nt", apg_nt, acg_nt))
        else
            throw(MeshPeriodicityError(ang1_nt, pair_tag))
        end
    end
    return node_pairs
end


function node_pairs(
    node_groups::Dict{String,Matrix{IF64}},
    bbox::NTuple{6,Float64},
)::Vector{FENodePair}
    x_min, y_min, z_min, x_max, y_max, z_max = bbox
    lx, ly, lz = x_max - x_min, y_max - y_min, z_max - z_min
    #
    parent_child_map = Dict(
            "xn" => [("xp", (lx, 0.0, 0.0))],
            "yn" => [("yp", (0.0, ly, 0.0))],
            "zn" => [("zp", (0.0, 0.0, lz))],
            "xnyn" => [
                ("xnyp", (0.0, ly, 0.0)),
                ("xpyn", (lx, 0.0, 0.0)),
                ("xpyp", (lx, ly, 0.0)),
            ],
            "ynzn" => [
                ("ynzp", (0.0, 0.0, lz)),
                ("ypzn", (0.0, ly, 0.0)),
                ("ypzp", (0.0, ly, lz)),
            ],
            "znxn" => [
                ("znxp", (lx, 0.0, 0.0)),
                ("zpxn", (0.0, 0.0, lz)),
                ("zpxp", (lx, 0.0, lz)),
            ],
            "xnynzn" => [
                ("xnynzp", (0.0, 0.0, lz)),
                ("xnypzn", (0.0, ly, 0.0)),
                ("xnypzp", (0.0, ly, lz)),
                ("xpynzn", (lx, 0.0, 0.0)),
                ("xpynzp", (lx, 0.0, lz)),
                ("xpypzn", (lx, ly, 0.0)),
                ("xpypzp", (lx, ly, lz)),
            ],
        )
    fe_node_pairs::Vector{FENodePair} = FENodePair[]
    for (a_parent, children) in parent_child_map
        for (a_child, distance) in children
            if (a_parent in keys(node_groups) &&
                a_child in keys(node_groups))
                append!(
                    fe_node_pairs,
                    node_pairs(
                        node_groups[a_parent],
                        node_groups[a_child],
                        distance,
                        a_parent * "_" * a_child,
                    ),
                )
            end
        end
    end
    return fe_node_pairs
end


function node_pairs(
    node_groups::Dict{String,Matrix{IF64}},
    bbox::NTuple{4,Float64},
)::Vector{FENodePair}
    x_min, y_min, x_max, y_max= bbox
    lx, ly= x_max - x_min, y_max - y_min
    #
    parent_child_map = Dict(
            "xn" => [("xp", (lx, 0.0))],
            "yn" => [("yp", (0.0, ly))],
            "xnyn" => [
                ("xnyp", (0.0, ly)),
                ("xpyn", (lx, 0.0)),
                ("xpyp", (lx, ly)),
            ],
        )
    fe_node_pairs::Vector{FENodePair} = FENodePair[]
    for (a_parent, children) in parent_child_map
        for (a_child, distance) in children
            if (a_parent in keys(node_groups) &&
                a_child in keys(node_groups))
                append!(
                    fe_node_pairs,
                    node_pairs(
                        node_groups[a_parent],
                        node_groups[a_child],
                        distance,
                        a_parent * "_" * a_child,
                    ),
                )
            end
        end
    end
    return fe_node_pairs
end


"""
    make_rve_node_pairs(ntags, ncoor)

 Returns a Vector of finite element node pairs (FENodePair).

"""
function make_rve_node_pairs(
    node_tags::Vector{Int},
    node_coordinates::Matrix{Float64};
    small_para::Float64=1e-06
)::Vector{FENodePair}
    if length(node_tags) != size(node_coordinates, 2)
        throw(
            DimensionMismatch(
                "Number of node tags and number of columns in coordinate materix are not equal. Please ensure that, size of the coordinate matrix is (3, num_tags).
                Ensure that each column of the coordinate matrix represents a node coordinates.
                This is to take the advantage of column major representation of Julia.",
            ),
        )
    end
    ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = extrema(node_coordinates, dims=2)
    rve_bbox::NTuple{6,Float64} = (x_min, y_min, z_min, x_max, y_max, z_max)
    rve_dim_tag = get_rve_dim_tag(rve_bbox)
    if rve_dim_tag == "XYZ"
        rve_node_groups = get_rve_node_groups(node_tags, node_coordinates, rve_bbox; small_par=small_para)
        return node_pairs(rve_node_groups, rve_bbox)
    elseif rve_dim_tag == "XY"
        bbox = (x_min, y_min, x_max, y_max)
        rve_node_groups = get_rve_node_groups(node_tags, node_coordinates, bbox; small_par=small_para)
        return node_pairs(rve_node_groups, bbox)
    end    
end


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



