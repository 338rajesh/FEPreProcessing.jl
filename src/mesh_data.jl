

"""
Check periodicity and get the node pairs to be contrained
Tasks:
    1. Get the groups of nodes on faces, edges and vertices
    2. Verify, if the number of nodes on opposite faces, edges match
    3. Verify, if the vertices groups has single node
    4. Verify, if a node on one face is pairing with morethan one node on opposite faces
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
    small_par::Float64 = 1e-06
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
    get_group_data(v1::Vector{Int}, v2::Vector{Float64}) = begin
        return IF64[reshape(v1, 1, :); reshape(v2, 3, :)]
    end
    #
    if ((x_min != x_max) && (y_min != y_max) && (z_min != z_max))
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
    elseif (x_min == x_max && y_min != y_max && z_min != z_max)  # yz-plane
        return Dict(
            "yp" => IF64[reshape(yp_ntags, 1, :); reshape(yp_ncoor, 3, :)],
            "yn" => IF64[reshape(yn_ntags, 1, :); reshape(yn_ncoor, 3, :)],
            "zp" => IF64[reshape(zp_ntags, 1, :); reshape(zp_ncoor, 3, :)],
            "zn" => IF64[reshape(zn_ntags, 1, :); reshape(zn_ncoor, 3, :)],
            "ypzn" => IF64[reshape(ypzn_ntags, 1, :); reshape(ypzn_ncoor, 3, :)],
            "ynzp" => IF64[reshape(ynzp_ntags, 1, :); reshape(ynzp_ncoor, 3, :)],
            "ypzp" => IF64[reshape(ypzp_ntags, 1, :); reshape(ypzp_ncoor, 3, :)],
            "ynzn" => IF64[reshape(ynzn_ntags, 1, :); reshape(ynzn_ncoor, 3, :)],
            "outer_nodes" => IF64[reshape(outer_ntags, 1, :); reshape(outer_ncoor, 3, :)],
        )
    elseif (x_min != x_max && y_min == y_max && z_min != z_max)  # zx-plane
        return Dict(
            "zp" => [reshape(zp_ntags, 1, :); reshape(zp_ncoor, 3, :)],
            "zn" => [reshape(zn_ntags, 1, :); reshape(zn_ncoor, 3, :)],
            "xp" => [reshape(xp_ntags, 1, :); reshape(xp_ncoor, 3, :)],
            "xn" => [reshape(xn_ntags, 1, :); reshape(xn_ncoor, 3, :)],
            "znxp" => [reshape(znxp_ntags, 1, :); reshape(znxp_ncoor, 3, :)],
            "zpxn" => [reshape(zpxn_ntags, 1, :); reshape(zpxn_ncoor, 3, :)],
            "zpxp" => [reshape(zpxp_ntags, 1, :); reshape(zpxp_ncoor, 3, :)],
            "znxn" => [reshape(znxn_ntags, 1, :); reshape(znxn_ncoor, 3, :)],
            "outer_nodes" => [reshape(outer_ntags, 1, :); reshape(outer_ncoor, 3, :)],
        )
    elseif (x_min != x_max && y_min != y_max && z_min == z_max)  # xy-plane
        return Dict(
            "xp" => [reshape(xp_ntags, 1, :); reshape(xp_ncoor, 3, :)],
            "xn" => [reshape(xn_ntags, 1, :); reshape(xn_ncoor, 3, :)],
            "yp" => [reshape(yp_ntags, 1, :); reshape(yp_ncoor, 3, :)],
            "yn" => [reshape(yn_ntags, 1, :); reshape(yn_ncoor, 3, :)],
            "xnyp" => [reshape(xnyp_ntags, 1, :); reshape(xnyp_ncoor, 3, :)],
            "xpyn" => [reshape(xpyn_ntags, 1, :); reshape(xpyn_ncoor, 3, :)],
            "xpyp" => [reshape(xpyp_ntags, 1, :); reshape(xpyp_ncoor, 3, :)],
            "xnyn" => [reshape(xnyn_ntags, 1, :); reshape(xnyn_ncoor, 3, :)],
            "outer_nodes" => [reshape(outer_ntags, 1, :); reshape(outer_ncoor, 3, :)],
        )
    else
        throw(FEPrepError("Node groups on lines are not being handled...!"))
    end
end



function node_pairs(
    ng1::Matrix{IF64},
    ng2::Matrix{IF64},
    shift::NTuple{3,Float64},
    pair_tag::String;
    ϵ::Float64 = 1e-06
)::Vector{FENodePair}
    δx, δy, δz = shift
    node_pairs::Vector{FENodePair} = []
    #
    # rounding to a place higher than tolerance order. This is to avoid the influence of machine
    # precision on the sorting process.
    ng1[2:end, :] = round.(ng1[2:end, :]; digits=(1 + ceil(Int, log10(1/ϵ))))  # rounding to a place higher than tolerance order
    ng2[2:end, :] = round.(ng2[2:end, :]; digits=(1 + ceil(Int, log10(1/ϵ))))  
    #
    sorted_ng1 = sortslices(ng1, dims = 2, by = x -> (x[2], x[3], x[4]), rev = false)
    sorted_ng2 = sortslices(ng2, dims = 2, by = x -> (x[2], x[3], x[4]), rev = false)
    if size(sorted_ng1, 2) != size(sorted_ng2, 2)
        throw(FEPrepError("Unequal number of nodes on opposite groups"))
    end
    num_nodes = size(ng1, 2)
    for i = 1:num_nodes
        ang1_nt, ang1_nx, ang1_ny, ang1_nz = sorted_ng1[:, i]
        ang2_nt, ang2_nx, ang2_ny, ang2_nz = sorted_ng2[:, i]  # TODO check if vectorization increases performance
        if (
            abs(abs(ang1_nx - ang2_nx) - δx) < ϵ &&
            abs(abs(ang1_ny - ang2_ny) - δy) < ϵ &&
            abs(abs(ang1_nz - ang2_nz) - δz) < ϵ
        )
            push!(node_pairs, FENodePair(pair_tag, ang1_nt, ang2_nt))
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
        "xp" => [("xn", (lx, 0.0, 0.0))],
        "yp" => [("yn", (0.0, ly, 0.0))],
        "zp" => [("zn", (0.0, 0.0, lz))],
        "ypzp" => [
            ("ypzn", (0.0, 0.0, lz)),
            ("ynzp", (0.0, ly, 0.0)),
            ("ynzn", (0.0, ly, lz)),
        ],
        "zpxp" => [
            ("znxn", (lx, 0.0, lz)),
            ("znxp", (0.0, 0.0, lz)),
            ("zpxn", (lx, 0.0, 0.0)),
        ],
        "xpyp" => [
            ("xnyn", (lx, ly, 0.0)),
            ("xpyn", (0.0, ly, 0.0)),
            ("xnyp", (lx, 0.0, 0.0)),
        ],
        "xpypzp" => [
            ("xnynzn", (lx, ly, lz)),
            ("xnynzp", (lx, ly, 0.0)),
            ("xnypzn", (lx, 0.0, lz)),
            ("xpynzn", (0.0, ly, lz)),
            ("xnypzp", (lx, 0.0, 0.0)),
            ("xpynzp", (0.0, ly, 0.0)),
            ("xpypzn", (0.0, 0.0, lz)),
        ],
    )
    fe_node_pairs::Vector{FENodePair} = FENodePair[]
    for a_pid in keys(parent_child_map)
        if a_pid in keys(node_groups)
            for (a_cid, shift) in parent_child_map[a_pid]
                append!(
                    fe_node_pairs,
                    node_pairs(
                        node_groups[a_pid],
                        node_groups[a_cid],
                        shift,
                        a_pid * "_" * a_cid,
                    ),
                )
            end
        end
    end
    return fe_node_pairs
end


"""
    make_rve_node_pairs(ntags, ncoor)

 Returns a tuple of finite element node pairs (FENodePair).

"""
function make_rve_node_pairs(
    node_tags::Vector{Int},
    node_coordinates::Matrix{Float64};
    small_para::Float64 = 1e-06
)::Vector{FENodePair}
    println("Making RVE node pairs...!")
    if length(node_tags) != size(node_coordinates, 2)
        throw(
            DimensionMismatch(
                "Number of node tags and number of columns in coordinate materix are not equal. Please ensure that, size of the coordinate matrix is (3, num_tags).
                Ensure that each column of the coordinate matrix represents a node coordinates.
                This is to take the advantage of column major representation of Julia.",
            ),
        )
    end
    ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = extrema(node_coordinates, dims = 2)
    rve_bbox::NTuple{6,Float64} = (x_min, y_min, z_min, x_max, y_max, z_max)
    rve_node_groups =
        get_rve_node_groups(node_tags, node_coordinates, rve_bbox; small_par = small_para)
    return node_pairs(rve_node_groups, rve_bbox)
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
     Element conne gmsh ids are used to identify element sets which will eventually be converted to ABAQUS convention of element naming.
 - **nodal_data**::Dict{Int64, Tuple{Vararg{Float64}}},    
 - **elset_tag**::Union{Int64, String};
 - **num_ip**::Int64=-1,
 - **gq_order**::Int64=-1,
 - **print_debug_info**::Bool=false,

"""


function make_finite_element_set(
    element_type::Int,
    element_connectivity::Matrix{Int},
    nodal_data::Dict{Int, Vector{Float64}};
    elset_tag::Union{String, Int} = "Element-Set",
    material::Union{Nothing, Material} = nothing,
    num_ip::Int64 = -1,
    gq_order::Int64 = -1,
)::FiniteElementSet
    #
    println("Making element set with ID ", elset_tag)
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
            push!(shfg_xyz, SMatrix{1 + dim, num_shape_func, Float64}(aip_shf_grad_xyz))
        end
        push!(
            elements,
            fele_data_type(
                ele_tag,
                SVector{num_node_per_fele, Int}(ele_nodes),
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
    element_connectivities::Dict{Int64, Matrix{Int64}},
    nodal_data::Dict{Int64, Vector{Float64}},
    ele_material::Material,
    elset_tag::Union{Int64,String};
    num_ip::Int64 = -1,
    gq_order::Int64 = -1,
)::Vector{FiniteElementSet}
    #
    return [make_finite_element_set(
        a_ele_type,
        a_ele_type_ele_connectivity,
        nodal_data;
        elset_tag=elset_tag,
        material=ele_material,
        num_ip=num_ip,
        gq_order=gq_order,
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
    nodal_dof::Tuple{Vararg{String}} = ("x", "y", "z")
)::Dict{Int64, Vector{Float64}}  # TODO make node sets object oriented
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
    print_info::Bool = false,
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
        println("Sum of volumes of ", element_counter," elements is  ", total_volume)    
    end
    return total_volume
end





# ================================================
#                   Archives
# ================================================

# function get_felement_sets_a(
#     element_connectivity::Dict{Int64,Matrix{Int64}},
#     nodal_data::Dict{Int64,Vector{Float64}},
#     elset_tag::Union{Int64,String};
#     num_ip::Int64 = -1,
#     gq_order::Int64 = -1,
#     print_debug_info::Bool = false
# )::Vector{AbstractFElementSet}
#     #
#     println("Preparing $elset_tag...!")
#     element_sets::Vector{AbstractFElementSet} = []
#     elements::Vector{AbstractFElement} = []
#     num_nodes::Int64 = 0
#     num_shf::Int64 = 0
#     num_felements::Int64 = 0
#     ivols::Vector{Float64} = Float64[]
#     shfg_xyz::Vector{SMatrix} = SMatrix[]
#     ele_tag::Int = 0
#     ele_nodes::Vector{Int} = Int[]
#     #
#     # TODO add element type assertion
#     tot_vol::Float64 = 0.0   # *****************
#     #
#     # Running over each element type, with their element_connectivity information.
#     for (ele_type, ele_con) in element_connectivity
#         empty!(elements)
#         num_nodes = size(ele_con)[1] - 1  # as the first row corresponds to element tag
#         num_felements = size(ele_con)[2]
#         num_shf = num_nodes
#         dim, fele_data_type, integration_points = get_felement_properties(
#             ele_type,
#             gq_order,
#             num_ip
#         )
#         #
#         # For the given element type, get shape functions and its gradients at
#         # gauss points in parent coordinate system
#         shf_data_parent_csys = get_shape_functions_info(
#             fele_data_type,
#             integration_points
#         )
#         #
#         if print_debug_info
#             println(
#                 "Element Data Type: $fele_data_type,
#                 Number of elements :: $(num_felements),
#                 Number of nodes/element :: $num_nodes,
#                 Integration points :: $integration_points,",
#             )
#         end
#         #
#         # Running over each element of the given type.
#         for a_col in eachcol(ele_con)
#             empty!(ivols)
#             empty!(shfg_xyz)
#             #
#             ele_tag, ele_nodes = a_col[1], a_col[2:end]
#             nodal_coordinates = (nodal_data[i][1:dim] for i in ele_nodes)  # NTuple{N, Vector{Float64}}
#             #
#             # run over each gauss point and get Jacobian and its determinant
#             for (a_ip, shfg_parent) in shf_data_parent_csys
#                 Jacobian = zeros(Float64, dim, dim)
#                 # run over pairs of coordinates and shape function data
#                 for (k2, a_node_coor) in enumerate(nodal_coordinates)
#                     Jacobian +=
#                         (shfg_parent[2:end, k2:k2] * reshape(collect(a_node_coor), 1, :))
#                 end
#                 Jacobian_det = det(Jacobian)

#                 @assert Jacobian_det > 0.0 "Determinant of a Jacobian must be >= 0.0 but the element $(ele_tag) has $Jacobian_det."

#                 push!(ivols, Jacobian_det * a_ip.wt)
#                 #
#                 aip_shf_grad_xyz =
#                     [shfg_parent[1:1, :]; inv(Jacobian) * shfg_parent[2:end, :]]
#                 push!(shfg_xyz, SMatrix{1 + dim,num_shf,Float64}(aip_shf_grad_xyz))
#             end
#             elements = [
#                 elements
#                 fele_data_type(
#                     tag = ele_tag,
#                     node_tags = SVector{num_nodes}(ele_nodes),
#                     ip_volumes = ivols,
#                     sfgrad_xyz = shfg_xyz,
#                 )
#             ]  
#             tot_vol += sum(ivols)  # *****************
#         end
#         #
#         push!(element_sets, FiniteElementSet(elset_tag, elements))
#     end 
#     println("Total Volume by brute sum ", tot_vol)   # *****************
#     return element_sets
# end


# function get_felement_sets(
#     element_connectivity::Dict{Int64,Matrix{Int64}},
#     nodal_data::Dict{Int64, Vector{Float64}},
#     elset_tag::Union{Int64,String};
#     num_ip::Int64 = -1,
#     gq_order::Int64 = -1,
#     print_debug_info::Bool = false,
# )::Vector{AbstractFElementSet}
#     #
#     println("Preparing $elset_tag...!")
#     element_sets::Vector{AbstractFElementSet} = AbstractFElementSet[]
#     elements::Vector{AbstractFElement} = AbstractFElement[]
#     num_nodes::Int64 = 0
#     num_shf::Int64 = 0
#     num_felements::Int64 = 0
#     ivols::Vector{Float64} = []
#     shfg_xyz::Vector{SMatrix} = []
#     #
#     #
#     for (ele_type, ele_con) in element_connectivity
#         # @assert ele_type in fe_libray_g "Element type $ele_type (GMSH id) is not implemented as of now..!"
#         empty!(elements)
#         num_nodes = size(ele_con)[1] - 1  # as the first row corresponds to element tag
#         num_felements = size(ele_con)[2]
#         dim, fele_data_type, integration_points =
#             get_felement_properties(ele_type, gq_order, num_ip)
#         num_ip = length(integration_points)
#         num_shf = num_nodes
#         if print_debug_info
#             println(
#                 "For $fele_data_type, \nNumber of elements :: $(num_felements), \nNumber of nodes/element :: $num_nodes, \nIntegration points :: $integration_points,",
#             )
#         end
#         #
#         # get shape functions and its gradients at gauss points in parent coordinate system
#         shf_data_parent_csys = get_shape_functions_info(fele_data_type, integration_points)

#         # get Jacobian and its determinant
#         for (ele_tag, ele_nodes...) in eachcol(ele_con)
#             nodal_coordinates = (nodal_data[i][1:dim] for i in ele_nodes)
#             empty!(ivols)
#             empty!(shfg_xyz)
#             # run over each gauss point
#             for (a_ip, shfg_parent) in shf_data_parent_csys
#                 Jacobian = zeros(Float64, dim, dim)
#                 # run over pairs of coordinates and shape function data

#                 for (k2, a_node_coor) in enumerate(nodal_coordinates)
#                     Jacobian +=
#                         (shfg_parent[2:end, k2:k2] * reshape(collect(a_node_coor), 1, :))
#                 end

#                 Jacobian_det = det(Jacobian)

#                 @assert Jacobian_det > 0.0 "Determinant of a Jacobian must be >= 0.0 but $Jacobian_det observed on element $(ele_tag)"

#                 push!(ivols, Jacobian_det * a_ip.wt)
#                 #
#                 aip_shf_grad_xyz =
#                     [shfg_parent[1:1, :]; inv(Jacobian) * shfg_parent[2:end, :]]
#                 push!(shfg_xyz, SMatrix{1 + dim,num_shf,Float64}(aip_shf_grad_xyz))
#             end
#             elements = [
#                 elements
#                 fele_data_type(
#                     tag = ele_tag,
#                     node_tags = SVector(ele_nodes...),
#                     ip_volumes = ivols,
#                     sfgrad_xyz = shfg_xyz,
#                 )
#             ]
#         end
#         #
#         push!(element_sets, FiniteElementSet(elset_tag, elements))
#     end
#     return element_sets
# end




