
const IF64 = Union{Int64, Float64}

# ===================================================================
#                       ERRORS
# ===================================================================

abstract type FEPrepException <: Exception end

struct FEPrepError <: FEPrepException
    msg::AbstractString
end

Base.showerror(io::IO, e::FEPrepError) = begin
    print(io, "FEPrep Error :: " * e.msg)
end 


struct MeshPeriodicityError <: FEPrepException
    node_tag::Int
    group::AbstractString
end

Base.showerror(io::IO, e::MeshPeriodicityError) = begin
    print(io, "MeshPeriodicityError: Unable to find the corresponding node
    for a node, with tag ", e.node_tag, ", located in group ", e.group)
end


# ===================================================================
#                       INTEGRATION POINTS
# ===================================================================

abstract type AbstractIntegrationPoint end

struct IntegrationPoint1D <: AbstractIntegrationPoint
    wt::Float64
    r::Float64
end


struct IntegrationPoint2D <: AbstractIntegrationPoint
    wt::Float64
    r::Float64
    s::Float64
end


struct IntegrationPoint3D <: AbstractIntegrationPoint
    wt::Float64
    r::Float64
    s::Float64
    t::Float64
end


"""It returns tuple of integration points for a given finite elemenet type. One needs to provide either the order of gauss quadrature or number of integration points or both depending on the finite element type.

1. `num_ip`    : Triangular or Tetrahedral regions
2. `gq_order`  : Rectangular or Brick regions
3. `num_ip` and `gq_order` : Prism regions

"""
function get_integration_points(
    fele_type::DataType;
    gq_order::Int64 = -1,
    num_ip::Int64 = -1
)::Tuple{Vararg{AbstractIntegrationPoint}}
    @assert gq_order > 0 || num_ip > 0 "Gauss Quadrature order or Number of integration points depending on the element type."
    if fele_type == T2D2
        @assert gq_order > 0 "Line elements requires order of integration > 1 and it depends on the desired level of accuracy but $gq_order is given"
        return get_1D_integration_points(gq_order)
    elseif fele_type == CPS3
        @assert num_ip > 0 "Triangular region requires nuber of integration points. Valid numbers are 1, 3, 4"
        return get_triangular_region_integration_points(num_ip)
    elseif fele_type == CPS4
        @assert gq_order > 0 "Rectangular regions requires order of integration > 1 and it depends on the desired level of accuracy but $gq_order is given"
        return get_rectangular_region_integration_points(gq_order)
    elseif fele_type == C3D4
        @assert num_ip > 0 "Tetrahdral region requires number of integration points. Valid numbers are 1, 4, 5"
        return get_tetrahedral_region_integration_points(num_ip)
    elseif fele_type == C3D8
        @assert gq_order > 0 "Brick element regions require order of integration > 1 and it depends on the desired level of accuracy but $gq_order is given"
        return get_brick_region_integration_points(gq_order)
    elseif fele_type == C3D6
        @assert (num_ip > 0 && gq_order > 0) "Wedge element regions require order of integration > 1 and number of integration points > 1 which depends on the desired level of accuracy. Given order of integration = $gq_order, number of integration points = $num_ip"
        return get_wedge_region_integration_points(gq_order, num_ip)
    else
        @error "Encounteres unknown finite element type :: $fele_type"
    end
end


"""
Source:: https://pomax.github.io/bezierinfo/legendre-gauss.html
It contains data up to an order of 64.
"""
function get_1D_integration_points(
    order::Int64
)::Tuple{Vararg{IntegrationPoint1D}}
    wt::Vector{Float64} = []
    ξ::Vector{Float64} = []
    if order == 1
        wt, ξ = ([2.0,], [0.0,])
    elseif order == 2
        wt, ξ = ([1.0, 1.0], [-1 / sqrt(3), 1 / sqrt(3)])
    elseif order == 3
        wt, ξ = ([5 / 9, 8 / 9, 5 / 9], [-sqrt(0.6), 0.0, sqrt(0.6)])
    elseif order == 4
        wt, ξ = ([0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538], [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526])
    end
    return tuple((IntegrationPoint1D(aw, aξ) for (aw, aξ) in zip(wt, ξ))...)
end


function get_rectangular_region_integration_points(
    order::Int64,
)::Tuple{Vararg{IntegrationPoint2D}}
    ip_1d::Tuple{Vararg{IntegrationPoint1D}} = get_1D_integration_points(order)
    return tuple(
        (IntegrationPoint2D(a1Dpnti.wt * a1Dpntj.wt, a1Dpnti.r, a1Dpntj.r) for a1Dpnti in ip_1d for a1Dpntj in ip_1d)...
    )
end


function get_brick_region_integration_points(
    order::Int64,
)::Tuple{Vararg{IntegrationPoint3D}}
    ip_1d::Tuple{Vararg{IntegrationPoint1D}} = get_1D_integration_points(order)
    return tuple(
        (IntegrationPoint3D(a1Dpnti.wt * a1Dpntj.wt * a1Dpntk.wt, a1Dpnti.r, a1Dpntj.r, a1Dpntk.r) for a1Dpnti in ip_1d for a1Dpntj in ip_1d for a1Dpntk in ip_1d)...
    )
end


function get_wedge_region_integration_points(
    order::Int64,
    num_ip::Int64,
)::Tuple{Vararg{IntegrationPoint3D}}
    return tuple(
        IntegrationPoint3D(1 / 6, 1 / 6, 1 / 6, -1 / sqrt(3.0)),
        IntegrationPoint3D(1 / 6, 1 / 6, 1 / 6, 1 / sqrt(3.0)),
        IntegrationPoint3D(1 / 6, 1 / 6, 2 / 3, -1 / sqrt(3.0)),
        IntegrationPoint3D(1 / 6, 1 / 6, 2 / 3, 1 / sqrt(3.0)),
        IntegrationPoint3D(1 / 6, 2 / 3, 1 / 6, -1 / sqrt(3.0)),
        IntegrationPoint3D(1 / 6, 2 / 3, 1 / 6, 1 / sqrt(3.0)),
    )
end


"""
returns integration Points of triangular regions, given the number of gauss/integration points.
Source:: Concepts and applications of Finite Element Analysis by Robert D. Cook, Fourth Edition, Page Number: 267
"""
function get_triangular_region_integration_points(
    num_ip::Int64,  # number of integration points
)::Tuple{Vararg{IntegrationPoint2D}}
    wt::Vector{Float64} = []
    ξ::Vector{Float64} = []
    η::Vector{Float64} = []
    if num_ip == 1
        return (IntegrationPoint2D(1 / 2, 1 / 3, 1 / 3),)
    elseif num_ip == 3
        return (
            IntegrationPoint2D(1 / 6, 2 / 3, 1 / 6),
            IntegrationPoint2D(1 / 6, 1 / 6, 1 / 6),
            IntegrationPoint2D(1 / 6, 1 / 6, 2 / 3),
        )
    elseif num_ip == 4
        return (
            IntegrationPoint2D(-27 / 96, 1 / 3, 1 / 3),
            IntegrationPoint2D(25 / 96, 3 / 5, 1 / 5),
            IntegrationPoint2D(25 / 96, 1 / 5, 1 / 5),
            IntegrationPoint2D(25 / 96, 1 / 5, 3 / 5),
        )
    end
end


"""
returns integration Points of Tetrahedral regions, given the number of gauss/integration points.
Source:: Concepts and applications of Finite Element Analysis by Robert D. Cook, Fourth Edition, Page Number: 268
"""
function get_tetrahedral_region_integration_points(
    num_ip::Int64,  # number of integration points
)::Tuple{Vararg{IntegrationPoint3D}}
    wt::Vector{Float64} = []
    ξ::Vector{Float64} = []
    η::Vector{Float64} = []
    if num_ip == 1
        return (
            IntegrationPoint3D(1 / 6, 1 / 4, 1 / 4, 1 / 4),
        )
    elseif num_ip == 4
        a::Float64 = (5 + (3 * sqrt(5))) / 20
        b::Float64 = (5 - sqrt(5)) / 20
        return (
            IntegrationPoint3D(1 / 24, a, b, b),
            IntegrationPoint3D(1 / 24, b, b, b),
            IntegrationPoint3D(1 / 24, b, b, a),
            IntegrationPoint3D(1 / 24, b, a, b),
        )
    elseif num_ip == 5
        return (
            IntegrationPoint3D(3 / 40, 1 / 2, 1 / 6, 1 / 6),
            IntegrationPoint3D(3 / 40, 1 / 6, 1 / 6, 1 / 6),
            IntegrationPoint3D(3 / 40, 1 / 6, 1 / 6, 1 / 2),
            IntegrationPoint3D(3 / 40, 1 / 6, 1 / 2, 1 / 6),
            IntegrationPoint3D(-2 / 15, 1 / 4, 1 / 4, 1 / 4),
        )
    end
end


# ===================================================================
#                       MESH
# ===================================================================


abstract type AbstractNode end
abstract type AbstractElement end

abstract type AbstractFENode <: AbstractNode end
abstract type AbstractFElement <:AbstractElement end

abstract type SolidFElement <: AbstractFElement end
abstract type PlaneFElement <: AbstractFElement end
abstract type LineFElement <: AbstractFElement end

struct FENodePair
    tag::Union{String,Int}
    p_nodetag::Int  # parent node
    c_nodetag::Int  # child node
end

abstract type AbstractFElementSet end

"""Structure for storing a particular type of finite elements"""
@with_kw struct FiniteElementSet <: AbstractFElementSet
    tag::Union{Int64,String}
    elements::Vector{AbstractFElement}
    material::Union{Material, Nothing} = nothing
end


function get_num_ele(fele_set::FiniteElementSet)::Dict{DataType, Int}
    fele_type_summary::Dict{DataType, Int} = Dict{DataType, Int}()
    for a_fele in fele_set.elements
        fele_type = typeof(a_fele)
        if fele_type in keys(fele_type_summary)
            fele_type_summary[fele_type] += 1
        else
            fele_type_summary[fele_type] = 1
        end
    end
    return fele_type_summary
end

function get_num_ele(fele_sets::Vector{FiniteElementSet})
    fele_type_summary::Dict{DataType, Int} = Dict{DataType, Int}()
    for a_fele_set in fele_sets
        aset_fele_type_summary = get_num_ele(a_fele_set)
        for (ak, av) in aset_fele_type_summary
            if ak in keys(fele_type_summary)
                fele_type_summary[ak] += av
            else
                fele_type_summary[ak] = av
            end
        end
    end
    return fele_type_summary
end 

""" PhaseFiniteElementConnectivity

Phase is the union of different types of elements, with same material
type"""
struct PhaseFiniteElementConnectivity
    tag::String
    ele_conn::Dict{Int, Matrix{Int}}
    material::Material
end

raw"""

    Elements local coordinates 

    source: gmsh documentation

                C3D4                                 C3D10

                        v
                      .
                    ,/
                   /
                3                                     3
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,7    '.   `6
        ,/       |     `\                     ,/       9     `\
      ,/         |       `\                 ,/         |       `\
     1-----------'.--------2 --> u         1--------5--'.--------2
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,10
            `\.   '. ,/                           `8.   '. ,/
               `\. |/                                `\. |/
                  `4                                    `4
                     `\.
                        ` w

         C3D8                       C3D20          Hexahedron27:

            v
     4----------3            3----13----2           3----13----2
     |\     ^   |\           |\         |\          |\         |\
     | \    |   | \          | 15       | 14        |15    24  | 14
     |  \   |   |  \         9  \       11 \        9  \ 20    11 \
     |   8------+---7        |   7----19+---6       |   7----19+---6
     |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
     1---+---\--2   |        0---+-8----1   |       0---+-8----1   |
      \  |    \  \  |         \  17      \  18       \ 17    25 \  18
       \ |     \  \ |         10 |        12|        10 |  21    12|
        \|      w  \|           \|         \|          \|         \|
         5----------6            4----16----5           4----16----5


              C3D6                 C3D15                    

                w
                ^
                |
                4                       3                      3
              ,/|`\                   ,/|`\                  ,/|`\
            ,/  |  `\               12  |  13              12  |  13
          ,/    |    `\           ,/    |    `\          ,/    |    `\
         5------+------6         4------14-----5        4------14-----5
         |      |      |         |      8      |        |      8      |
         |    ,/|`\    |         |      |      |        |    ,/|`\    |
         |  ,/  |  `\  |         |      |      |        |  15  |  16  |
         |,/    |    `\|         |      |      |        |,/    |    `\|
        ,|      |      |\        10     |      11       10-----17-----11
      ,/ |      1      | `\      |      0      |        |      0      |
     u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
         |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
         |,/         `\|         |,/         `\|        |,/         `\|
         2-------------3         1------9------2        1------9------2

"""

@with_kw struct C3D4 <: SolidFElement
    tag::Int64
    node_tags::SVector{4,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end

@with_kw struct C3D6 <: SolidFElement
    tag::Int64
    node_tags::SVector{6,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end

@with_kw struct C3D8 <: SolidFElement
    tag::Int64
    node_tags::SVector{8,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end


raw"""
 Triangle:               Triangle6:          Triangle9/10:          Triangle12/15:
 
 v
 ^                                                                   2
 |                                                                   | \
 3                       2                    2                      9   8
 |`\                     |`\                  | \                    |     \
 |  `\                   |  `\                7   6                 10 (14)  7
 |    `\                 5    `4              |     \                |         \
 |      `\               |      `\            8  (9)  5             11 (12) (13) 6
 |        `\             |        `\          |         \            |             \
 1----------2 --> u      0-----3----1         0---3---4---1          0---3---4---5---1
 
 
 CPS4                         CPS8                  CPS9
 
       v
       ^
       |
 3-----------2          4-----7-----3           4-----7-----3
 |     |     |          |           |           |           |
 |     |     |          |           |           |           |
 |     +---- | --> u    8           6           8     9     6
 |           |          |           |           |           |
 |           |          |           |           |           |
 0-----------1          1-----5-----2           1-----5-----2
 
"""

@with_kw struct CPS3 <: PlaneFElement
    tag::Int64
    node_tags::SVector{3,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end


@with_kw struct CPS4 <: PlaneFElement
    tag::Int64
    node_tags::SVector{4,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end


@with_kw struct CPS8 <: PlaneFElement
    tag::Int64
    node_tags::SVector{8,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end


@with_kw struct T2D2 <: LineFElement
    tag::Int64
    node_tags::SVector{2,Int64}
    ip_volumes::Tuple{Vararg{Float64}} = ()
    sfgrad_xyz::Tuple{Vararg{SMatrix}} = ()
end


""" Returns the GEOMETRIC DIMENSION of a finite element 
"""
dim(fe_dtype::DataType)::Int64 = begin
    if fe_dtype in (C3D4, C3D6, C3D8)
        return 3
    elseif fe_dtype in (CPS3, CPS4, CPS8,)
        return 2
    elseif fe_dtype in (T2D2, )
        return 1
    else
        throw(error("Encountered incorrect FE data type $fe "))
    end
end

"""
It returns finite elemenet data type (with ABAQUS convention) for a given gmsh elemenet naming convention.
"""
function get_felement_data_type(
    gmsh_ele_type::Int64
)::DataType
    if gmsh_ele_type == 4
        return C3D4
    elseif gmsh_ele_type == 5
        return C3D8
    elseif gmsh_ele_type == 6
        return C3D6
    elseif gmsh_ele_type == 2
        return CPS3
    elseif gmsh_ele_type == 3
        return CPS4
    elseif gmsh_ele_type == 16
        return CPS8
    elseif gmsh_ele_type == 1
        return T2D2
    end
end

function get_felement_properties(
    element_type::Int64,
    gq_order::Int64 = -1,
    num_ip::Int64 = -1,
)
    ele_data_type::DataType = get_felement_data_type(element_type)
    return (
        dim(ele_data_type),
        ele_data_type,
        get_integration_points(ele_data_type; num_ip = num_ip, gq_order = gq_order),
        )
end


function get_num_nodes_per_ele(e::DataType)
    num_nodes_per_ele::Int = if (e===C3D8) || (e===CPS8)
        8
    elseif e===C3D6
        6
    elseif (e===C3D4) || (e===CPS4)
        4
    elseif e===CPS3
        3
    else
        @error "Element datatype $e is invalid, please provide correct element data type!."
    end
    return num_nodes_per_ele
end


# =============================================================================
#                                 SHAPE FUNCTIONS
# =============================================================================

"""
Returns `SMatrix` of `2x2` shape with each column for a Node/Shape function. First row is shape function evaluations and last row is shape function gradient evaluations at a given integration point.
"""
function shf_T2D2(
    intPoint::IntegrationPoint1D,
)::SMatrix{2,2,Float64}
    ξ::Float64 = intPoint.r
    @assert -1 ≤ ξ ≤ 1 "parent coordinate ξ must be between 0 & 1 but $(ξ) is found."
    N1::Float64 = 0.5 * (1 - ξ)
    N2::Float64 = 0.5 * (1 + ξ)
    dN1::Float64 = -0.5
    dN2::Float64 = 0.5
    return SMatrix{2,2,Float64}([N1 N2; dN1 dN2])
end

"""
Returns `SMatrix` of `3x3` shape with each column for a Node/Shape function. First row is shape function evaluations and remaining two rows are shape function gradient evaluations at a given integration point.
"""
function shf_CPS3(
    intPoint::IntegrationPoint2D,
)::SMatrix{3,3,Float64}
    ξ::Float64 = intPoint.r
    η::Float64 = intPoint.s
    @assert all(0 .≤ (ξ, η) .≤ 1) "parent coordinates have to be between 0 & 1"
    @assert ξ + η ≤ 1.0 "Sum of parent coordinates must be ≤ 1 for Triangular element"
    #
    N1::Float64 = 1 - ξ - η
    N2::Float64 = ξ
    N3::Float64 = η
    dN1::Vector{Float64} = [-1.0; -1.0]
    dN2::Vector{Float64} = [1.0; 0.0]
    dN3::Vector{Float64} = [0.0; 1.0]
    #
    return SMatrix{3,3,Float64}([[N1 N2 N3]; [dN1 dN2 dN3]])
end


"""
Returns `SMatrix` of `3x4` shape with each column for a Node/Shape function. First row is shape function evaluations and remaining two rows are shape function gradient evaluations at a given integration point.
"""
function shf_CPS4(
    intPoint::IntegrationPoint2D,
)::SMatrix{3,4,Float64}
    ξ::Float64 = intPoint.r
    η::Float64 = intPoint.s
    @assert all(-1 .≤ (ξ, η) .≤ 1) "parent coordinates for quad FElement have to be between 0 & 1"
    #
    N1::Float64 = (1 - ξ) * (1 - η)
    N2::Float64 = (1 + ξ) * (1 - η)
    N3::Float64 = (1 + ξ) * (1 + η)
    N4::Float64 = (1 - ξ) * (1 + η)
    dN1::Vector{Float64} = [η - 1.0; ξ - 1.0]
    dN2::Vector{Float64} = [1.0 - η; -1.0 - ξ]
    dN3::Vector{Float64} = [1.0 + η; 1.0 + ξ]
    dN4::Vector{Float64} = [-1.0 - η; 1.0 - ξ]
    #
    return 0.25 * SMatrix{3,4,Float64}([[N1 N2 N3 N4]; [dN1 dN2 dN3 dN4]])
end


"""
Returns `SMatrix` of `4x4` shape with each column for a Node/Shape function. First row is shape function evaluations and remaining three rows are shape function gradient evaluations at a given integration point.
"""
function shf_C3D4(
    intPoint::IntegrationPoint3D,
)::SMatrix{4,4,Float64}
    ξ::Float64 = intPoint.r
    η::Float64 = intPoint.s
    ζ::Float64 = intPoint.t
    @assert all(0 .≤ (ξ, η, ζ) .≤ 1) "parent coordinates have to be between 0 & 1 but $ξ, $η, $ζ are supplied."
    @assert ξ + η + ζ ≤ 1.0 "Sum of parent coordinates must be ≤ 1"
    #
    N1::Float64 = 1 - ξ - η - ζ
    N2::Float64 = ξ
    N3::Float64 = η
    N4::Float64 = ζ
    dN1::Vector{Float64} = [-1.0; -1.0; -1.0]
    dN2::Vector{Float64} = [1.0; 0.0; 0.0]
    dN3::Vector{Float64} = [0.0; 1.0; 0.0]
    dN4::Vector{Float64} = [0.0; 0.0; 1.0]
    #
    return SMatrix{4,4,Float64}([[N1 N2 N3 N4]; [dN1 dN2 dN3 dN4]])
end


"""
Returns `SMatrix` of `4x6` shape with each column for a Node/Shape function. First row is shape function evaluations and remaining three rows are shape function gradient evaluations at a given integration point.
"""
function shf_C3D6(
    intPoint::IntegrationPoint3D,
)::SMatrix{4,6,Float64}
    ξ::Float64 = intPoint.r
    η::Float64 = intPoint.s
    ζ::Float64 = intPoint.t
    @assert all(0 .≤ (ξ, η) .≤ 1) "parent coordinates ξ, η have to be between 0 & 1"
    @assert -1 ≤ ζ ≤ 1 "Parent coordinate ζ along the axis of the wedge element must be between -1 and 1"
    @assert ξ + η ≤ 1.0 "Sum of parent coordinates must be ≤ 1"
    #
    kn::Float64, kp::Float64  = (1.0 - ζ) * 0.5, (1.0 + ζ) * 0.5
    N1::Float64 = (1.0 - ξ - η) * kn
    N2::Float64 = ξ * kn
    N3::Float64 = η * kn
    N4::Float64 = (1.0 - ξ - η) * kp
    N5::Float64 = ξ * kp
    N6::Float64 = η * kp
    dN1::Vector{Float64} = [-kn; -kn; (-1.0 + ξ + η) * 0.5]
    dN2::Vector{Float64} = [kn; 0.0; -ξ * 0.5]
    dN3::Vector{Float64} = [0.0; kn; -η * 0.5]
    dN4::Vector{Float64} = [-kp; -kp; (1.0 - ξ - η) * 0.5]
    dN5::Vector{Float64} = [kp; 0.0; ξ * 0.5]
    dN6::Vector{Float64} = [0.0; kp; η * 0.5]
    #
    return SMatrix{4,6,Float64}([[N1 N2 N3 N4 N5 N6]; [dN1 dN2 dN3 dN4 dN5 dN6]])
end


"""
Returns `SMatrix` of `4x8` shape with each column for a Node/Shape function. First row is shape function evaluations and remaining three rows are shape function gradient evaluations at a given integration point.

"""
function shf_C3D8(
    intPoint::IntegrationPoint3D,
)::SMatrix{4,8,Float64}
    ξ::Float64 = intPoint.r
    η::Float64 = intPoint.s
    ζ::Float64 = intPoint.t
    @assert all(-1 .≤ (ξ, η, ζ) .≤ 1) " All parent coordinates have to be between 0 & 1"
    N1::Float64 = (1.0 - ξ) * (1.0 - η) * (1.0 - ζ)
    N2::Float64 = (1.0 + ξ) * (1.0 - η) * (1.0 - ζ)
    N3::Float64 = (1.0 + ξ) * (1.0 + η) * (1.0 - ζ)
    N4::Float64 = (1.0 - ξ) * (1.0 + η) * (1.0 - ζ)
    N5::Float64 = (1.0 - ξ) * (1.0 - η) * (1.0 + ζ)
    N6::Float64 = (1.0 + ξ) * (1.0 - η) * (1.0 + ζ)
    N7::Float64 = (1.0 + ξ) * (1.0 + η) * (1.0 + ζ)
    N8::Float64 = (1.0 - ξ) * (1.0 + η) * (1.0 + ζ)
    dN1::Vector{Float64} = [(η - 1.0) * (1.0 - ζ); (ξ - 1.0) * (1.0 - ζ); (η - 1.0) * (1.0 - ξ)]
    dN2::Vector{Float64} = [(1.0 - η) * (1.0 - ζ); (ξ + 1.0) * (ζ - 1.0); (η - 1.0) * (1.0 + ξ)]
    dN3::Vector{Float64} = [(1.0 + η) * (1.0 - ζ); (ξ + 1.0) * (1.0 - ζ); (-η - 1.0) * (1.0 + ξ)]
    dN4::Vector{Float64} = [(1.0 + η) * (ζ - 1.0); (1.0 - ξ) * (1.0 - ζ); (η + 1.0) * (ξ - 1.0)]
    dN5::Vector{Float64} = [(η - 1.0) * (1.0 + ζ); (ξ - 1.0) * (1.0 + ζ); (1.0 - ξ) * (1.0 - η)]
    dN6::Vector{Float64} = [(1.0 - η) * (1.0 + ζ); (-ξ - 1.0) * (1.0 + ζ); (1.0 - η) * (1.0 + ξ)]
    dN7::Vector{Float64} = [(1.0 + η) * (1.0 + ζ); (ξ + 1.0) * (1.0 + ζ); (η + 1.0) * (1.0 + ξ)]
    dN8::Vector{Float64} = [(-1.0 - η) * (1.0 + ζ); (1.0 - ξ) * (1.0 + ζ); (1.0 - ξ) * (1.0 + η)]
    return 0.125 .* SMatrix{4,8,Float64}([[N1 N2 N3 N4 N5 N6 N7 N8]; [dN1 dN2 dN3 dN4 dN5 dN6 dN7 dN8]])
end


# ===================================>

"""
    get_shape_functions_info() -> Dict{AbstractIntegrationPoint, SMatrix{1+dim, num_shf, Float64}}

# Args
- `**element_type**: DataType`
- `**integration_points**:Tuple{Vararg{AbstractIntegrationPoint}}`

"""
function get_shape_functions_info(
    element_type::DataType,
    integration_points::Tuple{Vararg{AbstractIntegrationPoint}};
)::Dict{AbstractIntegrationPoint,SMatrix}
    #
    shape_functions_info = begin
        if element_type == CPS3
            shf_CPS3
        elseif element_type == CPS4
            shf_CPS4
        elseif element_type == C3D4
            shf_C3D4
        elseif element_type == C3D6
            shf_C3D6
        elseif element_type == C3D8
            shf_C3D8
        elseif element_type == T2D2
            shf_T2D2
        else
            throw(FEPrepError("Unable to get the shape functions for the element type $element_type"))
        end    
    end
    return Dict((
        a_int_pnt => shape_functions_info(a_int_pnt)
    ) for a_int_pnt in integration_points)
end




# ===================================================================
#                       ANALYSIS
# ===================================================================

abstract type AbstractFEAnalysis end

abstract type FEA_3Danalysis <: AbstractFEAnalysis end
abstract type FEA_2Danalysis <: AbstractFEAnalysis end
#
struct Elastic_3DFEA <: FEA_3Danalysis end
struct Piezo_Electric_3DFEA <: FEA_3Danalysis end
struct Piezo_Magnetic_3DFEA <: FEA_3Danalysis end
struct Magneto_Electro_Elastic_3DFEA <: FEA_3Danalysis end
struct Thermal_Conductivity_3DFEA <: FEA_3Danalysis end
#
struct Thermo_Elastic_3DFEA <: FEA_3Danalysis end
struct Thermo_Piezo_Electric_3DFEA <: FEA_3Danalysis end
struct Thermo_Piezo_Magnetic_3DFEA <: FEA_3Danalysis end
struct Thermo_Magneto_Electro_Elastic_3DFEA <: FEA_3Danalysis end
#
abstract type Elastic_2DFEA <: FEA_2Danalysis end
struct PlaneStress_2DFEA <: Elastic_2DFEA end
struct PlaneStrain_2DFEA <: Elastic_2DFEA end
struct Thermo_Elastic_PlaneStrain_2DFEA <: Elastic_2DFEA end
"""
Returns dof associated with each node for a given analysis and dimension
"""
function get_nodal_dof(
    analysis_type::DataType,
)::Tuple{Vararg{Symbol}}
    if analysis_type in (Elastic_3DFEA, Thermo_Elastic_3DFEA)
        return (:ux, :uy, :uz,)
    elseif analysis_type <: Elastic_2DFEA
        return (:ux, :uy,)
    elseif analysis_type === Thermal_Conductivity_3DFEA
        return (:T,)
    elseif analysis_type === Piezo_Electric_3DFEA
        return (:ux, :uy, :uz, :ϕe,)
    elseif analysis_type === Magneto_Electro_Elastic_3DFEA
        return (:ux, :uy, :uz, :ϕe, :ϕm,)
    end
end

"""
    get_material_matrix_size(analysis_type) -> Tuple(Int, Int)
    
"""
function get_material_matrix_size(
    analysis_type::DataType,
)
    if analysis_type in (Elastic_3DFEA, Thermo_Elastic_3DFEA, )
        return (6, 6)
    elseif (analysis_type <: Elastic_2DFEA) || (analysis_type <: Thermal_Conductivity_3DFEA)
        return (3, 3)
    elseif analysis_type == Piezo_Electric_3DFEA
        return (9, 9)
    elseif analysis_type == Magneto_Electro_Elastic_3DFEA
        return (12, 12)
    else
        @error "Encountered an unknown analysis type for returning material matrix size"
    end
end


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                         Archives/Dev code             
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# abstract type AbstractMesh end
# abstract type AbstractFEMesh <: AbstractMesh end

# struct Mesh <: AbstractMesh
#     nodes::StructVector{AbstractNode}
#     elements::StructVector{AbstractElement}
# end

# struct FEMesh <: AbstractFEMesh
#     nodes::StructVector{AbstractFENode}
#     elements::StructVector{AbstractFElement}
# end


# @with_kw struct FENode <: AbstractFENode
#     tag::Int
#     x::Float64
#     y::Float64 = 0.0
#     z::Float64 = 0.0
# end

# struct C3D{N} <: SolidFElement
#     tag::Int
#     nodes::SVector{N, FENode}
# end

# struct FiniteElement{N} <: AbstractFElement
#     tag::Int64
#     nodes::SVector{N, FENode}
#     ip_volumes::SVector{N, Float64}
#     shf_grad_xyz::SVector{N, SMatrix}
# end

# """
#     get_felement_data_type() -> DataType

#  Returns finite elemenet data type converted to a common data type used in the
#  present code.
 
#  # Args
#  - `element_ID::Union{Int, String, Symbol}`
#  - `convention::AbstractString`, can be either "GMSH" or "ABAQUS"
 
#  # KWArgs
#  - `sos::String="PS"`


# """
# function get_felement_data_type(
#     element_ID::Union{Int, String, Symbol},
#     convention::String;  # 
#     sos::String="PS",
# )::DataType
#     convention = uppercase(convention)

#     if convention=="GMSH"
#         if element_ID == 4
#             return C3D{4}
#         elseif element_ID == 5
#             return C3D{8}
#         elseif element_ID == 6
#             return C3D{6}
#         elseif element_ID == 2
#             if uppercase(sos)=="PS"
#                 return CPS{3}
#             elseif uppercase(sos)=="PE"
#                 return CPE{3}
#             end
#         elseif element_ID == 3
#             if uppercase(sos)=="PS"
#                 return CPS{4}
#             elseif uppercase(sos)=="PE"
#                 return CPE{4}
#             end
#         elseif element_ID == 16
#             if uppercase(sos)=="PS"
#                 return CPS{8}
#             elseif uppercase(sos)=="PE"
#                 return CPE{8}
#             end
#         elseif element_ID == 1
#             return T2D{2}
#         end
#     elseif convention=="ABAQUS"
#         if element_ID == "C3D4"
#             return C3D{4}
#         elseif element_ID == "C3D8"
#             return C3D{8}
#         elseif element_ID == "C3D6"
#             return C3D{6}
#         elseif element_ID == "CPS3"
#             return CPS{3}
#         elseif element_ID == "CPE3"
#             return CPE{3}
#         elseif element_ID == "CPS4"
#             return CPS{4}
#         elseif element_ID == "CPE4"
#             return CPE{4}
#         elseif element_ID == "CPS8"
#             return CPS{8}
#         elseif element_ID == "CPE8"
#             return CPE{8}
#         elseif element_ID == "T2D2"
#             return T2D{2}
#         end
#     else
#         throw(FEPrepError("At present, only gmsh and ABAQUS based element codes are handled "))
#     end    
# end

# #
# gmsh_ele_type(et::T2D2)::Int = 1
# #
# gmsh_ele_type(et::CPS{3})::Int = 2
# gmsh_ele_type(et::CPS{4})::Int = 3
# gmsh_ele_type(et::CPS{8})::Int = 16
# #
# gmsh_ele_type(et::C3D{4})::Int = 4
# gmsh_ele_type(et::C3D{6})::Int = 6
# gmsh_ele_type(et::C3D{8})::Int = 5



# """
# Returns SMatrix of adjacent local node labels for each local node of the
# given element type.

# Each column holds a local node neighbours in order

# """
# function adjacent_local_node_indices(
#     el_dtype::DataType
# )::SMatrix
#     if el_dtype == T2D{2}
#         return SMatrix{1,2,Int64}([2 1])
#     elseif el_dtype in (CPS{3}, CPE{3}, )
#         return SMatrix{2,3,Int64}([[2, 3] [3, 1] [1, 2]])
#     elseif el_dtype in (CPS{4}, CPE{4}, )
#         return SMatrix{2,4,Int64}([[2, 4] [3, 1] [4, 2] [1, 3]])
#     elseif el_dtype in (CPS{8}, CPE{8}, )
#         return SMatrix{2,8,Int64}(
#             [[5, 8] [5, 6] [6, 7] [7, 8] [1, 2] [2, 3] [3, 4] [4, 1]]
#         )
#     elseif el_dtype == C3D{4}
#         return SMatrix{3,4,Int64}([[2, 3, 4] [3, 1, 4] [1, 2, 4] [1, 2, 3]])
#     elseif el_dtype == C3D{6}
#         return SMatrix{3,6,Int64}(
#             [[2, 3, 4] [1, 3, 5] [1, 2, 6] [1, 5, 6] [2, 4, 6] [3, 4, 5]]
#         )
#     elseif el_dtype == C3D{8}
#         return SMatrix{3,8,Int64}(
#             [[2, 4, 5] [3, 1, 6] [7, 2, 4] [1, 8, 3] [1, 8, 6] [2, 5, 7] [3, 6, 8] [4, 5, 7]]
#         )
#     end
# end
