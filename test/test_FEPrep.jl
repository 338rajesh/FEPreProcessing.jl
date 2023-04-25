using UnitCellGenerator
using UnitCellModelling
const UCG = UnitCellGenerator
const UCM = UnitCellModelling
using FEPreProcessing
using Materials

println("\t\t Testing FEPreProcessing module \n")

# ============================================
#           UNIT CELL GENERATION
# ============================================
r_mean::Float64 = 2.0
r_std::Float64 = 0.0
rve_size::Float64 = 30.0
#
volume_fractions = Dict{String, Float64}(
    "Circles"       => 0.1,
    "Ellipses"      => 0.1,
    "nLobes"        => 0.05,
    "Reg_Polygons"  => 0.05,
    "Rectangles"    => 0.05,
    "Capsules"      => 0.05,
)

#
ruc_bounds::NTuple{4, Float64} = ((-1.0, -1.0, 1.0, 1.0) .* (0.5 * rve_size))  #TODO getting RUC size
const ruc_info = RUC_data(
    bbox=UCG.BBox2D(ruc_bounds...,),
    ssd_ratio=0.07,
    inclusion_distr=RANDOM,
    periodicity=true,
)
# ===================
#   Inclusions data
# ===================
const circ_inc = Inclusion_data(
    volume_fraction=volume_fractions["Circles"],
    shape=Circle,
    size_params=Dict(:RADIUS => Normal(r_mean, r_std),),
)
const caps_inc = Inclusion_data(
    volume_fraction=volume_fractions["Capsules"],
    shape=Capsule,
    size_params=Dict(:SMJRX => Normal(1.42, 0.0), :SMNRX => Normal(0.925, 0.0),),
)
const lobular_inc = Inclusion_data(
    volume_fraction=volume_fractions["nLobes"],
    shape=nLobeShape,
    size_params=Dict(:NLOBES => 2, :EQRAD => Normal(r_mean, r_std), :LOBE_DIST_FACTOR => 0.5,)
)

const ell_inc = Inclusion_data(
    volume_fraction=volume_fractions["Ellipses"],
    shape=Ellipse,
    size_params=Dict(:SMJRX => Normal(2.0, 0.0), :SMNRX => Normal(1.0, 0.0),),
)
const rect_inc = Inclusion_data(
    volume_fraction=volume_fractions["Rectangles"],
    shape=Rectangle,
    size_params=Dict(:SMJRX => Normal(2.0, 0.0), :SMNRX => Normal(1.0, 0.0), :CRAD => Normal(0.2, 0.0),),
)
inclusions_data = generate_unit_cell(ruc_info, (circ_inc, caps_inc, ell_inc, rect_inc, lobular_inc,), adjust_ruc_bbox=false,)


# ============================================
#       MODELLING and MESHING in GMESH
# ============================================
rve_bounds =(
    inclusions_data["bbox"][1], inclusions_data["bbox"][2], -0.5*r_mean,
    inclusions_data["bbox"][3], inclusions_data["bbox"][4], 0.5*r_mean,
)
rve_inclusions_data = Dict(k => v for (k,v) in inclusions_data if k!="bbox")
udc_3d = UCM.UDC3D(UCG.BBox3D(rve_bounds...,), rve_inclusions_data,)

ruc_model_data = make_unit_cell_model(
    udc_3d,
    mesh_periodicity=true,
    element_types=(:C3D6, :C3D8),
    geom_export_paths=(),
    extr_dir_num_ele=Int64[3,],
    extr_dir_cum_heights=Float64[1.0,],
    extr_dir_recombine_ele=true,
    min_ele_size_factor=1/8,  #FIXME
    max_ele_size_factor=1/4,
    mesh_opt_algorithm="Netgen",
    show_mesh_stats=true,
    show_rve=true,
    node_renum_algorithm="",
)
#
# ============================================
#             FE Pre-Processing
# ============================================


matrix_phase_info = FEPreProcessing.PhaseFiniteElementConnectivity(
    "Matrix Element Set",
    ruc_model_data["mesh_data"]["matrix_element_connectivity"],
    Materials.aluminium_matrix,
    )

inclusion_phase_info = FEPreProcessing.PhaseFiniteElementConnectivity(
    "Inclusions Element Set",
    ruc_model_data["mesh_data"]["inclusions_element_connectivity"],
    Materials.boron_fibre,
)

uc_femodel = make_unit_cell_FEA_model(
    ruc_model_data["mesh_data"]["all_node_tags"],
    ruc_model_data["mesh_data"]["all_node_coordinates"],
    [matrix_phase_info, inclusion_phase_info]
)

total_ele_vol = FEPreProcessing.get_total_volume_of_elements(uc_femodel["element_sets"]; print_info=true)
