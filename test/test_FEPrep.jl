using UnitCellGenerator
using UnitCellModelling
using FEPreProcessing
using Materials
const UCG = UnitCellGenerator
const UCM = UnitCellModelling

working_folder = mkpath(joinpath(homedir(), "WorkOuts", "scratch"))
println("\t\t Finding Effective properties of composite material using RUC of circular fibres \n")

# ============================================
#           UNIT CELL GENERATION
# ============================================
const ruc_bounds = ((-1.0, -1.0, 1.0, 1.0) .* 10.0)  #TODO getting RUC size
const ruc_info = RUC_data(
    bbox = UCG.BBox2D(ruc_bounds...,),
    ssd_ratio = 0.07,
    inclusion_distr = RANDOM,
    periodicity = true,
)
const circ_inc = Inclusion_data(
    volume_fraction = 0.65,
    shape = Circle,
    size_params = Dict(:RADIUS => Normal(2.0, 0.0),),
)
const caps_inc = Inclusion_data(
    volume_fraction = 0.2,
    shape = Capsule,
    size_params = Dict(:SMJRX => Normal(1.42, 0.0), :SMNRX => Normal(0.925, 0.0),),
)
const elliptical_inclusions_data = Inclusion_data(
    volume_fraction = 0.20,
    shape = Ellipse,
    size_params = Dict(:SMJRX => Normal(2.0, 0.0), :SMNRX => Normal(1.0, 0.0),),
)
const rectangular_inclusions_data = Inclusion_data(
    volume_fraction = 0.30,
    shape = Rectangle,
    size_params = Dict(:SMJRX => Normal(2.0, 0.0), :SMNRX => Normal(1.0, 0.0), :CRAD => Normal(0.2, 0.0),),
)
inclusions_data = generate_unit_cell(ruc_info, (circ_inc,),)

# ============================================
#       MODELLING and MESHING in GMESH
# ============================================
udc_3d = UCM.UDC3D(
    UCG.BBox3D(
        ruc_bounds[1], ruc_bounds[2], -0.5,
        ruc_bounds[3], ruc_bounds[4], 0.5,
        ),
    inclusions_data,) 
ruc_model_data = make_unit_cell_model(
    udc_3d,
    mesh_periodicity=true, 
    element_types = (:C3D6, :C3D8),
    geom_export_paths=(joinpath(homedir(), "test.msh"),),
    extr_dir_num_ele = Int64[3,],
    extr_dir_cum_heights = Float64[1.0,],
    extr_dir_recombine_ele = true,
    min_ele_size_factor = 0.20,  #FIXME
    max_ele_size_factor = 0.24,
    mesh_opt_algorithm = "Netgen",
    show_mesh_stats = true,
    show_rve = false,
)
#
# ============================================
#             FE Pre-Processing
# ============================================


matrix_phase_info = FEPreProcessing.PhaseFiniteElementConnectivity(
    "Matrix Element Set ",
    ruc_model_data["mesh_data"]["matrix_element_connectivity"],
    Materials.aluminium_matrix,
    )

inclusion_phase_info = FEPreProcessing.PhaseFiniteElementConnectivity(
    "Inclusions Element Set ",
    ruc_model_data["mesh_data"]["inclusions_element_connectivity"],
    Materials.boron_fibre,
)

uc_femodel = make_unit_cell_FEA_model(
    ruc_model_data["mesh_data"]["all_node_tags"],
    ruc_model_data["mesh_data"]["all_node_coordinates"],
    [matrix_phase_info, inclusion_phase_info]
)

total_ele_vol = FEPreProcessing.get_total_volume_of_elements(uc_femodel["element_sets"]; print_info=true)
