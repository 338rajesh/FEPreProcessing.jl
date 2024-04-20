using NPZ
using FEPreProcessing

ncoor = npzread(joinpath(dirname(@__FILE__), "sample_nodes.npy"))
# transpose nodes, and make it matrix of Float64
ncoor = convert(Matrix{Float64}, ncoor')
num_nodes = size(ncoor, 2)
ntags = collect(1:num_nodes)

uc_node_set = FEPreProcessing.NodeSet(tags=ntags, coor=ncoor, tag="uc_all_nodes")

# FEPreProcessing.print_node_set(uc_node_set)

node_pairs = FEPreProcessing.make_node_pairs(uc_node_set; small_par=1.0e-6)
FEPreProcessing.print_node_pairs(node_pairs)
