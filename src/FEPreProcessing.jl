module FEPreProcessing
    __precompile__()
    using LinearAlgebra
    #
    using StructArrays
    using StaticArrays
    using SparseArrays
    using Parameters
    #
    using Materials


    const IF64 = Union{Int64, Float64}

    include("FEPrepBase.jl")
    include("mesh_data.jl")
    include("assembly.jl")
    include("prepare_material_data.jl")
    #
    include("UnitCellFEModel.jl")
    export make_unit_cell_FEA_model
    #
end
