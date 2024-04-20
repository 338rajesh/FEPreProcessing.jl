module FEPreProcessing
    __precompile__()
    
    using Printf

    using LinearAlgebra
    using StructArrays
    using StaticArrays
    using SparseArrays
    using Parameters
    using Materials
    
    include("FEPrepBase.jl") # Defines various custom type definitions and implements functions on these types.
    include("prn_utils.jl")  # Implements printing functions
    include("mesh_data.jl")  # Implements Periodicity checking, making node pairs and finite element and node sets.
    include("node_pairing.jl")  # Implements node pairing
    include("assembly.jl")  # Implements strain displacement matrices, constraint matrices (CSC form)
    include("prepare_material_data.jl")  # Prepares material matrices, according to the specified analysis type.
    include("UnitCellFEModel.jl")  # Makes finite element model of the unit cell, making used of the above steps.
    
    export make_unit_cell_FEA_model
end
