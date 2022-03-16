module FEPreProcessing
    using LinearAlgebra
    using StructArrays
    using StaticArrays
    using Parameters

    const IF64 = Union{Int64, Float64}

    include("FEPrepBase.jl")
    include("mesh_data.jl")
    include("assembly.jl")
    include("UnitCellFEModel.jl")
        
    export unit_cell_FEModel

end
