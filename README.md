# Finite Element pre-processing of Unit Cell

This Julia module contains methods capable of creating Finite Element model which should be ready for solving. Though, we are trying to write the FE pre-processing code as generic as possible the primary motivation for this module is for Unit Cell Homogenization. This is inevitable, given the vast applications of Finite element analysis.

## Installation

From Julia REPL, execute the following,

```julia
julia> ]add https://github.com/338rajesh/FEPreProcessing.jl
```

## Documentation

`FEPreProcessing.jl` exports a single function `unit_cell_FEModel()` which takes the following set of arguments.

* Positional Arguments
  * `ntags::Vector{Int}`, Node Tags
  * `ncoor::Matrix{Float64}`, Nodal Coordinates matrix. Here each column of `ncoor` should corresponds to a *finite element* node.
  * `phases_data::Dict`, data corresponding to each phase of the unit cell can be given supplied as the following key-value pairs.
    * **key**: phase *finite element* set ID.
    * **value**: A tuple of the folloing two values
      * **phase *finite element* connectivity**, of `Matrix{Int}` type where each column corresponds to a *finite element*. In each of these columns, first element is *finite element* tag and the remaining are *finite element* node tags.
      * ***Material***,

* Keyword Arguments
  * `ϵ::Float64 = 1e-06`, small parameter/tolerance to be used in the modelling process.
  * `numIP::Int = 4`, Number of integration points to be used.
  * `gqOrder::Int = 2`, Gass quadrature order

### Analysis Types

* `AbstractFEAnalysis`
  * `FEA_3Danalysis`
    * `Elastic_3DFEA`
    * `Piezo_Electric_3DFEA`
    * `Piezo_Magnetic_3DFEA`
    * `Magneto_Electro_Elastic_3DFEA`
    * `Thermal_Conductivity_3DFEA`
    * `Thermo_Elastic_3DFEA`
    * `Thermo_Piezo_Electric_3DFEA`
    * `Thermo_Piezo_Magnetic_3DFEA`
    * `Thermo_Magneto_Electro_Elastic_3DFEA`
  * `FEA_2Danalysis`
    * `Elastic_2DFEA`
      * `PlaneStress_2DFEA`
      * `PlaneStrain_2DFEA`

### Usage

* Nodal degrees of freedom for an ananlysis can be obtained using,

```julia
julia> using FEPreProcessing
julia> FEPreProcessing.get_nodal_dof(Elastic_3DFEA)
(:ux, :uy, :uz)
julia> FEPreProcessing.get_nodal_dof(Thermal_Conductivity_3DFEA)
(:T,)
```

* node pairs of the unit cells, assuming model has periodic mesh, can be obtained using `make_rve_node_pairs` method. It takes node tags of `Vector{Int}` type and node coordinates of `Matrix{Float64}` type. Note that, each column of node coordinates matrix corresponds to a node.

```julia
julia> using FEPreProcessing
julia> node_pairs = FEPreProcessing.make_rve_node_pairs(all_node_tags, all_node_coor)
julia> pc_node_tags = [(i.n1, i.n2) for i in node_pairs]
```

### Code (dependency) structure

* unit_cell_FEModel
  * `make_rve_node_pairs()`
    * `get_rve_node_groups()`
    * `node_pairs(::Dict{String,Matrix{IF64}}, ::NTuple{6,Float64})`
      * `node_pairs(::Matrix{IF64}, ::Matrix{IF64}, ::NTuple{3,Float64}, ::String; ::Float64)`
  * `make_fe_node_set()`
  * `get_felement_sets()`
