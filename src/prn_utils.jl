function print_node_set(nset::NodeSet)
    println("\n" * repeat("-", 50))
    println("\nNode Set: ", nset.tag)
    println("\nTags: ", nset.tags)
    println("\nCoordinates: ", nset.coor)
    println("\n" * repeat("-", 50) * "\n")
end

function print_node_sets_summary(n_sets::Vector{NodeSet})
    println("\n" * repeat("-", 80) * "\nNode Pairs:")
    @printf(  # Print the header
        "%20s ; %20s ; %20s\n", 
        "Tag", 
        "#Tags", 
        "Coor Size",
        )
    println("\n" * repeat("-", 80))
    for a_n_set in n_sets  # Print the data rows
        @printf(
            "%20s; %20s; %20s;\n", 
            a_n_set.tag, 
            length(a_n_set.tags),
            size(a_n_set.coor),
        )
    end
    println("\n" * repeat("-", 80) * "\n")
end

function print_node_pair(n_pair::NodePair)
    println("\n" * repeat("-", 50))
    println("\nNode Pair: ", n_pair.tag)
    println("\nParent Tag: ", n_pair.parent_tag)
    println("\nChild Tag: ", n_pair.child_tag)
    println("\nParent Coordinates: ", n_pair.parent_coor)
    println("\nChild Coordinates: ", n_pair.child_coor)
    println("\n" * repeat("-", 50) * "\n")    
end


function print_node_pairs(n_pairs::Vector{NodePair}; tag_width::Int=20, coor_width::Int=15)
    println("\n" * repeat("-", 80) * "\nNode Pairs:")
    @printf(  # Print the header
        "%20s ; %10s ; %10s ; %15s ; %15s\n", 
        "Tag", 
        "Parent Tag", 
        "Child Tag", 
        "Parent Coor", 
        "Child Coor"
        )
    println("\n" * repeat("-", 80))
    for a_pair in n_pairs  # Print the data rows
        @printf(
            "%20s; %10s; %10s; %5f %5f %5f ; %5f %5f %5f\n", 
            a_pair.tag, a_pair.parent_tag, a_pair.child_tag, 
            a_pair.parent_coor[1], a_pair.parent_coor[2], a_pair.parent_coor[3],
            a_pair.child_coor[1], a_pair.child_coor[2], a_pair.child_coor[3],
        )
    end
    println("\n" * repeat("-", 80) * "\n")
end