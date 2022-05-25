# edges = [[1,2],[2,3],[3,1]]
# edges = [[1,2],[1,2],[1,2]]
# Counter(edges)

# Assumes that a connected graph is given as input
function Counter(edges::Vector{Vector{Int64}})
    # find the indicies of the vertices
    vertices = unique(reduce(vcat,edges))
    
    # find the number of edges attached to each vertex
    edges_per_vertex = [sum([v in e for e in edges]) for v in vertices]
    
    # compute degree of KC
    deg_KC = [-2+e for e in edges_per_vertex]
    
    # print what we computed thus far
    print("vertices: " * string(vertices) * "\n")
    print("edges per vertex: " * string(edges_per_vertex) * "\n")
    print("degree of KC: " * string(deg_KC) * "\n\n")
    
    # Make binary choice: each edge can either be blown up or not
    total_possibilities = []
    for k in 0:length(edges)
        combinations = collect(Combinatorics.combinations(1:length(edges),k))
        for c in combinations
            # compute the new degrees from these blowups
            new_degs = [d for d in deg_KC]
            for edge_index in c
                new_degs[edges[edge_index][1]] = new_degs[edges[edge_index][1]] - 1
                new_degs[edges[edge_index][2]] = new_degs[edges[edge_index][2]] - 1
            end
            
            # check if all degrees are divisible by two
            if all(y->iseven(y),new_degs)
                new_edges = []
                for i in 1:length(edges)
                    if !(i in c)
                        push!(new_edges, edges[i])
                    end
                end
                push!(total_possibilities,[Vector{Int64}(vertices), Vector{Vector{Int64}}(new_edges), Vector{Int64}([div(d,2) for d in new_degs])])
            end
        end
    end

    # Inform the user about the setups to be analyzed
    print("\n######################\n")
    print("Setups to be analyzed:\n")
    print("######################\n\n")
    for p in total_possibilities
        print("Vertices: " * string(p[1]) * "\n")
        print("Remaining edges: " * string(p[2]) * "\n")
        print("Degrees of line bundles: " * string(p[3]) * "\n\n")
    end
    print("\n######################\n\n")

    # Now compute h0 for each of these configuration under generic assumptions
    h0s = [H0(p[1], p[2], p[3]) for p in total_possibilities]
    return [total_possibilities, h0s]
end
export Counter


# Assumes that a connected graph is given as input
function H0(vertices::Vector{Int64}, edges::Vector{Vector{Int64}}, degrees::Vector{Int64})
    # Find the number of edges attached to each vertex
    edges_per_vertex = [sum([v in e for e in edges]) for v in vertices]

    # Identify isolated components
    h0 = 0
    for i in 1:length(vertices)
        if (edges_per_vertex[i] == 0) && (degrees[i] >= 0)
            h0 = h0 + (degrees[i]+1)
        end
    end

    # Treat remaining components


    # Return h0
    return h0
end
export H0
