# Assumes that a connected graph is given as input
function Counter(edges::Vector{Vector{Int64}})
    # find the indicies of the vertices
    vertices = unique(reduce(vcat,edges))
    
    # find the number of edges attached to each vertex
    edges_per_vertex = [sum([v in e for e in edges]) for v in vertices]
    
    # compute degree of KC
    deg_KC = [-2+e for e in edges_per_vertex]

    # Make binary choice: each edge can either be blown up or not
    h0s = Int64[]
    remaining_nodes = Int64[]
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
                new_degrees = Vector{Int64}([div(d,2) for d in new_degs])
                h0 = H0(Vector{Int64}(vertices), Vector{Vector{Int64}}(new_edges), Vector{Int64}([div(d,2) for d in new_degs]))
                push!(h0s,h0)
                push!(remaining_nodes,length(new_edges))
            end
        end
    end

    # Compute statistics in multiples of 2^b1:
    max = maximum(h0s)
    b1 = length(edges) + 1 - length(vertices)
    total = 0
    res_matrix = zero_matrix(ZZ, length(edges)+2, max+1)
    for i in 1:length(h0s)
        total += 2^b1
        res_matrix[length(edges)+2,h0s[i]+1] += 2^b1
        res_matrix[remaining_nodes[i]+1,h0s[i]+1] += 2^b1
    end

    # Check that we found exactly as many roots as expected
    if total != 2^(2*b1)
        error("Found more or less roots than exist!")
    end

    # Return the global sections
    return res_matrix
end
export Counter


# Compute generic h0
function H0(vertices::Vector{Int64}, edges::Vector{Vector{Int64}}, degrees::Vector{Int64})
    # Find indices of the components with non-negative degree
    indices = findall(>=(0), degrees)

    # Construct a zero matrix
    number_of_local_sections = sum([degrees[i] + 1 for i in indices])
    matrix = zero_matrix(ZZ, length(edges), number_of_local_sections)

    # Identify index in matrix for non-trivial sections on each component
    indices_dict = Dict{Int64, Int64}()
    if length(indices) > 0
        indices_dict[indices[1]] = 1
        index = 1
        for i in 2:length(indices)
            index += degrees[indices[i]] + 1
            indices_dict[indices[i]] = index
        end
    end

    # For each edge/node add entries to the matrix
    for i in 1:length(edges)
        e = edges[i]
        d1 = degrees[e[1]]
        d2 = degrees[e[2]]
        if d1 >= 0
            pos_v = rand(-20:20)
            for j in 0:d1
                matrix[i, indices_dict[e[1]]+j] = pos_v^j
            end
        end
        if d2 >= 0
            pos_v = rand(-20:20)
            for j in 0:d2
                matrix[i, indices_dict[e[2]]+j] = (-1) * pos_v^j
            end
        end
    end

    # Return h0
    return nullspace(matrix)[1]
end
export H0
