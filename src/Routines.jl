# Assumes that a connected graph is given as input
function Counter(edges::Vector{Vector{Int64}})
    # Make dictionary with index for each vertex
    vertices = unique(reduce(vcat,edges))
    vert_dict = Dict{Int64, Int64}()
    for i in 1:length(vertices)
        vert_dict[vertices[i]] = i
    end

    # Compute degree of KC
    nodes_per_vertex = fill(0,length(vertices))
    for e in edges
        nodes_per_vertex[vert_dict[e[1]]] += 1
        nodes_per_vertex[vert_dict[e[2]]] += 1
    end
    deg_KC = [-2+e for e in nodes_per_vertex]

    # Initialize variables
    res_matrix = zero_matrix(ZZ, length(edges)+1, sum([deg_KC[k]+1 for k in findall(>=(0), deg_KC)])+1)
    mult = 2^(length(edges) + 1 - length(vertices))
    total = 0

    # Iterate over all roots (below, c are the edges/nodes that we blow up)
    for k in 0:length(edges)
        combinations = Oscar.Hecke.subsets([e for e in 1:length(edges)],k)
        for c in combinations

            # Compute the new degrees from the blowups
            new_degs = [d for d in deg_KC]
            for edge_index in c
                new_degs[vert_dict[edges[edge_index][1]]] = new_degs[vert_dict[edges[edge_index][1]]] - 1
                new_degs[vert_dict[edges[edge_index][2]]] = new_degs[vert_dict[edges[edge_index][2]]] - 1
            end

            # Check if all degrees are divisible by two
            if all(y->iseven(y),new_degs)

                # Compute h0
                new_edges = Vector{Int64}[]
                for i in 1:length(edges)
                    if !(i in c)
                        push!(new_edges, edges[i])
                    end
                end
                h0 = H0(vert_dict, new_edges, [div(d,2) for d in new_degs])

                # Update result matrix and total
                res_matrix[length(new_edges)+1,h0+1] += mult
                total += mult
            end
        end
    end

    # Check that we found exactly as many roots as expected
    if total != mult * mult
        error("Found more or less roots than exist!")
    end

    # Return the result
    return res_matrix
end
export Counter


# Compute generic h0
function H0(vert_dict::Dict{Int64, Int64}, edges::Vector{Vector{Int64}}, degrees::Vector{Int64})

    # Find indices of the components with non-negative degree
    indices = findall(>=(0), degrees)

    # Construct a zero matrix
    number_of_local_sections = sum([degrees[i]+1 for i in indices])
    matrix = zero_matrix(ZZ, length(edges), number_of_local_sections)

    # Identify column index for non-trivial sections on each component
    indices_dict = Dict{Int64, Int64}()
    if length(indices) > 0
        indices_dict[indices[1]] = 1
        index = 1
        for i in 2:length(indices)
            index += degrees[indices[i]] + 1
            indices_dict[indices[i]] = index
        end
    end

    # For each node add entries to the matrix
    for i in 1:length(edges)
        e = edges[i]
        if degrees[vert_dict[e[1]]] >= 0
            pos_v = rand(-20:20)
            for j in 0:degrees[vert_dict[e[1]]]
                matrix[i, indices_dict[vert_dict[e[1]]]+j] += pos_v^j
            end
        end
        if degrees[vert_dict[e[2]]] >= 0
            pos_v = rand(-20:20)
            for j in 0:degrees[vert_dict[e[2]]]
                matrix[i, indices_dict[vert_dict[e[2]]]+j] += (-1) * pos_v^j
            end
        end
    end

    # Return h0
    return nullspace(matrix)[1]
end
export H0
