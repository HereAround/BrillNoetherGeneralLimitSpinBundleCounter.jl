# edges = [[1,2],[2,3],[3,1]]
# Counter(edges)

# Assumes that a connected graph is given as input
function Counter(edges::Vector{Vector{Int64}})
    # find the indicies of the vertices
    vertices = unique(reduce(vcat,edges))
    
    # find the number of edges attached to each vertex
    edges_per_vertex = [sum([v in e for e in edges]) for v in vertices]
    
    # compute degree of KC
    deg_KC = [-2+e for e in edges_per_vertex]
    
    # Make binary choice: each edge can either be blown up or not
    
    
    return [vertices, edges_per_vertex, deg_KC]
end
export Counter