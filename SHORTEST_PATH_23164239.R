# Dijksta's Algorithm
# Define the adjacency matrix for the graph
# Inf denotes no direct path between the nodes
graph <- matrix(c(0, 6, Inf, 1, Inf,    # A to others
                  Inf, 0, 5, 2, Inf,  # B to others
                  Inf, Inf, 0, Inf, 5,  # C to others
                  Inf, 2, Inf, 0, 1,    # D to others
                  Inf, Inf, Inf, Inf, 0), # E to others
                nrow = 5, byrow = TRUE, dimnames = list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E")))

# Dijkstra's algorithm implementation in R
dijkstra <- function(graph, start_vertex) {
  num_vertices <- nrow(graph)
  dist <- rep(Inf, num_vertices)  # Initialize distances to infinity
  dist[start_vertex] <- 0         # Distance to self is zero
  previous <- rep(NA, num_vertices)  # Initialize previous vertices
  unvisited <- 1:num_vertices    # All vertices are initially unvisited
  
  # Priority queue operation using a simple vector and extraction by min distance
  while(length(unvisited) > 0) {
    current <- unvisited[which.min(dist[unvisited])]  # Find the vertex with the smallest distance
    unvisited <- unvisited[unvisited != current]      # Remove it from unvisited
    
    # Relaxation step
    neighbors <- which(!is.infinite(graph[current, ]))  # Find all neighbors (connected vertices)
    for (vertex in neighbors) {
      if (vertex %in% unvisited) {  # Only consider unvisited vertices
        alt <- dist[current] + graph[current, vertex]
        if (alt < dist[vertex]) {
          dist[vertex] <- alt
          previous[vertex] <- current
        }
      }
    }
  }
  
  list(distance = dist, previous = previous)
}

# Run Dijkstra from vertex 'A' (index 1)
result <- dijkstra(graph, start_vertex = 1)

# Output results
cat("Distances from A:\n")
print(result$distance)
cat("\nPrevious vertices in shortest path from A:\n")
print(result$previous)


## BELLMAN FORD
# Define the adjacency matrix for the graph
graph <- matrix(c(0, 6, Inf, 1, Inf,    # A to B, A to D
                  Inf, 0, 5, -2, Inf,   # B to C, B to D
                  Inf, Inf, 0, Inf, -5, # C to E
                  Inf, Inf, Inf, 0, 2,   # D to E
                  Inf, Inf, Inf, Inf, 0), # E
                nrow = 5, byrow = TRUE, 
                dimnames = list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E")))

# Bellman-Ford algorithm function
bellman_ford <- function(graph, source) {
  num_vertices <- nrow(graph)
  dist <- rep(Inf, num_vertices)
  pred <- rep(NA, num_vertices)
  dist[source] <- 0
  
  for (i in 1:(num_vertices-1)) {
    for (u in 1:num_vertices) {
      for (v in 1:num_vertices) {
        if (graph[u, v] != Inf && dist[u] + graph[u, v] < dist[v]) {
          dist[v] <- dist[u] + graph[u, v]
          pred[v] <- u
        }
      }
    }
  }
  
  # Check for negative weight cycles
  for (u in 1:num_vertices) {
    for (v in 1:num_vertices) {
      if (graph[u, v] != Inf && dist[u] + graph[u, v] < dist[v]) {
        cat("Graph contains a negative weight cycle\n")
        return(list(distances = NULL, predecessors = NULL))
      }
    }
  }
  
  return(list(distances = dist, predecessors = pred))
}

# Run Bellman-Ford from vertex 'A' (index 1)
result <- bellman_ford(graph, source = 1)

# Output results
if (!is.null(result$distances)) {
  cat("Distances from A:\n")
  print(result$distances)
  cat("\nPaths from A:\n")
  for (i in seq_along(result$distances)) {
    if (!is.infinite(result$distances[i])) {
      path <- get_path(result$predecessors, i)
      cat(paste(names(result$distances)[i], ":", paste(names(result$distances)[path], 
                                                       collapse = " -> "), "\n"))
    } else {
      cat(paste(names(result$distances)[i], ": Unreachable\n"))
    }
  }
} else {
  cat("The graph contains a negative weight cycle and no solution is possible.\n")
}

## FLOYD WARSHALL ALGORITHM
# Number of nodes
V <- 4
# Initialize the distance matrix with Inf where there is no direct edge
Dist <- matrix(Inf, nrow = V, ncol = V)

# Distance from a node to itself is 0
diag(Dist) <- 0

# Define the edges with their weights
edges <- list(c(1, 3, -2), c(2, 1, 4), c(2, 3, 3), c(4, 2, -1), c(3, 4, 2))

# Update the matrix with the given weights for direct edges
for (edge in edges) {
  i <- edge[1]
  j <- edge[2]
  weight <- edge[3]
  Dist[i, j] <- weight
}
# Floyd-Warshall algorithm
for (k in 1:V) {
  for (i in 1:V) {
    for (j in 1:V) {
      if (Dist[i, j] > Dist[i, k] + Dist[k, j]) {
        Dist[i, j] <- Dist[i, k] + Dist[k, j]
      }
    }
  }
}

# Print the final distance matrix
print(Dist)



## Dijkstra's Algorithm
# Define the adjacency matrix for the graph
# Here, Inf represents no direct connection between nodes
graph <- matrix(c(0, 6, Inf, 1, Inf,
                  Inf, 0, 5, 2, Inf,
                  Inf, Inf, 0, Inf, -5,
                  Inf, Inf, 8, 0, 2,
                  Inf, Inf, Inf, Inf, 0),
                nrow = 5, ncol = 5, byrow = TRUE)

rownames(graph) <- colnames(graph) <- c("A", "B", "C", "D", "E")

# Dijkstra's algorithm function
dijkstra <- function(graph, start_vertex) {
  num_vertices <- nrow(graph)
  dist <- rep(Inf, num_vertices)  # Distance from start
  previous <- rep(NA, num_vertices)  # Previous nodes in optimal path
  dist[start_vertex] <- 0  # Distance to self is 0
  unvisited <- list(vertices = 1:num_vertices)  # All nodes initially unvisited
  
  # Simple priority queue operation (extract min)
  extract_min <- function() {
    min_vertex <- unvisited$vertices[which.min(dist[unvisited$vertices])]
    unvisited$vertices <- unvisited$vertices[unvisited$vertices != min_vertex]
    return(min_vertex)
  }
  
  while (length(unvisited$vertices) > 0) {
    u <- extract_min()
    # Look at each vertex v in the graph adjacent to u
    for (v in 1:num_vertices) {
      if (graph[u, v] != Inf) {  # There is a direct connection
        alt <- dist[u] + graph[u, v]
        # Relax (u, v)
        if (alt < dist[v]) {
          dist[v] <- alt
          previous[v] <- u
        }
      }
    }
  }
  
  list(distance = dist, previous = previous)
}

# Running Dijkstra's algorithm from vertex 1 (A)
result <- dijkstra(graph, start_vertex = 1)

# Print results
print("Distances:")
print(result$distance)
print("Paths (Previous Nodes):")
print(result$previous)

