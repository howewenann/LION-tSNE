# Helper functions -----------------------------------------------------------------------

Radius_Calculator <- function(data, percentile){
   
    # Calculate Euclidean for nearest neighbout (1-NN)
    one_NN_list = RANN::nn2(data, k = 2)
    dist_old <- one_NN_list$nn.dists[, 2]
    
    # Get global radius percentile
    radius1 <- quantile(dist_old, probs = percentile)
    
    return(radius1)
}


# Helper function to perform IDW interpolation
Local_IDW_interpolation <- function(X_new, X_train, Y_train, neighbor_indices, p){
    
    # Get neighbours
    X_train_neighbour <- X_train[neighbor_indices, ]
    Y_train_neighbour <- Y_train[neighbor_indices, ]
    
    # Calculate ||x - xi||
    X_dist <- as.numeric(as.matrix(pdist::pdist(X = X_train_neighbour, Y = X_new)))
    
    # Calculate weights wi = ||x - xi|| ^ -p / sum(||x - xi|| ^ -p)
    wi <- (X_dist ^ (-p)) / sum(X_dist ^ (-p))
    
    # Get interpolated Y
    Y_interpolated <- colSums(wi * Y_train_neighbour)
    
    # convert to dataframe before output
    Y_interpolated <- data.table::data.table(t(Y_interpolated))
    
    # Return interpolated Y
    return(Y_interpolated)
    
}


# Helper function to perform single neighbour placement
Single_neighbour_placement <- function(Y_train, y_radius_close, neighbor_indices){
    
    # Get Y neighbours
    Y_train_neighbour <- Y_train[neighbor_indices, ]
    
    # Generate random magnitude
    random_dist = runif(n = 1, min = 0, max = y_radius_close)
    
    # Generate random angles depending on dimension and create randompoint y
    if(dim(Y_train_neighbour)[2] == 2){
        
        random_angle = runif(n = 1, min = 0, max = 2 * pi)
        
        Y_x1 <- random_dist * cos(random_angle)
        Y_x2 <- random_dist * sin(random_angle)
        
        Y_new_random <- c(Y_x1, Y_x2)
        
    } else if (dim(Y_train)[2] == 3){
        
        random_angle1 = runif(n = 1, min = 0, max = 2 * pi)
        random_angle2 = runif(n = 1, min = 0, max = 2 * pi)
        
        Y_x1 <- random_dist * cos(random_angle1) * cos(random_angle2)
        Y_x2 <- random_dist * sin(random_angle1)
        Y_x3 <- random_dist * cos(random_angle1) * sin(random_angle2)
        
        Y_new_random <- c(Y_x1, Y_x2, Y_x3)
        
    }
    
    Y_new <- Y_train_neighbour + Y_new_random
    Y_new <- data.table::data.table(Y_new)
    
    return(Y_new)
}


# Helper function to pre compute outlier positions - cells which do not contain points
PreCompute_Outlier_Positions <- function(Y_train, y_radius){
    
    y_min <- unlist(lapply(Y_train, min))
    y_max <- unlist(lapply(Y_train, max))
    
    original_cell_nums <- list()
    adjusted_cell_sizes <- list()
    
    # Number of cells per dimension and cell sizes
    for(i in 1:dim(Y_train)[2]){
        
        original_cell_nums_temp <- floor((y_max[i] - y_min[i]) / (2 * y_radius)) 
        adjusted_cell_sizes_temp <- (y_max[i] - y_min[i]) / original_cell_nums_temp
        
        original_cell_nums[[i]] <- original_cell_nums_temp
        adjusted_cell_sizes[[i]] <- adjusted_cell_sizes_temp
    }
    
    cell_df <- expand.grid(lapply(original_cell_nums, function(x) (seq(x) - 1)))
    cell_list <- split(cell_df, seq(nrow(cell_df)))
    
    # Calculate bounds for each dimension
    available_cells <- list()
    adjusted_cell_sizes_vec <- unlist(adjusted_cell_sizes)
    
    for(i in 1:length(cell_list)){
        
        cell_bounds_min <- unlist(y_min + cell_list[[i]] * adjusted_cell_sizes_vec)
        cell_bounds_max <- unlist(y_min + (cell_list[[i]] + 1) * adjusted_cell_sizes_vec)
        
        # Check if there are samples within cells
        samples_in_cell <- rep(TRUE, dim(Y_train)[1])
        
        for(j in 1:length(cell_bounds_min)){
            
            samples_in_cell <- samples_in_cell & Y_train[, j, with = F] >= cell_bounds_min[j] & Y_train[, j, with = F] <= cell_bounds_max[j]
            
        }
        
        if(all(!samples_in_cell) == T){
            
            available_cells[[i]] <- cell_list[[i]]
            
        }
        
    }
    
    # Remove NULLs from list
    available_cells <- data.table::rbindlist(available_cells)
    available_cells <- split(available_cells, seq(nrow(available_cells)))
    
    return(available_cells)
}


# Helper function to group outliers according to distance using hclust
Group_Outliers <- function(outliers1, x_radius){
    
    # Get hclust tree
    hclust_tree <- fastcluster::hclust.vector(outliers1)
    
    # Get outlier group by cutting the tree at a given height (dist between points)
    outlier_groups <- cutree(hclust_tree, h = x_radius)
    
    return(outlier_groups)
    
}


# Helper function to generate vector for grid creation
Grid_Modify <- function(x, layer1){
    orig_layer <- seq(x) - 1
    if(layer1 > 0){
        new_layer <- c(sort(-seq(layer1)),
                       orig_layer,
                       max(orig_layer) + sort(seq(layer1)))
    } else {
        new_layer <- orig_layer
    }
    return(new_layer)
}


# Helper function to calculate cell centers
Calculate_Cell_Centers <- function(cell, original_cell_nums, adjusted_cell_sizes, y_radius, y_min, y_max){
    
    cell_center <- cell
    
    for(i in 1:length(original_cell_nums)) {
        
        if(cell[[i]] < 0){
            
            cell_center[[i]] <- y_min[[i]] + cell[[i]] * 2 * y_radius + y_radius
            
        } else if (cell[[i]] >= original_cell_nums[[i]]) {
            
            cell_center[[i]] <- y_max[[i]] + (cell[[i]] - original_cell_nums[[i]]) * 2 * y_radius + y_radius
            
        } else {
            
            cell_center[[i]] <- y_min[[i]] + cell[[i]] * adjusted_cell_sizes[[i]] + adjusted_cell_sizes[[i]] / 2
            
        }
        
    }
    
    return(cell_center)
}