# Main Engine -----------------------------------------------------------------------
LION_tSNE <- function(X_new,
                      X_train,
                      Y_train,
                      p = 20,
                      r_x_percentile = 0.99,
                      r_y_percentile = 0.99,
                      r_y_close_percentile = 0.1,
                      subsample = 1){
    
    # X_new <- data_X_test[1:300, ]
    # Y_new <- data_Y_test[1:300, ]
    # X_train <- data_X_train
    # Y_train <- data_Y_train
    # 
    # # Test parameters
    # p = 20
    # r_x_percentile = 0.5
    # r_y_percentile = 0.5
    # r_y_close_percentile = 0.1
    # r_y_coeff = 1
    # subsample = 0.1
    
    # Initialize Y_output matrix
    Y_output <- data.table::data.table(
        matrix(as.numeric(NA), 
               ncol = ncol(Y_train),
               nrow = nrow(X_new)
        )
    )
    
    # Step 0: Subsample training data if required (To calculate radius only!)
    subsample_rows <- sample(c(1:nrow(X_train)),
                             round(nrow(X_train) * subsample))
    
    X_train_subsample <- X_train[subsample_rows, ]
    Y_train_subsample <- Y_train[subsample_rows, ]
    
    cat('\nCalculating global radius')
    
    # Get global radius
    x_radius <- Radius_Calculator(X_train_subsample, r_x_percentile)
    y_radius <- Radius_Calculator(Y_train_subsample, r_y_percentile)
    
    cat('\nCalculating pair-wise distances')
    
    # Get neighbours for each x_new with x_radius
    X_train_X_new_nn <- RANN::nn2(X_train, X_new, k = nrow(X_train), searchtype = 'radius', radius = x_radius)
    X_train_X_new_nn_id <- data.table::data.table(t(X_train_X_new_nn$nn.idx))
    dist_new <- data.table::data.table(t(X_train_X_new_nn$nn.dists))
    neighbours_list <- lapply(X_train_X_new_nn_id, function(x) x[x > 0])
    
    # 1. Get exact match
    dist_min <- as.numeric(apply(dist_new, 2, min))
    dist_min_index <- as.numeric(apply(dist_new, 2, which.min))
    
    exact_match_index_new <- which(dist_min == 0 & is.na(Y_output[, 1]))
    exact_match_index_old <- dist_min_index[exact_match_index_new]
    
    if(length(exact_match_index_new) > 0) {
        
        cat('\nExact Match:', 'Processing', length(exact_match_index_new), 'observations')
        
        Y_output[exact_match_index_new, ] <- Y_train[exact_match_index_old, ]
        
    }
    
    # 2. Get Local IDW interpolation    
    num_of_neighbours_in_x_radius <- unlist(lapply(neighbours_list, length))
    IDW_index_new <- which(num_of_neighbours_in_x_radius > 1 & is.na(Y_output[, 1]))
    
    if(length(IDW_index_new) > 0){
        
        cat('\nIDW Interpolaton:', 'Processing', length(IDW_index_new), 'observations')
        
        Y_new_interpolated <- data.table::rbindlist(
            lapply(as.list(IDW_index_new),
                   function(x) Local_IDW_interpolation(X_new[x, ], X_train, Y_train, neighbours_list[[x]], p)
            )
        )
        
        Y_output[IDW_index_new, ] <- Y_new_interpolated
        
    }
    
    # 3. Get single neighbour
    single_neighbour_index <- which(num_of_neighbours_in_x_radius == 1  & is.na(Y_output[, 1]))
    
    if(length(single_neighbour_index) > 0){
        
        # Get the X_train indexes corresponding to the single neighbours
        X_train_index <- unlist(neighbours_list[single_neighbour_index])
        
        # Get distance of single neighbour to its nearest neighbour
        dist1 <- as.matrix(pdist::pdist(X = X_train, Y = X_train[X_train_index, ]))
        min_dist1 <- apply(dist1, 2, function(x) min(x[x > 0]))
        
        # Get xi which is not outlier
        xi_not_outlier_index <- single_neighbour_index[which(min_dist1 <= x_radius)]
        
        # Get xi which is outlier
        xi_outlier_index <- single_neighbour_index[which(min_dist1 > x_radius)]
        
    }
    
    
    # If xi is an outlier, perform single neighbour placement method else perform outlier placement method
    if(length(xi_outlier_index) > 0){
        
        cat('\nSingle Neighbour Placement:', 'Processing', length(xi_outlier_index), 'observations')
        
        # Precompute r_y_close
        y_radius_close <- Radius_Calculator(Y_train_subsample, r_y_close_percentile)
        
        # Perform single neighbour placement
        Y_new_single_neighbour_placement <- data.table::rbindlist(
            lapply(as.list(xi_outlier_index),
                   function(x) Single_neighbour_placement(Y_train, y_radius_close, neighbours_list[[x]])
            )
        )
        
        Y_output[xi_outlier_index, ] <- Y_new_single_neighbour_placement
        
    }
    
    
    # Compile outlier placement index
    outlier_placement_index_temp <- which(num_of_neighbours_in_x_radius == 0  & is.na(Y_output[, 1]))
    outlier_placement_index <- c(outlier_placement_index_temp, xi_not_outlier_index)
    
    if(length(outlier_placement_index) > 0){
        
        cat('\nOutlier Placement:', 'Processing', length(outlier_placement_index), 'observations')
        
        # Precompute r_y_close
        y_radius_close <- Radius_Calculator(Y_train_subsample, r_y_close_percentile)
        
        # Precompute outlier positions
        available_cells <- PreCompute_Outlier_Positions(Y_train, y_radius)
        
        # Get outliers
        outliers1 <- X_new[outlier_placement_index, ]
        
        # Get outlier groups
        outlier_groups <- Group_Outliers(outliers1, x_radius)
        
        # Create index-group lookup
        index_group_lookup <- data.table::data.table(outlier_index = outlier_placement_index,
                                                     group = outlier_groups)
        
        # Get centroid of outlier groups 
        # for groups with only 1 point, centroid will be itself
        outlier_centroid <- aggregate(outliers1, list(outlier_groups), mean)
        outlier_counts <- as.numeric(table(outlier_groups))
        
        # find nearest neighbour of centroid from training data
        get_knn_train <- RANN::nn2(X_train, outlier_centroid[, -1], k = 1)
        get_knn_train_index <- get_knn_train$nn.idx
        get_knn_train_Y <- Y_train[get_knn_train_index, ]
        
        # Place outlier representatives one-by-one
        current_outer_layers = 0
        current_available_cells = available_cells
        
        # Number of cells per dimension and cell sizes (Recaculated as it was not saved)
        y_min <- unlist(lapply(Y_train, min))
        y_max <- unlist(lapply(Y_train, max))
        
        original_cell_nums <- list()
        adjusted_cell_sizes <- list()
        
        for(i in 1:dim(Y_train)[2]){
            original_cell_nums_temp <- floor((y_max[i] - y_min[i]) / (2 * y_radius)) 
            adjusted_cell_sizes_temp <- (y_max[i] - y_min[i]) / original_cell_nums_temp
            
            original_cell_nums[[i]] <- original_cell_nums_temp
            adjusted_cell_sizes[[i]] <- adjusted_cell_sizes_temp
        }
        
        # Calculate cell centers
        cell_centers <- lapply(current_available_cells, 
                               function(x) Calculate_Cell_Centers(x, original_cell_nums, adjusted_cell_sizes, y_radius, y_min, y_max))
        
        # Convert cell centers and current available cells to data frame to avoid indexing problem
        cell_centers <- data.table::rbindlist(cell_centers)
        current_available_cells <- data.table::rbindlist(current_available_cells)
        
        # Get cell grid (df)
        cell_df <- expand.grid(lapply(original_cell_nums, function(x) Grid_Modify(x, current_outer_layers)))
        
        for(i in sort(unique(outlier_groups))){
            
            # If out of cells, request new ones [new layer]
            if(nrow(current_available_cells) == 0){
                
                current_outer_layers <- current_outer_layers + 1
                
                # Get updated grid
                cell_df_update <- expand.grid(lapply(original_cell_nums, function(x) Grid_Modify(x, current_outer_layers)))
                
                # Get anti-join
                data.table::setDT(cell_df_update)
                data.table::setDT(cell_df)
                
                cell_df_new <- fsetdiff(cell_df_update, cell_df, all = TRUE)
                
                # Update current available cells
                current_available_cells <- split(cell_df_new, seq(nrow(cell_df_new)))
                
                # Calculate cell centers
                cell_centers <- lapply(current_available_cells, 
                                       function(x) Calculate_Cell_Centers(x, original_cell_nums, adjusted_cell_sizes, y_radius, y_min, y_max))
                
                # Convert cell centers and current available cells to data frame to avoid indexing problem
                cell_centers <- data.table::rbindlist(cell_centers)
                current_available_cells <- data.table::rbindlist(current_available_cells)
                
                # Update grid (replace)
                cell_df <- data.table::data.table(cell_df_update)
                
            }
            
            # Get the Y (embedded value) of the nearest neighbour to outlier centroid / single
            outlier_y_train <- get_knn_train_Y[i, ]
            
            # Get nearest cell to the nearest neighbour above
            nearest_cell <- RANN::nn2(cell_centers, outlier_y_train, k = 1)
            
            # Get outlier index
            outlier_index <- as.numeric(subset(index_group_lookup, group == i, select = 'outlier_index')[, 1])
            
            # Place outlier cluster around nearest cell center
            if(outlier_counts[i] > 1){
                
                # Perform single neighbour placement
                Y_out_single_neighbour_placement <- data.table::rbindlist(
                    lapply(as.list(outlier_index),
                           function(x) Single_neighbour_placement(cell_centers[as.numeric(nearest_cell$nn.idx), ], y_radius_close, 1)
                    )
                )
                
                # overwrite Y_output 
                Y_output[outlier_index, ] <- Y_out_single_neighbour_placement
                
            } else {
                
                Y_output[outlier_index, ] <- cell_centers[as.numeric(nearest_cell$nn.idx), ]
                
            }
            
            # Remove used cell from cell_centers and current_available_cells
            cell_centers <- cell_centers[-as.numeric(nearest_cell$nn.idx), ]
            current_available_cells <- current_available_cells[-as.numeric(nearest_cell$nn.idx), ]
            
        }
        
    }
    
    # Create Outputs
    index <- list(as.numeric(exact_match_index_new), 
                  as.numeric(IDW_index_new),
                  as.numeric(xi_outlier_index),
                  as.numeric(outlier_placement_index))
    
    names(index) <- c('exact_match', 'IDW', 'single_neighbour', 'outlier')
    
    out <- list(index, Y_output)
    names(out) <- c('index', 'Y')
    
    cat('\nLION-tSNE Complete')
    
    return(out)
    
}