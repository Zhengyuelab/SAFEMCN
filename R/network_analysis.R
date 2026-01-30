#' @export
analyze_network_topology <- function(work_dir,
                                     otu_file,
                                     r_threshold = 0.6,
                                     p_threshold = 0.05,
                                     cor_method = "spearman",
                                     start_size = 5,
                                     end_size = NULL,
                                     step_size = 1,
                                     replicates = 50,
                                     use_parallel = FALSE,
                                     output_prefix = "network_parameters",
                                     plot_topology = FALSE,
                                     plot_start_size = 11,
                                     fit_formula = FALSE,
                                     fit_start_size = 11) {
  
  setwd(work_dir)
  
  otu_table <- read.csv(otu_file, header = TRUE, row.names = 1, check.names = FALSE)
  otu_table <- t(otu_table)
  
  total_samples <- nrow(otu_table)
  
  if (is.null(end_size)) {
    end_size <- total_samples
  }
  
  sample_sizes <- seq(start_size, end_size, by = step_size)
  
  set.seed(123)
  
  get_network_params <- function(sample_size, otu_table, r_thresh, p_thresh, method) {
    samples <- sample(rownames(otu_table), sample_size)
    otu_sub <- otu_table[samples, ]
    
    otu_sub <- otu_sub[, colSums(otu_sub) > 0]
    
    if (ncol(otu_sub) < 2) {
      return(data.frame(
        sample_size = sample_size,
        node_num = NA,
        edge_num = NA,
        average_degree = NA,
        clustering_coefficient = NA,
        modularity = NA,
        network_density = NA,
        average_path_length = NA,
        network_diameter = NA
      ))
    }
    
    cor_matrix <- Hmisc::rcorr(as.matrix(otu_sub), type = method)
    R <- cor_matrix$r
    P <- cor_matrix$P
    
    adjacency_matrix <- ifelse(!is.na(R) & abs(R) > r_thresh & P < p_thresh, 1, 0)
    diag(adjacency_matrix) <- 0
    
    if (sum(adjacency_matrix, na.rm = TRUE) == 0) {
      return(data.frame(
        sample_size = sample_size,
        node_num = ncol(otu_sub),
        edge_num = 0,
        average_degree = 0,
        clustering_coefficient = 0,
        modularity = 0,
        network_density = 0,
        average_path_length = 0,
        network_diameter = 0
      ))
    }
    
    net <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", diag = FALSE)
    
    node_num <- igraph::vcount(net)
    edge_num <- igraph::ecount(net)
    average_degree <- mean(igraph::degree(net))
    clustering_coefficient <- igraph::transitivity(net, type = "average")
    fc <- igraph::cluster_fast_greedy(net)
    modularity_value <- igraph::modularity(fc)
    network_density <- igraph::edge_density(net)
    average_path_length <- igraph::mean_distance(net)
    network_diameter <- igraph::diameter(net)
    
    data.frame(
      sample_size = sample_size,
      node_num = node_num,
      edge_num = edge_num,
      average_degree = average_degree,
      clustering_coefficient = clustering_coefficient,
      modularity = modularity_value,
      network_density = network_density,
      average_path_length = average_path_length,
      network_diameter = network_diameter
    )
  }
  
  if (use_parallel) {
    
    calculate_parallel <- function(n, otu_table, replicates, r_thresh, p_thresh, method) {
      
      cl <- parallel::makeCluster(parallel::detectCores())
      
      parallel::clusterExport(cl, c("get_network_params", "otu_table", "n", "r_thresh", "p_thresh", "method"), envir = environment())
      
      parallel::clusterEvalQ(cl, {
        library(igraph)
        library(Hmisc)
      })
      
      results <- parallel::parLapply(cl, 1:replicates, function(rep) {
        result <- tryCatch({
          params <- get_network_params(n, otu_table, r_thresh, p_thresh, method)
          params$replicate <- rep
          params
        }, error = function(e) {
          data.frame(sample_size = n,
                     node_num = NA,
                     edge_num = NA,
                     average_degree = NA,
                     clustering_coefficient = NA,
                     modularity = NA,
                     network_density = NA,
                     average_path_length = NA,
                     network_diameter = NA,
                     replicate = rep)
        })
        return(result)
      })
      
      parallel::stopCluster(cl)
      
      do.call(rbind, results)
    }
    
    final_results <- list()
    for (n in sample_sizes) {
      cat("Processing sample size:", n, "\n")
      res <- calculate_parallel(n, otu_table, replicates, r_threshold, p_threshold, cor_method)
      final_results[[as.character(n)]] <- res
    }
    
    final_df <- do.call(rbind, final_results)
    
  } else {
    
    results_list <- list()
    for (n in sample_sizes) {
      cat("Processing sample size:", n, "\n")
      for (rep in 1:replicates) {
        params <- get_network_params(n, otu_table, r_threshold, p_threshold, cor_method)
        params$replicate <- rep
        results_list[[length(results_list) + 1]] <- params
      }
    }
    
    final_df <- do.call(rbind, results_list)
  }
  
  write.csv(final_df, paste0(output_prefix, "_raw_results.csv"), row.names = FALSE)
  
  params <- c("node_num", "edge_num", "average_degree", "clustering_coefficient", 
              "modularity", "network_density", "average_path_length", "network_diameter")
  
  summary_df <- data.frame()
  
  for (n in sample_sizes) {
    sub_df <- final_df[final_df$sample_size == n, ]
    
    if (nrow(sub_df) == 0) next
    
    row_data <- data.frame(sample_size = n)
    
    for (param in params) {
      values <- sub_df[[param]]
      values <- values[!is.na(values)]
      
      if (length(values) > 0) {
        mean_val <- mean(values)
        sd_val <- sd(values)
      } else {
        mean_val <- NA
        sd_val <- NA
      }
      
      row_data[[paste0(param, "_mean")]] <- mean_val
      row_data[[paste0(param, "_sd")]] <- sd_val
    }
    
    summary_df <- rbind(summary_df, row_data)
  }
  
  write.csv(summary_df, paste0(output_prefix, "_summary.csv"), row.names = FALSE)
  
  if (plot_topology) {
    
    summary_df_data <- summary_df
    
    topo_params <- colnames(summary_df_data)[grepl("_mean$", colnames(summary_df_data))]
    topo_params <- gsub("_mean$", "", topo_params)
    
    filtered_data <- summary_df_data[summary_df_data$sample_size >= plot_start_size, ]
    
    for (param in topo_params) {
      mean_col <- paste0(param, "_mean")
      sd_col <- paste0(param, "_sd")
      
      if (!(mean_col %in% colnames(filtered_data) && sd_col %in% colnames(filtered_data))) {
        next
      }
      
      p <- ggplot2::ggplot(filtered_data, ggplot2::aes(x = sample_size, y = .data[[mean_col]])) +
        ggplot2::geom_point(size = 2, color = "#2E86AB") +
        ggplot2::geom_line(linewidth = 0.8, color = "#2E86AB", alpha = 0.7) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[mean_col]] - .data[[sd_col]], 
                          ymax = .data[[mean_col]] + .data[[sd_col]]),
                      width = 0.8, color = "#2E86AB", alpha = 0.6) +
        ggplot2::labs(x = "Sample Size", y = param,
             title = paste(param, "vs Sample Size"),
             subtitle = paste("From sample size =", plot_start_size)) +
        ggplot2::scale_x_continuous(breaks = seq(min(filtered_data$sample_size), 
                                        max(filtered_data$sample_size), by = 5)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1))
      
      print(p)
      
      ggplot2::ggsave(paste0("topology_", param, "_from_", plot_start_size, ".png"), 
             p, width = 8, height = 6, dpi = 300)
    }
  }
  
  if (fit_formula) {
    
    topo_params <- colnames(summary_df)[grepl("_mean$", colnames(summary_df))]
    topo_params <- gsub("_mean$", "", topo_params)
    
    for (param in topo_params) {
      
      dat <- data.frame(
        sample_size = summary_df$sample_size,
        param_mean = summary_df[[paste0(param, "_mean")]]
      )
      
      colnames(dat)[2] <- "mean_value"
      
      dat_filtered <- dat[dat$sample_size >= fit_start_size, ]
      
      if (nrow(dat_filtered) < 3) {
        warning(paste("Parameter", param, "has insufficient data points, skipping fit"))
        next
      }
      
      x.pan <- dat_filtered$sample_size
      y.pan <- dat_filtered$mean_value
      
      tryCatch({
        
        fit_result <- ggtrendline::ggtrendline(x.pan, y.pan, model = "exp3P", 
                              CI.fill = "darkgreen",
                              CI.color = "darkgreen",
                              xlab = paste("Sample Size (>=", fit_start_size, ")"),
                              ylab = param) +
          ggplot2::ggtitle(paste(param, "vs Sample Size (", fit_start_size, "+ samples)")) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
          )
        
        print(fit_result)
        
        ggplot2::ggsave(paste0("trendline_", param, "_from_", fit_start_size, ".png"), 
               fit_result, width = 8, height = 6, dpi = 300)
        
        cat("Completed fit and plot for parameter:", param, "\n")
        
      }, error = function(e) {
        warning(paste("Fit failed for parameter", param, ":", e$message))
      })
    }
  }
  
  return(list(raw_results = final_df, summary = summary_df))
}