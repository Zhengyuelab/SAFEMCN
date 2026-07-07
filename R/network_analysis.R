#' Analyze Network Topology Parameters with Rarefaction
#'
#' Calculate network topology parameters from OTU tables with customizable
#' correlation thresholds, parallel processing options, and visualization
#' capabilities including topology parameter plots, exponential formula fitting,
#' prediction of future sample sizes, and AR1 analysis.
#'
#' @param work_dir Working directory path
#' @param otu_file OTU table file name
#' @param r_threshold Correlation coefficient threshold (default: 0.6)
#' @param p_threshold P-value threshold (default: 0.05)
#' @param cor_method Correlation method: "spearman" or "pearson" (default: "spearman")
#' @param start_size Starting sample size (default: 5)
#' @param end_size Ending sample size (default: NULL, uses total samples)
#' @param step_size Step size for sample sizes (default: 1)
#' @param replicates Number of replicates (default: 50)
#' @param use_parallel Use parallel processing (default: FALSE)
#' @param output_prefix Output file prefix (default: "network_parameters")
#' @param plot_topology Whether to plot topology parameters vs sample size (default: FALSE)
#' @param plot_start_size Starting sample size for topology plots (default: 11)
#' @param fit_formula Whether to perform exponential formula fitting (default: FALSE)
#' @param fit_start_size Starting sample size for formula fitting (default: 11)
#' @param predict_end_size End sample size for prediction (default: NULL, no prediction)
#' @param run_ar1 Whether to run AR1 analysis (default: FALSE)
#' @param ar1_input_file Input file name for AR1 analysis (default: NULL, auto-select)
#' @param ar1_windows Sliding window sizes for AR1 analysis (default: c(15))
#' @param ar1_start_size Starting sample size for AR1 analysis. If NULL, plot_start_size is used.
#' @return A list containing raw results, summary statistics, fit parameters, and predictions
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
                                     fit_start_size = 11,
                                     predict_end_size = NULL,
                                     run_ar1 = FALSE,
                                     ar1_input_file = NULL,
                                     ar1_windows = c(15)
                                     ar1_start_size = NULL) {

  setwd(work_dir)

  otu_table <- read.csv(otu_file, header = TRUE, row.names = 1, check.names = FALSE)
  otu_table <- t(otu_table)

  total_samples <- nrow(otu_table)

  if (is.null(end_size)) {
    end_size <- total_samples
  }

  sample_sizes <- seq(start_size, end_size, by = step_size)

  set.seed(123)

  # =========================================================
  # Internal: compute network parameters for one subsample
  # =========================================================
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

net <- igraph::graph_from_adjacency_matrix(
  adjacency_matrix,
  mode = "undirected",
  diag = FALSE
)

node_num <- igraph::vcount(net)
edge_num <- igraph::ecount(net)
average_degree <- mean(igraph::degree(net))
clustering_coefficient <- igraph::transitivity(net, type = "average")
fc <- igraph::cluster_fast_greedy(net)
modularity_value <- igraph::modularity(fc)
network_density <- igraph::edge_density(net)
average_path_length <- igraph::mean_distance(net)
network_diameter <- igraph::diameter(net)
    

return(data.frame(
  sample_size = sample_size,
  node_num = node_num,
  edge_num = edge_num,
  average_degree = average_degree,
  clustering_coefficient = clustering_coefficient,
  modularity = modularity_value,
  network_density = network_density,
  average_path_length = average_path_length,
  network_diameter = network_diameter
))
  }

  # =========================================================
  # Step 1: compute raw topology results (parallel or serial)
  # =========================================================
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

  # =========================================================
  # Step 2: compute summary (mean / sd for each sample size)
  # =========================================================
  params_names <- c("node_num", "edge_num", "average_degree", "clustering_coefficient",
                    "modularity", "network_density", "average_path_length", "network_diameter")

  summary_df <- data.frame()

  for (n in sample_sizes) {
    sub_df <- final_df[final_df$sample_size == n, ]

    if (nrow(sub_df) == 0) next

    row_data <- data.frame(sample_size = n)

    for (param in params_names) {
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

  # =========================================================
  # Step 3 (optional): plot topology parameters vs sample size
  # =========================================================
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

  # =========================================================
  # Step 4 (optional): exp3P fitting with nls()
  # =========================================================
  fit_summary <- NULL

  if (fit_formula) {

# Internal function: fit y = a * exp(b * x) + c using nls()
fit_exp3P <- function(x, y) {

  y_range <- diff(range(y))
  if (y_range == 0) {
    warning("Fit failed: y has no variation.")
    return(NULL)
  }

  # Determine whether the trajectory is increasing or decreasing
  trend_cor <- suppressWarnings(cor(x, y, method = "spearman"))

  if (is.na(trend_cor)) {
    warning("Fit failed: cannot determine the trend of y.")
    return(NULL)
  }

  if (trend_cor >= 0) {
    # Increasing saturating curve
    # Equivalent to y = c - |a| * exp(-|b| * x)
    # In the original formula, this means a < 0 and b < 0
    c_init <- max(y) + 0.01 * y_range
    a_init <- min(y) - c_init
    b_init <- -0.05

    lower_bounds <- c(a = -Inf, b = -Inf, c = max(y))
    upper_bounds <- c(a = 0,    b = 0,    c = Inf)

  } else {
    # Decreasing saturating curve
    # Equivalent to y = c + |a| * exp(-|b| * x)
    # In the original formula, this means a > 0 and b < 0
    c_init <- min(y) - 0.01 * y_range
    if (c_init < 0 && all(y > 0)) c_init <- min(y) * 0.9

    a_init <- max(y) - c_init
    b_init <- -0.05

    lower_bounds <- c(a = 0, b = -Inf, c = -Inf)
    upper_bounds <- c(a = Inf, b = 0, c = min(y))
  }

  tryCatch({
    fit <- nls(
      y ~ a * exp(b * x) + c,
      start = list(a = a_init, b = b_init, c = c_init),
      algorithm = "port",
      lower = lower_bounds,
      upper = upper_bounds,
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    )

    a_coef <- coef(fit)["a"]
    b_coef <- coef(fit)["b"]
    c_coef <- coef(fit)["c"]

    y_pred <- predict(fit)
    RSS <- sum((y - y_pred)^2)
    TSS <- sum((y - mean(y))^2)
    R2 <- 1 - RSS / TSS

    formula_str <- sprintf(
      "y = %.4f * exp(%.4f * x) + %.4f   (R^2 = %.4f)",
      a_coef, b_coef, c_coef, R2
    )

    return(list(
      fit = fit,
      params = coef(fit),
      formula = formula_str,
      R2 = R2
    ))

  }, error = function(e) {
    warning("Fit failed: ", e$message)
    return(NULL)
  })
}
    topo_params <- colnames(summary_df)[grepl("_mean$", colnames(summary_df))]
    topo_params <- gsub("_mean$", "", topo_params)

    fit_summary <- data.frame(
      parameter = character(),
      formula = character(),
      a = numeric(),
      b = numeric(),
      c = numeric(),
      R2 = numeric(),
      stringsAsFactors = FALSE
    )

    for (param in topo_params) {

      dat <- data.frame(
        sample_size = summary_df$sample_size,
        param_mean = summary_df[[paste0(param, "_mean")]]
      )
      colnames(dat)[2] <- "mean_value"

      dat_filtered <- dat[dat$sample_size >= fit_start_size & !is.na(dat$mean_value), ]

      if (nrow(dat_filtered) < 4) {
        warning(paste("Parameter", param, "has insufficient data points (", nrow(dat_filtered), "), skipping fit"))
        next
      }

      x_vals <- dat_filtered$sample_size
      y_vals <- dat_filtered$mean_value

      fit_result <- fit_exp3P(x_vals, y_vals)

      if (is.null(fit_result)) {
        warning(paste("Fit failed for parameter", param, ", skipping plot"))
        next
      }

      cat("Parameter:", param, "\n")
      cat("Fit formula:", fit_result$formula, "\n")
      cat("Fit parameters:\n")
      print(fit_result$params)
      cat("\n")

      fit_summary <- rbind(
        fit_summary,
        data.frame(
          parameter = param,
          formula = fit_result$formula,
          a = fit_result$params[["a"]],
          b = fit_result$params[["b"]],
          c = fit_result$params[["c"]],
          R2 = fit_result$R2,
          stringsAsFactors = FALSE
        )
      )

      # Generate prediction curve for the plot
      x_pred <- seq(min(x_vals), max(x_vals), length.out = 100)
      y_pred <- predict(fit_result$fit, newdata = list(x = x_pred))
      pred_df_plot <- data.frame(x = x_pred, y = y_pred)

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(data = dat_filtered,
                            ggplot2::aes(x = sample_size, y = mean_value),
                            size = 2, color = "#2E86AB") +
        ggplot2::geom_line(data = pred_df_plot,
                           ggplot2::aes(x = x, y = y),
                           color = "darkgreen", linewidth = 0.8) +
        ggplot2::labs(x = "Sample Size", y = param,
                      title = paste(param, "vs Sample Size"),
                      subtitle = paste0("Exponential fit (exp3P) from sample size = ", fit_start_size)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1))

      # Add formula annotation
      x_pos <- max(x_vals) - 0.1 * diff(range(x_vals))
      y_pos <- max(y_vals) - 0.1 * diff(range(y_vals))
      p <- p + ggplot2::annotate("text", x = x_pos, y = y_pos,
                                  label = fit_result$formula,
                                  hjust = 1, vjust = 1, size = 3, color = "black")

      print(p)

      ggplot2::ggsave(paste0("trendline_", param, "_from_", fit_start_size, ".png"),
                      p, width = 8, height = 6, dpi = 300)
    }

    # Save all fit parameters
    write.csv(fit_summary, "exp3P_fit_parameters.csv", row.names = FALSE)
    cat("All fit parameters saved to exp3P_fit_parameters.csv\n")
  }

  # =========================================================
  # Step 5 (optional): predict future sample sizes
  # =========================================================
  pred_df <- NULL

  if (!is.null(predict_end_size) && !is.null(fit_summary) && nrow(fit_summary) > 0) {

    new_sample_sizes <- (end_size + 1):predict_end_size

    if (length(new_sample_sizes) > 0) {
      pred_df <- data.frame(sample_size = new_sample_sizes)

      for (i in 1:nrow(fit_summary)) {
        param <- fit_summary$parameter[i]
        a <- fit_summary$a[i]
        b <- fit_summary$b[i]
        c_val <- fit_summary$c[i]
        pred_values <- a * exp(b * new_sample_sizes) + c_val
        pred_df[[param]] <- pred_values
      }

      pred_file <- paste0("predicted_topology_params_", end_size + 1, "_to_", predict_end_size, ".csv")
      write.csv(pred_df, pred_file, row.names = FALSE)
      cat("Predicted values saved to", pred_file, "\n")
    }
  }

  # =========================================================
  # Step 6: generate combined and original CSVs per parameter
  # =========================================================
  summary_file <- paste0(output_prefix, "_summary.csv")
  raw_df <- read.csv(summary_file, check.names = FALSE)

  mean_cols <- grep("_mean$", names(raw_df), value = TRUE)

  # --- Original CSVs (observed only) ---
  for (col in mean_cols) {
    param_name <- gsub("_mean$", "", col)
    out_data <- raw_df[, c("sample_size", col)]
    names(out_data) <- c("year", "area")
    out_file <- paste0(param_name, "_original.csv")
    write.csv(out_data, out_file, row.names = FALSE)
    cat("Generated file:", out_file, "\n")
  }
  cat("All original parameter tables generated.\n")

  # --- Combined CSVs (observed + predicted) ---
  if (!is.null(pred_df) && nrow(pred_df) > 0) {

    param_clean <- gsub("_mean$", "", mean_cols)

    # Build observed long format
    raw_long <- data.frame()
    for (col in mean_cols) {
      param_name <- gsub("_mean$", "", col)
      tmp <- data.frame(
        sample_size = raw_df$sample_size,
        parameter = param_name,
        area = raw_df[[col]],
        source = "observed",
        stringsAsFactors = FALSE
      )
      raw_long <- rbind(raw_long, tmp[!is.na(tmp$area), ])
    }

    # Build predicted long format
    pred_long <- data.frame()
    pred_param_cols <- setdiff(names(pred_df), "sample_size")
    for (param_name in pred_param_cols) {
      tmp <- data.frame(
        sample_size = pred_df$sample_size,
        parameter = param_name,
        area = pred_df[[param_name]],
        source = "predicted",
        stringsAsFactors = FALSE
      )
      pred_long <- rbind(pred_long, tmp)
    }

    combined_long <- rbind(raw_long, pred_long)
    combined_long <- combined_long[order(combined_long$parameter, combined_long$sample_size), ]

    for (param in unique(combined_long$parameter)) {
      param_data <- combined_long[combined_long$parameter == param, c("sample_size", "area")]
      names(param_data) <- c("year", "area")
      param_data <- param_data[order(param_data$year), ]
      out_file <- paste0(param, "_combined.csv")
      write.csv(param_data, out_file, row.names = FALSE)
      cat("Generated file:", out_file, "\n")
    }
    cat("All combined parameter tables generated.\n")
  }

  # =========================================================
  # Step 7 (optional): AR1 analysis (pure R implementation)
  # =========================================================
  if (run_ar1) {

    # Auto-select input file if not specified
    if (is.null(ar1_input_file)) {
      ar1_input_file <- "node_num_combined.csv"
      if (!file.exists(ar1_input_file)) {
        ar1_input_file <- "node_num_original.csv"
      }
    }

    cat("\n", paste(rep("=", 40), collapse = ""), "\n")
    cat("AR1 Analysis using:", ar1_input_file, "\n")
    cat(paste(rep("=", 40), collapse = ""), "\n")

    if (!file.exists(ar1_input_file)) {
      warning(paste("AR1 input file not found:", ar1_input_file))
    } else {

# --- AR1 helper functions ---
pop_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  sqrt(mean((x - mean(x))^2))
}

AR_func <- function(a1, a2) {
  ok <- is.finite(a1) & is.finite(a2)
  a1 <- a1[ok]
  a2 <- a2[ok]

  if (length(a1) == 0 || length(a2) == 0) return(NA_real_)

  std_a1 <- pop_sd(a1)
  std_a2 <- pop_sd(a2)

  if (is.na(std_a1) || is.na(std_a2) || std_a1 == 0 || std_a2 == 0) {
    return(0)
  }

  mean(((a1 - mean(a1)) / std_a1) * ((a2 - mean(a2)) / std_a2))
}

CalAR1 <- function(var, wlen) {
  n <- length(var)

  if (n <= wlen) {
    return(numeric(0))
  }

  AR1_result <- numeric(n - wlen)

  for (i in seq_len(n - wlen)) {
    a1 <- var[i:(i + wlen - 1)]
    a2 <- var[(i + 1):(i + wlen)]
    AR1_result[i] <- AR_func(a1, a2)
  }

  return(AR1_result)
}

CalVarn <- function(var, wlen) {
  n <- length(var)

  if (n < wlen) {
    return(numeric(0))
  }

  Varn_result <- numeric(n - wlen + 1)

  for (i in seq_len(n - wlen + 1)) {
    a1 <- var[i:(i + wlen - 1)]
    Varn_result[i] <- var(a1)
  }

  return(Varn_result)
}

      # --- Plot settings ---
      color1_line <- "#2c4ca0"
      color2_line <- "#aa4843"
      color1_fill <- "#2c4ca0"
      color2_fill <- "pink"

# Read data
df <- read.csv(ar1_input_file, check.names = FALSE)

if (!("year" %in% names(df))) {
  names(df) <- c("year", "value")[1:ncol(df)]
}

# Make AR1 start consistent with topology plots by default
if (is.null(ar1_start_size)) {
  ar1_start_size <- plot_start_size
}

# Filter AR1 input data
df <- df[df$year >= ar1_start_size, , drop = FALSE]

if (nrow(df) <= max(ar1_windows)) {
  stop(
    "Not enough data points for AR1 analysis after filtering by ar1_start_size. ",
    "Please reduce ar1_windows or lower ar1_start_size."
  )
}

time1 <- as.integer(df$year)
val11 <- as.numeric(df[[2]])

      # Detrend
      linear_detrend <- function(x) {
  idx <- seq_along(x)
  ok <- is.finite(x)

  if (sum(ok) < 2) {
    return(x - mean(x, na.rm = TRUE))
  }

  fit <- lm(x[ok] ~ idx[ok])
  out <- rep(NA_real_, length(x))
  out[ok] <- residuals(fit)
  out
}

val11de <- linear_detrend(val11)

      # Store all AR1 minimum results
      all_min_results <- data.frame(
        file = character(),
        window = integer(),
        Nmin = integer(),
        min_AR1 = numeric(),
        stringsAsFactors = FALSE
      )

      # Create PDF
      pdf_filename <- paste0(gsub("\\.csv$", "", ar1_input_file), "_results.pdf")
      grDevices::pdf(pdf_filename, width = 12, height = 4)

      for (w1 in ar1_windows) {
        nmin_year <- NA
        min_ar1_val <- NA
        
        # Detrended data for AR1
        val1 <- val11de

        AR1Val <- CalAR1(val1, w1)
        VarVal <- CalVarn(val1, w1)

# Time axis for AR1: center position of each sliding window
half_w <- w1 %/% 2

wt_ar1 <- time1[seq(
  from = half_w + 1,
  length.out = length(AR1Val)
)]

        # Save AR1 values to CSV
        output_ar1 <- data.frame(year = wt_ar1, AR1 = AR1Val)
        output_filename <- paste0(gsub("\\.csv$", "", ar1_input_file), "_window", w1, "_AR1.csv")
        write.csv(output_ar1, output_filename, row.names = FALSE)
        cat("Saved AR1 values to", output_filename, "\n")

# Find AR1 minimum after the window threshold
if (length(AR1Val) > 0) {
  
  # Only consider AR1 points whose sample-size position is >= window size
  valid_idx <- which(wt_ar1 >= w1 & !is.na(AR1Val))
  
  if (length(valid_idx) > 0) {
    min_ar1_idx <- valid_idx[which.min(AR1Val[valid_idx])]
    nmin_year <- wt_ar1[min_ar1_idx]
    min_ar1_val <- AR1Val[min_ar1_idx]
    
    all_min_results <- rbind(all_min_results, data.frame(
      file = ar1_input_file,
      window = w1,
      Nmin = nmin_year,
      min_AR1 = min_ar1_val,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf(
      "AR1 minimum after window threshold at year %d (value=%.4f, window=%d)\n",
      nmin_year, min_ar1_val, w1
    ))
  } else {
    warning(sprintf(
      "No valid AR1 points found after window threshold for window = %d",
      w1
    ))
  }
}

        # Plot: Original data + AR1
        par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

        # Left panel: Original data
        plot(time1, val11, type = "l", lwd = 1.5, col = color1_line,
             xlab = "Sample Size", ylab = "Original Value",
             main = paste0(gsub("\\.csv$", "", ar1_input_file), " (Window=", w1, ")"))

        # Right panel: AR1 coefficient
        plot(wt_ar1, AR1Val, pch = 16, col = color1_line, cex = 0.8,
             xlab = "Sample Size", ylab = "AR1 Coefficient",
             main = paste0("AR1 (Window=", w1, ")"))
        # Add vertical line for Nmin selected after the window threshold
        if (exists("nmin_year") && !is.na(nmin_year)) {
          abline(v = nmin_year, col = color2_line, lwd = 1.5, lty = 3)
          points(nmin_year, min_ar1_val, pch = 17, col = color2_line, cex = 1.2)
        }
        if (length(AR1Val) > 1) {
          # Add trend line for AR1
          ar1_lm <- lm(AR1Val ~ wt_ar1)
          abline(ar1_lm, col = color2_line, lwd = 1.5, lty = 2)
        }
      }

      grDevices::dev.off()
      cat("Created PDF:", pdf_filename, "\n")

      # Save AR1 minimum summary
      if (nrow(all_min_results) > 0) {
        write.csv(all_min_results, "AR1_min.csv", row.names = FALSE)
        cat("Saved AR1 minimum years summary to AR1_min.csv\n")
      }
    }

    cat("\nAR1 processing complete.\n")
  }

  # =========================================================
  # Return results
  # =========================================================
  result <- list(
    raw_results = final_df,
    summary = summary_df
  )

  if (!is.null(fit_summary)) {
    result$fit_parameters <- fit_summary
  }

  if (!is.null(pred_df)) {
    result$predictions <- pred_df
  }

  return(result)
}
