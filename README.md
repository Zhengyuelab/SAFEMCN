# *SAFEMCN* <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/logo.png" width="80" align="right" />

***SAFEMCN*** : *SAFEMCN* (A Sample-Size-Aware Framework for Evaluating Microbial Co-occurrence Networks) is dependent on R ≥ 3.5, the other dependencies in the R environment are igraph, Hmisc, parallel, ggplot2, grDevices, graphics, rlangstats and utils.


# **Overview**  
The *SAFEMCN* R package is a standardized framework designed to determine the minimum sample size (N<sub>min</sub>) required for robust microbial co-occurrence network inference, which systematically integrates three core functional modules: (1) Network construction through automated, parallelized subsampling across designated sample-size gradients; (2) Determination of the N<sub>min</sub> via lag-1 autoregressive (AR1) coefficient assessments; and (3) N<sub>min</sub> estimation via exponential model extrapolation, which provides a predictive, three-parameter exponential fitting strategy to prospectively forecast dataset requirements under constrained sample sizes.

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/Workflow-SAFEMCN.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>


# **Main Features**  
* Traces topological structure changes across a sample-size gradient using a rarefaction-curve approach.
* The lag-1 autocorrelation (AR1) coefficient is calculated to determine the minimum effective sample size (N<sub>min</sub>) for network structure stability.
* Extrapolates topological trajectories beyond the sample-size range via a three-parameter exponential model and uses lag-1 autocorrelation (AR1) coefficients on fitted sequences to predict N<sub>min</sub>.
* Supports parallel computation for large-scale datasets by distributing iterations across multiple CPU cores.


# **Install R/RStudio**  
If you do not already have R/RStudio installed, follow these steps:   
1. Install [R](https://www.r-project.org/)
2. Install [RStudio](https://posit.co/download/rstudio-desktop/)

Open RStudio -> Tools -> Global Options -> Packages, select the appropriate mirror in Primary CRAN repository.

# **Install *SAFEMCN*** 
Install *SAFEMCN*  package from CRAN.
```r
install.packages("SAFEMCN")
```
Install the latest development version from Github.
```r
#If devtools package is not installed, first install it.
install.packages("devtools")
devtools::install_github("Zhengyuelab/SAFEMCN")
```




# **Usage**  
### **1. Parameter Settings**
* work_dir: Specifies the working directory for reading input files and saving output results.  
* otu_file: Input OTU/ASV table (rows: features, cols: samples).  
* r_threshold: The correlation coefficient threshold (|R|).   
* p_threshold: The significance level(filters non-significant edges). 
* cor_method: Correlation algorithm ("spearman" or "pearson").  
* start_size & end_size: Defines the range of the sample-size gradient (e.g., from 5 to 58 samples).  
* step_size: Increment of gradient (1 = single sample steps).  
* replicates: The number of repeated sampling for each sample size (e.g., 50 repeats).
* use_parallel: Enable multi-core parallel processing.
* output_prefix: Prefix for exported result files. 
* plot_topology: Logical; if TRUE, generates plots showing the trend of topological parameters across the sample-size gradient.  
* plot_start_size: Start sample size for plotting.  
* fit_formula: Logical; if TRUE, Apply exponential fitting to parameters.  
* fit_start_size:  Start sample size for fitting.  
* predict_end_size: Target sample size for extrapolation.  
* run_ar1: Logical; if TRUE, performs lag-1 Autoregressive (AR1) analysis to calculate the AR1 coefficient at each sample size, serving as an indicator of network stability.  
* ar1_input_file: The input file for AR1 analysis, containing network topological parameters (e.g., node number, edge number, or density) calculated across the sample-size gradient.  
* ar1_windows: Defines the sliding window sizes (e.g.,15) for calculating rolling statistics like standard deviation or autocorrelation in the AR1 stability analysis.

### **2. Case Studies**
#### **2.1 Determination of N<sub>min</sub>**
Based on the long-term time-series dataset from Station ALOHA, a gradient of microbial co-occurrence networks was constructed by varying the sample size from 11 to 58 with a step size of 1. As illustrated in Figure 2, key network topological features, including network nodes, edges, density, diameter, average path length, and average degree, exhibited distinct dynamic patterns with expanding sample sizes, based on which the minimum sample size (N<sub>min</sub>) required for reliable network construction was successfully determined via lag-1 autoregressive (AR1) coefficients.
```r
library(SAFEMCN)

results <- analyze_network_topology(
  otu_file = "otu.csv",
  r_threshold = 0.6,
  p_threshold = 0.05,
  cor_method = "spearman",
  start_size = 5,
  end_size = 28,
  step_size = 1,
  replicates = 50,
  use_parallel = TRUE,
  output_prefix = "network_parameters",
  plot_topology = TRUE,
  plot_start_size = 11,
  fit_formula = TRUE,
  fit_start_size = 11,
  predict_end_size = 58,
  run_ar1 = TRUE,
  ar1_input_file = "node_num_combined.csv",
  ar1_windows = c(10, 15)
)
```

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/ALOHA-test.png" 
       width="65%" 
       alt="Variations in network topological parameters" />
</p>
<p align="center">
  <b>Figure 2. Variations in network topological parameters</b>
</p>

#### 📌 **Result Interpretation:**
With increasing sample sizes, node numbers expanded toward saturation, whereas edge counts dramatically decreased before leveling off. Across these trajectories, the network topology parameters stabilized in a clear sequential order: node first, followed by edge, density, average degree, and lastly, diameter and average path length, reflecting a gradient of sensitivity to sample size.  

The N<sub>min</sub> threshold was further identified by the AR1 module of the *SAFEMCN* package. The rapid decrease and then rapid increase in the AR1 value indicated the emergence of the topological plateau, providing a basis to determine the critical point corresponding to N<sub>min</sub>. To provide a robust estimate, the final Nmin was calculated as the average of the N<sub>min</sub> derived from all six topological parameters, resulting in a value of 31 ± 5. This threshold marked the earlier sampling depth at which inferred topology became reproducible across resampling replicates and could be regarded as meeting a minimum level of statistical validity.

#### **2.2 Prediction of N<sub>min</sub>**
To enable topology prediction when initial sampling was limited, the *SAFEMCN* package integrated a predictive module that coupled the AR1 analysis with curve fitting of network topology parameters. This functionality enables users to model the trajectories of network topological parameters as nonlinear functions varying with sample size using an exponential fitting model (). By fitting this model to observed data, the package extrapolates topological trajectories beyond the current sample-size range and computes AR1 coefficients on the fitted sequences to determine the predicted N<sub>min</sub>.


#### 📌 **Result Interpretation:**
By extrapolating node and edge trajectories via the fitted exponential model (Figs 3a–b), the AR1-based N<sub>min</sub> was successfully predicted under limited sampling. Crucially, the predicted N<sub>min</sub> tightly matched the empirical thresholds derived from the complete 58-sample dataset (node: 25 vs. 27 samples; edge: 29 vs. 28 samples). This minimal deviation validates the framework's reliability. Furthermore, networks constructed with sample sizes exceeding these predicted thresholds exhibited convergent and stable topological properties, such as network density and average path length (Figs 3c–d).  

Consequently, combining mechanistic extrapolation with the AR1 stability criterion provides a robust, practical solution for optimized sampling design when large datasets are unfeasible.
