# *SAFEMCN* <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/SAFEMCN-logo.png" width="130" align="right" />
**Introduction** : Microbial co-occurrence networks are now widely used to infer putative ecological associations, identify keystone taxa, and compare community organization. However, the reliability and comparability of such networks are often limited by an underappreciated methodological issue, namely the dependence of inferred network topology on sample size. Because sample sizes differ substantially among ecological studies, networks constructed from insufficient or uneven sampling may yield unstable topological patterns and reduce the reproducibility of ecological interpretations.  

To address this issue, we developed *SAFEMCN*, an open-source R package that provides a sample-size-aware framework for evaluating microbial co-occurrence networks. *SAFEMCN* integrates rarefaction-based repeated network reconstruction with lag-1 autocorrelation analysis to identify the minimum sample size required for stable network inference. In addition, the package incorporates exponential curve fitting to predict this threshold from limited pilot data, allowing researchers to evaluate sampling adequacy before or during experimental design. 

***SAFEMCN*** : *SAFEMCN* (A Sample-Size-Aware Framework for Evaluating Microbial Co-occurrence Networks) is dependent on R ≥ 3.5.0, the other dependencies in the R environment are igraph, Hmisc, parallel, ggplot2, grDevices, graphics, rlangstats and utils.

# **Overview**  
The *SAFEMCN* package begins with an ASV/OTU table derived from raw amplicon sequencing data (Steps 1-2).   

Microbial co-occurrence networks are then constructed using user-defined correlation methods and |R| and p-value thresholds (Step 3).   

Next, network topological trajectories are generated along a rarefaction curve determined by the start sample size, end sample size, step size, and replicates (Step 4).   

For N<sub>min</sub> (Step 5), the package either identifies N<sub>min</sub> by applying lag-1 autoregressive (AR1) coefficient analysis to the observed topology (Step 5a), or predicts Nmin using a three-parameter exponential fitting model ($y = ae^{bx} + c$) and AR1 analysis on the fitted topology when empirical sample sizes are limited (Step 5b).   

Finally, the package outputs the N<sub>min</sub> (Step 6).

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/workflow-SAFEMCN.png" 
       width="95%" 
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
* **otu_file:** Input OTU/ASV table in CSV format (rows: features, cols: samples).  
* **r_threshold:** The correlation coefficient threshold (|R|).   
* **p_threshold:** The significance level (*p*,filters non-significant edges). 
* **cor_method:** Correlation method ("spearman" or "pearson").  
* **start_size & end_size:** Defines the range of the sample-size gradient (e.g., from 5 to 58 samples).  
* **step_size:** Increment of gradient (1 = single sample steps).  
* **replicates:** The number of repeated sampling for each sample size (e.g., 50 repeats).
* **use_parallel:** If TRUE, enable multi-core parallel processing.
* **plot_topology:** If TRUE, generates plots showing the trend of topological parameters across the sample-size gradient.  
* **plot_start_size:** Start sample size for plotting.  
* **fit_formula:** If TRUE, apply exponential fitting to parameters.  
* **fit_start_size:**  Start sample size for fitting.  
* **predict_end_size:** Target sample size for extrapolation.  
* **run_ar1:** if TRUE, performs lag-1 Autoregressive (AR1) analysis to calculate the AR1 coefficient at each sample size.  
* **ar1_input_file:** The input file can be either node_num_original.csv or node_num_combined.csv which is the result of combining the fitting data.
* **ar1_windows**: Sliding window size(e.g.,15) for computing rolling AR1 coefficients and variances.

### **2. Case Studies**
#### **2.1 Determination of N<sub>min</sub>**
In this part, the long-term time-series dataset from Station ALOHA is used as an example data set to determine the required N<sub>min</sub> for this dataset.
```r
library(SAFEMCN)

results <- analyze_network_topology(
  otu_file = "otu.csv",
  r_threshold = 0.6,
  p_threshold = 0.05,
  cor_method = "spearman",
  start_size = 5,
  end_size = 58,
  step_size = 1,
  replicates = 50,
  use_parallel = TRUE,
  output_prefix = "network_parameters",
  plot_topology = TRUE,
  plot_start_size = 11,
  fit_formula = FALSE,
  fit_start_size = 11,
  predict_end_size = 58,
  run_ar1 = TRUE,
  ar1_input_file = "node_num_combined.csv",
  ar1_windows = c(15)
)
```

<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/ALOHA-test.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>


#### 📌 **Result Interpretation:**
The *SAFEMCN* package visualizes the variation of network topology parameters with sample size as a line graph, while also outputting the corresponding topological parameter values for each sample size (e.g., node_num_original.csv, edge_num_original.csv).   

Furthermore, the AR1 module of the *SAFEMCN*  package takes these output files (e.g., node_num_original.csv) as input to identify the N<sub>min</sub> threshold: the lowest AR1 value (marked by a rapid decrease followed by a rapid increase) is determined as N<sub>min</sub>, indicating the emergence of a topological plateau.   

Finally, the AR1 module outputs a CSV file containing the AR1 coefficient for each sample size (e.g., node_num_original_window15_AR1.csv) and also generates an AR1_min.csv file, which records the determined N<sub>min</sub> value and its corresponding AR1 value.  

#### **2.2 Prediction of N<sub>min</sub>**
This part uses the ALOHA dataset as the sample data and employs the prediction module of *SAFEMCN* to predict Nmin.
```r
library(SAFEMCN)

results <- analyze_network_topology(
  otu_file = "otu.csv",
  r_threshold = 0.6,
  p_threshold = 0.05,
  cor_method = "spearman",
  start_size = 5,
  end_size = 30,
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
  ar1_windows = c(15)
)

```
<p align="center">
  <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/ALOHA-predicted.png" 
       width="65%" 
       alt="Variations in network topological parameters and the lag-1 autocorrelation (AR1) coefficient" />
</p>

#### 📌 **Result Interpretation:**
First, the *SAFEMCN* package outputs a line-point plot of network topology parameters as a function of sample size under finite sample sizes. When `fit_formula = TRUE`, the package fits the existing trend using a three-parameter exponential formula to obtain the fitted equation.  

After setting `predict_end_size = 58`, the package calculates the corresponding network topology parameter values for sample sizes ranging from 31 to 58 based on the fitted equation, and integrates the observed values with the fitted values for output into files such as node_num_combined.csv and edge_num_combined.csv. These files are then used as inputs to the AR1 model for the determination of N<sub>min</sub>.  

Similarly, this fitting part will also output the node_num_combined_window15_AR1.csv and the corresponding AR1_min.csv files.
