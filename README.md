# *SAFEMCN* <img src="https://github.com/Zhengyuelab/SAFEMCN/blob/main/figure/logo.png" width="80" align="right" />

***SAFEMCN*** : *SAFEMCN* (A Sample-Size-Aware Framework for Evaluating Microbial Co-occurrence Networks) is dependent on R ≥ 3.5, the other dependencies in the R environment are igraph, Hmisc, parallel, ggplot2, grDevices, graphics, rlangstats and utils.


# **Overview**  
The *SAFEMCN* R package is a standardized framework designed to determine the minimum sample size (N<sub>min</sub>) required for robust microbial co-occurrence network inference, which systematically integrates three core functional modules: (1) Variations in network topological properties across diverse sample sizes are characterized via a rarefaction-curve method; (2) Determination of the N<sub>min</sub> via lag-1 autoregressive (AR1) coefficient assessments; and (3) N<sub>min</sub> estimation via exponential model extrapolation, which provides a predictive, three-parameter exponential fitting strategy to prospectively forecast dataset requirements under constrained sample sizes.

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
* **otu_file:** Input OTU/ASV table in CSV format (rows: features, cols: samples).  
* **r_threshold:** The correlation coefficient threshold (|R|).   
* **p_threshold:** The significance level (*p*,filters non-significant edges). 
* **cor_method:** Correlation method ("spearman" or "pearson").  
* **start_size & end_size:** Defines the range of the sample-size gradient (e.g., from 5 to 58 samples).  
* **step_size:** Increment of gradient (1 = single sample steps).  
* **replicates:** The number of repeated sampling for each sample size (e.g., 50 repeats).
* **use_parallel:** Enable multi-core parallel processing.
* **plot_topology:** if TRUE, generates plots showing the trend of topological parameters across the sample-size gradient.  
* **plot_start_size:** Start sample size for plotting.  
* **fit_formula:** if TRUE, apply exponential fitting to parameters.  
* **fit_start_size:**  Start sample size for fitting.  
* **predict_end_size:** Target sample size for extrapolation.  
* **run_ar1:** if TRUE, performs lag-1 Autoregressive (AR1) analysis to calculate the AR1 coefficient at each sample size, serving as an indicator of network stability.  
* **ar1_input_file:** The input file for AR1 analysis, The input file can be either node_num_original.csv or node_num_combined.csv which is the result of combining the fitting data.
* **ar1_windows**: Defines the sliding window sizes (e.g.,15) for calculating rolling statistics like standard deviation or autocorrelation in the AR1 stability analysis.

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
The *SAFEMCN* package visualizes the variation of network topology parameters with sample size as a line graph, providing an intuitive representation of topological changes, while also outputting the corresponding topological parameter values for each sample size (e.g., node_num_original.csv, edge_num_original.csv) and saving them in a designated folder.   

Furthermore, the AR1 module of the *SAFEMCN*  package takes these output files (e.g., node_num_original.csv) as input to identify the N<sub>min</sub> threshold: the lowest AR1 value (marked by a rapid decrease followed by a rapid increase) is determined as N<sub>min</sub>, indicating the emergence of a topological plateau.   

Additionally, the AR1 module outputs a CSV file containing the AR1 coefficient for each sample size (e.g., node_num_original_window15_AR1.csv) and generates a scatter plot showing how the AR1 coefficient changes with sample size.   

Finally, the module also generates an AR1_min.csv file, which records the determined N<sub>min</sub> value and its corresponding AR1 value.  


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
