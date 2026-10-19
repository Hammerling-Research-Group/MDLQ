# MDLQ: Multisource detection, localization, and quantification of methane emissions

This repository contains the R code and data for the multisource detection, localization, and quantification (MDLQ) model described in:

> Daniels, W. S., Nychka, D. W., and Hammerling, D. M. (2026). A Bayesian hierarchical model for methane emission source apportionment. *The Annals of Applied Statistics*, 20(3), 2195–2213. https://doi.org/10.1214/26-AOAS2202

## Overview

MDLQ performs source apportionment on oil and gas sites: given methane concentration measurements from a network of point-in-space continuous monitoring system (CMS) sensors and a set of potential emission sources whose locations are known, it estimates which sources are emitting and at what rate. Emissions from different sources are allowed to overlap in time, a key difference from single-source inversion methods, which assume that only one source emits at a time.

The model treats the background-removed, minute-frequency CMS sensor observations as a linear combination of concentrations simulated by an atmospheric transport model, one set of simulations per potential source. The emission rate of each source has a spike-and-slab prior: a point mass at zero (the source is not emitting) and an exponential distribution on the positive reals (the source is emitting). The model errors follow an AR(1) process, and an inverse-gamma prior on the error variance makes the errors Student's t-distributed a priori. The posterior is sampled with a Metropolis-within-Gibbs algorithm. To accommodate emission rates that change over time, the model is applied independently to consecutive 30-minute inversion windows, producing a piecewise-constant emission time series for every source. Each window returns, for each source with sufficient information, a posterior mean emission rate and a 95% credible interval (kg/hr) as well as the posterior probability that the source is emitting. The paper describes the model, the background removal procedure, and the evaluation in full; the derivations of the conditional distributions used in the sampler are in Section S3 of the supplementary material (https://doi.org/10.1214/26-AOAS2202SUPPA).

The code here applies MDLQ to the ADED 2024 controlled release experiment conducted at the Methane Emissions Technology Evaluation Center (METEC) in Fort Collins, Colorado. The experiment consists of 337 controlled releases from five pieces of equipment, measured by 10 CMS sensors, between February 5 and April 30, 2024. This repository contains the MDLQ model, the script that applies it to the METEC data, and the script that produces the main results and figures for the METEC experiment (Section 4.1 of the paper). It does not contain the code for the simulation study (Section 4.2), for the comparison methods other than the two MDLQ variants described below (Section 4.3), or for the sensitivity studies in the supplementary material.

**Note:** this README was generated with Claude Sonnet 5.5.

## Repository structure

```
MDLQ/
├── code/
│   ├── MAIN_1_MDLQ.R                                  # runs MDLQ on the CMS data
│   ├── MAIN_2_analyze_results.R                       # evaluates the output against ground truth, makes figures
│   ├── HELPER_MCMC.R                                  # the MDLQ Gibbs sampler, run.mdlq.mcmc()
│   ├── HELPER_MCMC_variant1_independent_errors.R      # sampler variant: independent (non-autocorrelated) errors
│   ├── HELPER_MCMC_variant2_laplace_errors.R          # sampler variant: independent Laplace errors
│   └── HELPER_functions.R                             # background removal
├── input_data/
│   ├── sensor_locations.csv
│   ├── source_locations.csv
│   ├── input_data_SAMPLE.RData
│   └── ground_truth.RData
├── output_data/
│   ├── MDLQ_output_30min_interval_30min_step_SAMPLE.RData
│   ├── MDLQ_output_30min_interval_30min_step.RData
│   ├── forward_simulations_SAMPLE.RData
│   ├── forward_simulations.RData
│   └── pregenerated_*.RData                           # ground truth and information masks, saved to save time
├── figures/                                           # populated by MAIN_2_analyze_results.R
├── LICENSE
└── README.md
```

## Data

### What MDLQ needs

Running MDLQ requires two kinds of data: (1) methane concentration observations from the CMS sensors, in ppm, and (2) simulated concentrations at each sensor location from an atmospheric transport model, with one set of simulations for every potential emission source on the site. The simulations in this repository were produced with the Gaussian puff model of Jia et al. (2025), https://doi.org/10.1038/s41598-025-99491-x, using the code at https://github.com/Hammerling-Research-Group/FastGaussianPuff. They are generated for a unitary emission rate of 1 g/s, so the estimated coefficients are emission rates in g/s, which `MAIN_1_MDLQ.R` converts to kg/hr. The locations and heights of the sensors and sources are needed to produce these simulations but not to run MDLQ itself, and are provided in `input_data/sensor_locations.csv` and `input_data/source_locations.csv`.

### Input files

| File | Contents |
|---|---|
| `input_data/sensor_locations.csv` | Latitude, longitude, and height (2 m) of the 10 CMS sensors. |
| `input_data/source_locations.csv` | Latitude, longitude, and height (m) of the five potential emission sources: West Wellhead, West Separator, Tanks, East Wellhead, East Separator. |
| `input_data/input_data_SAMPLE.RData` | A little over seven days (10,321 minutes) of minute-frequency CMS observations from mid-February 2024 and the corresponding forward simulations. |
| `input_data/ground_truth.RData` | Start time, end time, and release rate (kg/hr) of each controlled release for each piece of equipment, with one column per equipment ID (`EG4W` = West Wellhead, `EG4S` = West Separator, `EG4T` = Tanks, `EG5W` = East Wellhead, `EG5S` = East Separator). Missing values denote equipment that was not releasing. This file is not needed to run MDLQ, only to evaluate it in `MAIN_2_analyze_results.R`. |

Due to the proprietary nature of the CMS sensor technology used in this study, only a sample of the concentration data from the experiment can be shared publicly, which is why `input_data_SAMPLE.RData` covers one week rather than the full experiment. The complete set of CMS measurements can be made available for research purposes upon request to the corresponding author (see Contact).

`input_data_SAMPLE.RData` holds a named list with the elements `times` (minute-frequency time stamps), `obs` (a matrix of concentration observations with one row per minute and one column per sensor), `WD` and `WS` (wind direction and speed), and then one matrix of simulated concentrations per source, with the same dimensions as `obs` and named after the source. The sensor columns of `obs` and of each simulation matrix are the 10 sensors listed in `sensor_locations.csv`. The simulation matrices must come after the first four elements, and `first.sim.ind` in `MAIN_1_MDLQ.R` (set to 5) tells the script where they start. The wind data are carried along in this object but are not used by MDLQ itself; wind enters through the forward simulations. To apply MDLQ to a different site, build an object with this structure and adjust `first.sim.ind` if necessary. `MAIN_1_MDLQ.R` works for any number of sources and sensors, but `MAIN_2_analyze_results.R` is specific to the METEC experiment: it assumes the five sources are ordered West Wellhead, West Separator, Tanks, East Wellhead, East Separator, which determines the labels and colors in its figures.

Note that although these files use the `.RData` extension, they were written with `saveRDS()` and must be read with `readRDS()`, not `load()`.

### Output files

`MDLQ_output_30min_interval_30min_step_SAMPLE.RData` contains the MDLQ results for the sample week. `MDLQ_output_30min_interval_30min_step.RData` (about 67 MB) contains the results for the entire experiment, which cannot be recomputed from the public data but can be passed through `MAIN_2_analyze_results.R` to reproduce the main results of the paper. The accompanying `forward_simulations*.RData` files hold the forward simulations used to produce each output file, and are read by `MAIN_2_analyze_results.R`. The files starting with `pregenerated_` hold the ground truth on the 30-minute inversion windows and the information masks (see "Implementation details"), which are slow to compute; they can be regenerated with the toggles at the top of `MAIN_2_analyze_results.R`.

The `MDLQ_output` file is a list with elements `times`, `source.names`, and `out`, where `out` has one entry per inversion window. For a window that was run, that entry is a list with the following elements.

| Element | Description |
|---|---|
| `q.hat`, `q.hat.lower`, `q.hat.upper` | Posterior mean and 95% credible interval of the emission rate (kg/hr) for each source |
| `zs` | Posterior mean of the emission indicator for each source, i.e., the posterior probability that the source is emitting |
| `rates` | Retained MCMC samples of the emission rates (kg/hr) for the sources that were included in the inversion |
| `rs` | Retained MCMC samples of the AR(1) coefficient |
| `interval.start`, `interval.end` | First and last time stamps in the window |

Sources that were not part of a given inversion, or for which the sampler did not return a result, carry a text label rather than a number: `"no info"` (no sensor downwind of the source, see below) or `"DNC"` / `"broke on betas"` (the sampler failed). A window in which any sensor observation remained missing after interpolation is stored as the single string `"NA y"`. `MAIN_2_analyze_results.R` converts all of these to `NA`.

## Requirements

The code is written in R and uses the following packages. `MAIN_1_MDLQ.R` uses `lubridate`, `foreach`, `doParallel`, and `zoo`. `HELPER_MCMC.R` uses `Matrix`, `TruncatedNormal`, and `R.utils`; the independent-errors variant additionally uses `MASS` and `tmvtnorm`. `MAIN_2_analyze_results.R` uses `ald`, `fields`, `scales`, `lubridate`, and `zoo`, plus `rstudioapi` when run in RStudio.

Both main scripts also source helper functions directly from GitHub at run time, so an internet connection is required. Both load `HELPER_spike_detection_algorithm.R` from https://github.com/wsdaniels/DLQ, which contains the spike detection algorithm of Daniels, Jia, and Hammerling (2024), https://doi.org/10.1525/elementa.2023.00110, used for background removal. `MAIN_2_analyze_results.R` additionally loads `helper_functions.R` from https://github.com/wsdaniels/CMS-durations, which is used to construct the information mask.

## Usage

All paths in the scripts are relative to the `code` directory, so run them with `code` as the working directory. `MAIN_2_analyze_results.R` sets this automatically when run from RStudio; otherwise, use `setwd()` or run from a terminal inside `code`.

**1. Run MDLQ.** Open `code/MAIN_1_MDLQ.R` and review the settings in the user input block at the top.

| Setting | Default | Meaning |
|---|---|---|
| `num.cores.to.use` | 7 | Number of cores used to process inversion windows in parallel; set to 1 to run in serial |
| `interval.length` | 30 | Length of each inversion window in minutes |
| `step.size` | 30 | Number of minutes the window advances; set equal to `interval.length` for non-overlapping windows |
| `input.file.path` | `../input_data/input_data_SAMPLE.RData` | Input data object |
| `output.file.path` | `../output_data/MDLQ_output_30min_interval_30min_step_SAMPLE.RData` | Where results are saved |
| `sim.save.path` | `../output_data/forward_simulations_SAMPLE.RData` | Where the forward simulations are saved |
| `spike.slab.regression.path` | `../code/HELPER_MCMC.R` | Which sampler to use (see below) |
| `first.sim.ind` | 5 | Index of the first forward simulation in the input data object |

Then run the script. On the sample week, it takes about 20 minutes with 7 cores on an Apple MacBook Pro with an M1 chip. `MAIN_1_MDLQ.R` does not set a random seed, so the MCMC output will differ slightly from run to run, including from the output files provided. Running MDLQ on the full experiment requires the complete CMS data (see Data); in that case, change the three file paths in the table above so that the sample files are not overwritten.

**2. Analyze the results.** Open `code/MAIN_2_analyze_results.R`, set `run.sample <- T` to analyze the sample week or `run.sample <- F` to analyze the full experiment from the output files provided, and run the script. With the default settings (`generate.truth <- F` and `generate.info.mask <- F`) it reads the pregenerated ground truth and information masks rather than recomputing them. It writes three figures to `figures/` and prints the summary statistics to the console. If you changed `interval.length` or `step.size` in step 1, change them here to match; if the window length is longer than the step size, the estimates from overlapping windows are averaged.

To reproduce the results in the paper, use `run.sample <- F`. The outputs correspond to the paper as follows.

| Output | Corresponds to |
|---|---|
| `figures/loc_quant_results.png` | Figure 6: source- and site-level inventories with 95% credible intervals, the distribution of site-level quantification errors, and the number of correct localization estimates per window |
| `figures/loc_quant_results_info_adjusted.png` | The same three panels computed only on the inversion windows in which every source has information for more than 95% of the window, i.e., the analysis restricted to windows with a downwind sensor for all sources (Section 4.1) |
| `figures/result_time_series.png` | Figure 5: true and estimated emission state over time (for the full experiment, a five-day period as in the paper; for the sample, the whole week) |
| Console output | The inventory errors, mean site-level quantification error and its inner 95% range, coverage of the site-level credible intervals, and the detection metrics in Table 2 |

## Implementation details

**Background removal.** MAIN_1 first linearly interpolates gaps of up to 10 minutes in the sensor observations and then removes the methane background from each sensor using the gradient-based spike detection algorithm described in Section 2.3 of the paper. Observations outside of detected spikes are set to zero, and within a spike the background is the mean of the observations immediately before and after it. Negative values are set to zero.

**Windowing and information content.** For each inversion window, the sensor observations and simulations are concatenated across sensors to form the vectors and design matrix described in Section 3.1. A source is included in the inversion only if it has information, meaning that more than four of its simulated concentrations in the window exceed 0.25; otherwise its result is `"no info"`. If there are fewer than five positive observations in the window and all of them are below 1 ppm, the emission rate of every included source is set to zero without running the sampler.

**Sampler.** The sampler in `HELPER_MCMC.R`, `run.mdlq.mcmc()`, is run for 1,750 iterations with the first 250 discarded as burn-in. It is initialized at the ordinary least squares estimate, with nonpositive coefficients replaced by a small positive value. The hyperparameters default to the values used in the paper (1 for all of `a.vec`, `b.vec`, `c.vec`, `d.vec`, `alpha1`, and `alpha2`) and can be changed through the arguments of `run.mdlq.mcmc()`. The sampled emission rates are in g/s and are multiplied by 3.6 to give kg/hr.

**Information mask.** The information mask used in `MAIN_2_analyze_results.R` is separate from the check above and is computed after the fact from the forward simulations. For each source, `create.info.mask()` flags the minutes at which at least one sensor is predicted to observe an enhancement from that source, clustered into events lasting at least 15 minutes. The fraction of flagged minutes in each inversion window is then used to restrict the analysis to windows in which all sources have information.

**Inventories and alerts.** `MAIN_2_analyze_results.R` builds source- and site-level inventories by sampling from the retained posterior draws in each window, multiplying by the window length, and summing over windows (reported in metric tons). For windows in which a source has no information or the sampler failed, draws are taken from that source's estimates in other windows (Section 3.5). A source is classified as emitting in a window if its posterior mean emission indicator exceeds that source's average across windows, and is classified as not emitting if its estimated rate is exactly zero.

## MCMC variants

Two alternative samplers are included for comparison with the default model. They expose the same `run.mdlq.mcmc()` interface, so to use one, point `spike.slab.regression.path` in `MAIN_1_MDLQ.R` to the file and set `output.file.path` so that the results are not written over the default output.

- `HELPER_MCMC_variant1_independent_errors.R` fixes the autocorrelation coefficient at zero, giving a model with independent errors. The paper's supplementary material (Section S8) explains why modeling the autocorrelation matters.
- `HELPER_MCMC_variant2_laplace_errors.R` replaces the Gaussian error model with independent Laplace errors, motivated by the heavy-tailed errors that arise when a concentration enhancement is predicted but not observed, or vice versa. This is the "MDLQ Laplace errors" column of Table 3 and is defined in Section S9 of the supplementary material.

## Citation

If you use this code in your research, please cite the accompanying paper:

Daniels, W. S., Nychka, D. W., and Hammerling, D. M. (2026). A Bayesian hierarchical model for methane emission source apportionment. *The Annals of Applied Statistics*, 20(3), 2195–2213. https://doi.org/10.1214/26-AOAS2202

```bibtex
@article{daniels_bayesian_2026,
	title = {A {Bayesian} hierarchical model for methane emission source apportionment},
	volume = {20},
	issn = {1932-6157, 1941-7330},
	doi = {10.1214/26-AOAS2202},
	number = {3},
	journal = {The Annals of Applied Statistics},
	publisher = {Institute of Mathematical Statistics},
	author = {Daniels, William S. and Nychka, Douglas W. and Hammerling, Dorit M.},
	month = sep,
	year = {2026},
	pages = {2195--2213},
}
```

The Gaussian puff atmospheric dispersion model used to generate the forward simulations is described in Jia, M., Fish, R., Daniels, W. S., Sprinkle, B., and Hammerling, D. (2025). A fast and lightweight implementation of the Gaussian puff model for near-field atmospheric transport of trace gasses. *Scientific Reports*, 15, 18710. https://doi.org/10.1038/s41598-025-99491-x

## License

This project is released under the MIT License; see [LICENSE](LICENSE).

## Contact

William Daniels (wdanie16@jh.edu). Please also contact the author to request access to the complete CMS concentration data from the METEC experiment for research purposes.
