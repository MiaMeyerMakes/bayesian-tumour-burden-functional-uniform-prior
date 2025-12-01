# Bayesian study on tumour burden using functional uniform priors in nonlinear mixed-effects models

**Author:** Mia Meyer

**Institution:** Stellenbosch University, Department of Statistics and Actuarial Science

**Degree:** Master of Commerce in Statistics and Data Science 

**Year:** 2025

**Supervisor(s):** Prof PJ Mostert and E Lesaffre


## Abstract

The choice of the noninformative prior for the model parameters in a Bayesian analysis of nonlinear (mixed) models has received significant attention in the literature. This thesis considers the use of a functional uniform prior (FUP) within nonlinear (mixed) models, specifically in dose-response and tumour growth inhibition (TGI) model applications. Traditional noninformative priors like uniform and Jeffreys’ priors are widely used in the pharmaceutical industry, however, they can be quite informative in nature when mapping them onto a nonlinear functional space. Additionally, Jeffreys' prior is dependent on the full data structure being available when deriving it in the context of clinical trials. Bornkamp (2012) derived the FUP for a few nonlinear regression models, including exponential, power and hyperbolic-Emax models, but did not consider nonlinear mixed models. An extensive Bayesian simulation study is conducted to evaluate the operating characteristics of the FUP when compared with these standard traditional priors. The Bayesian simulation study is extended to mixed-effects models, specifically the exponential one-parameter models and the two-parameter TGI model. Finally, the performance of the FUP is explored when analysing oncology data on colorectal cancer. 

The FUP has the theoretical advantages of being transformation invariant and satisfies the likelihood-principle. While the FUP approximates Jeffreys’ prior, it also has the advantage of being specified prior to data collection in contrast to Jeffreys’ prior.

## Repository Structure

```
.
|   README.md
+---Fixed-effects models
|       Emax-FUP.R
|       Emax-Jeffreys.R
|       Emax-Unif.R
|       Emax-uniform-coverage.R
|       Exp-FUP.R
|       Exp-Jeffreys.R
|       Exp-Unif.R
|       Exp-uniform-coverage.R
|       Power-Jeffreys.R
|       Power-MH.R
|       Power-Unif.R
+---Literature review
|       Literature-review-Puromycin.Rmd
|       Literature-review-Theophylline.Rmd
+---Mixed-effects models
|       growth_stan_mixed_FUP.R
|       growth_stan_mixed_JEFF.R
|       growth_stan_mixed_UNIF.R
|       mixed-tgi-fup.R
|       mixed-tgi-jeffreys.R
|       mixed-tgi-uniform-MEDIUM.R
|       mixed-tgi-uniform-NARROW.R
|       mixed-tgi-uniform-WIDE.R
|       shrinkage_stan_mixed_FUP.R
|       shrinkage_stan_mixed_JEFF.R
|       shrinkage_stan_mixed_UNIF.R
|       tgi-mixed-analyse-results.R
|       tgi-mixed-filter-chains.R
\---Prime study
        create-prime-subset-sample.R
        mixed-tgi-MHstan-fup-PRIME.R
        mixed-tgi-MHstan-jeff-PRIME.R
        mixed-tgi-MHstan-unif-PRIME.R
        plot-tgi-model-fit.R
        prime-tgi-results-analyses.R
```

## Installation

### Requirements

- R version 4.4.1 or higher
- Key dependencies: rstan, nimble, ggplot2, nmle

### Setup

```bash
# Clone the repository
git clone https://github.com/[username]/[repository-name].git
cd [repository-name]
```

## Usage

### Reproducing Main Results

[Provide clear instructions for reproducing key results from your thesis]

```bash
# Example command to run main experiment
python experiments/main_experiment.py --config config/default.yaml
```

### Key Scripts

- `src/algorithms/[algorithm_name].py`: Implementation of [Algorithm Name] described in Section X. X
- `experiments/[experiment_name].py`: Generates results for Figure X.X / Table X.X
- `src/models/[model_name].py`: Mathematical model from Chapter X

## Citation

If you use this code in your research, please cite:

```bibtex
@mastersthesis{Meyer2025,
  author = {Mia Meyer},
  title = {Bayesian study on tumour burden using functional uniform priors in nonlinear mixed-effects models},
  school = {University of Stellenbosch},
  year = {2025},
  type = {Master's thesis},
}
```


## Contact

Mia Meyer - miameyer15@gmail.com

[The thesis will be linked here when published]

## Acknowledgments

I extend my heartfelt appreciation to my supervisor, Prof. Paul Mostert, for his thoughtful guidance, unwavering support, and deep knowledge. His patience as I navigated the challenging journey of (trying to) master Bayesian methods, combined with his commitment to mentorship and excellence, has fostered my development both academically and personally. He has not only encouraged me to pursue scholarly stewardship but also motivated me to embrace every opportunity for learning and growth. I am particularly grateful for how he facilitated connections with a remarkable network of researchers spanning both academic and industry settings, including my co-supervisor Prof.  Emmanuel Lesaffre, Shuting Tang, and the team at Roche. The expert guidance I received from the researchers Dr. David Dejardin and Dr. Francois Mercier throughout my thesis journey, alongside numerous late-night Teams meetings focused on biostatistics, pharmacokinetic modeling, and tumour dynamics, demonstrated not only their technical expertise but also their unwavering dedication to their field. Their constructive feedback significantly influenced many critical elements of the modelling approaches and simulations presented in this research. Every challenge I encountered instilled in me humility, discipline, and an appreciation for diligent work, effective communication, and meaningful collaboration.

I extend my gratitude to all parties who facilitated my research visit to Leuven in 2025, namely, the Stellenbosch University International Office for financing my flights, Prof. Geert Verbeke and his colleagues at KU Leuven for covering my accommodation expenses and ensuring I felt welcomed, and the SU Department of Statistics and Actuarial Science for providing the bursary that supported my daily living costs during my time abroad.

Special recognition goes to the South African Reserve Bank (SARB) for enabling me to pursue my master's degree through the provision of a full scholarship.

I wish to thank my parents for their support in every conceivable way - from the countless _eat-and-run_ meals to always telling me that I can achieve any dream if I give it my all.

Above all, I am grateful that His grace is sufficient in my weakness (2 Corinthians 12:9) and that every good gift comes from Him (James 1:17); though I have failed many times throughout this journey, I am driven by the truth that _whatever you do, work heartily, as for the Lord and not for men_ (Colossians 3:23).

**Soli Deo Gloria**
```
