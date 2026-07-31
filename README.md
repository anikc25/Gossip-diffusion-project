# ECO613M Course Project — Gossip vs. Elder Seeding for Information Diffusion

Course project for **ECO613M**, replicating and extending analysis from *"Using Gossips to Spread Information: Theory and Evidence from Two Randomized Controlled Trials"* (Banerjee, Chandrasekhar, Duflo & Jackson).

The project examines how the choice of "seed" households (gossips vs. village elders vs. random) affects how far information spreads through a village network, using data from the paper's two RCTs.

## Contents

```
├── Submission_Code/
│   └── ECO613M_Submission_Code.ipynb   # Main analysis notebook
├── PPT/
│   └── ECO613M_PPT.pptx.pdf            # Project presentation slides
└── rdz008_supplementary_data/          # Original replication data & appendix
    └── Supplementary/
        ├── Final_Data/                 # Data tables used in the analysis
        ├── appendixfigures/
        ├── appendixtables/
        └── MS24087-SupplementaryAppendix.pdf
```

## Analysis Notebook

`Submission_Code/ECO613M_Submission_Code.ipynb` reproduces the paper's key figures and tables and adds an original extension:

- **Figures 1–3** — diffusion patterns by treatment arm (gossip / elder / random)
- **Tables 1, 2, 3, 6–9** — regression and comparative results across seeding strategies
- **Extension — Predicting Optimal Seeds Without Network Data**: trains a machine-learning model to predict each household's diffusion centrality from cheap-to-collect observables (demographics, self-reported interactions) instead of full network data, then tests whether ML-predicted seeds outperform randomly chosen seeds.

Methods used: OLS, IV (2SLS), quantile regression, and gradient-boosted trees (XGBoost), via `statsmodels`, `linearmodels`, and `scikit-learn`.

## Data Source

The data in `rdz008_supplementary_data/` is the publicly released replication package accompanying the original paper (submission ID rdz008 / MS24087), included here for reproducibility. Credit for the underlying data belongs to the original authors.

## Requirements

```
numpy
pandas
scipy
statsmodels
linearmodels
scikit-learn
xgboost
```

Install with:
```bash
pip install numpy pandas scipy statsmodels linearmodels scikit-learn xgboost
```

## Author

Anik Chakraborty, 
Shamoitree Pal
