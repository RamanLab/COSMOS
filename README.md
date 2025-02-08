# COSMOS - COmmunity and Single Microbe Optimization System

COSMOS provides a computational framework for comparing **monocultures** and **co-cultures** in identical environments to identify the optimal microbial system. The analysis spans diverse environments, exploring factors such as **medium composition**, **community interactions**, and **carbon sources**.

---

## Key Scripts  
- **Masterscripts for Four Environments (effect of medium composition and oxygen availability):**
  - `aerobicRich.m`  
  - `aerobicMinimal.m`  
  - `anaerobicRich.m`  
  - `anaerobicMinimal.m`
 Initial inputs like models and parameters (environment and kinetic) can be changed by manipulating these files.

- **Effect of Interactions:**  
  `InteractionAnalysis.m`  

- **Effect of Carbon Source:**  
  `CSanalysis.m`  

- **dynamic Flux Balance Analysis (dFBA):**  
  `dFBACom.m`

- **Productivity Comparison:**  
  `ProductivityAnalysis.m`  

- **Yield Comparison:**  
  `YieldAnalysis.m`  
  All the results in the paper are based on productivity comparison. Alternatively one can use yield comparison if required.

- **Biomass Ratio Optimization:**  
  `optBiomassRatio.m`  

- **Experimental Validation:**  
  `validation_SoKp.m`
  
  To validate the algorithm by comparison with the experiments in https://doi.org/10.1186/s13068-023-02304-4

- **Collating data from all four environments for statistical analysis:**  
  `StatDataAnalysis.m` and `CS_StatDataAnalysis.m`
---

## Models Used  
All models were obtained from the **BiGG Models** and **BioModels** databases. Some annotations were adjusted to ensure consistent formatting across the dataset.

---

## Prerequisites  
- **MATLAB**
- **COBRA Toolbox**  
- **IBM ILOG CPLEX**

All simulations were conducted using **MATLAB R2018a**, **COBRA Toolbox v3.0**, and **IBM ILOG CPLEX 12.8**.

Please initialise the COBRA Toolbox and the solver before running COSMOS
We recommend the ibm_cplex solver as COSMOS uses *fastFVA* for efficiency.

---

## How to Use  
1. Clone the repository:  
   ```bash
   git clone https://github.com/COSMOS.git

