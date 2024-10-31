# COSMOS - COmmunity and Single Microbe Optimization System

COSMOS provides a computational framework for comparing **monocultures** and **co-cultures** in identical environments to identify the optimal microbial system. The analysis spans diverse environments, exploring factors such as **medium composition**, **community interactions**, and **carbon sources**.

---

## Key Scripts  
- **Masterscripts for Four Environments (effect of medium composition and oxygen availability):**
  - `aerobicRich.m`  
  - `aerobicMinimal.m`  
  - `anaerobicRich.m`  
  - `anaerobicMinimal.m`  

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

- **Biomass Ratio Optimization:**  
  `optBiomassRatio.m`  

- **Experimental Validation:**  
  `validation_SoKp.m`  

---

## Models Used  
All models were obtained from the **BiGG Models** and **BioModels** databases. Some annotations were adjusted to ensure consistent formatting across the dataset.

---

## Prerequisites  
- **MATLAB**
- **COBRA Toolbox**  
- **IBM ILOG CPLEX**

All simulations were conducted using **MATLAB R2018a**, **COBRA Toolbox v3.0**, and **IBM ILOG CPLEX 12.8**.

Please use the *fluxVariability* function instead of *fastFVA* for other solvers.

---

## How to Use  
1. Clone the repository:  
   ```bash
   git clone https://github.com/COSMOS.git

