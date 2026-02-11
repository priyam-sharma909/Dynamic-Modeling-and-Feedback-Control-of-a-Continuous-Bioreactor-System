# Control for Bioreactor (Chemostat) System

This repository contains MATLAB files for:

- Open-loop simulation of the bioreactor  
- FOPTD model building  
- Closed-loop control (P/PI/PID) of 3 output variables  

---

## 1. Open-Loop Simulation

**Folder:** `Open Loop Model`

Run the files in this folder to simulate the nonlinear dynamics of the bioreactor (chemostat).

**Outputs:**

- Biomass concentration `X`  
- Substrate concentration `S`  
- Product concentration `P`  

These responses are used as data for model identification.

---

## 2. FOPTD Models

**Folder:** `FOPTD Models`

Builds three types of FOPTD models for the bioreactor:

1. Nonlinear fmincon-based FOPTD model  
2. Deep Neural Network (DNN) based FOPTD model  
3. Reaction Curve method-based FOPTD model (most accurate and used for controller tuning)

---

## 3. Closed-Loop Control

**Folder:** `Closed Loop Model`

Uses the selected FOPTD models to design and simulate:

- P controller  
- PI controller  
- PID controller  

for all three output variables (`X`, `S`, `P`).

Closed-loop control is implemented to regulate the bioreactor states and study the effect of different controller structures on system stability and performance.

---

## Results

Simulation and control results are documented in the attached presentation with detailed explanations, including:

- Comparison of open-loop and closed-loop responses  
- Controller performance for each output variable  
- Stability behavior and transient response analysis  
- Effect of controller tuning on system dynamics  




---
