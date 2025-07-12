# MOSFET Modeling Tutorial – IRF740 (MATLAB/Simulink)

This repository contains a **complete hands-on tutorial** on the modeling of the **IRF740 power MOSFET** using **MATLAB/Simulink**.  
It is designed for **first-year Master’s students in Electronic Engineering**, but is also suitable for researchers and educators interested in power device modeling.

---

## About the Tutorial

The tutorial is based on the scientific article:

> **Modeling Power Electronics in Simulink. A Tutorial for MSc Electronic Engineering Students**  
> Leonardo Motta, Nicolò Busi — University of Genoa  

The tutorial guides learners through:
- Static modeling of the IRF740 (I-V characteristics)
- Dynamic modeling with parasitic capacitances
- Thermal modeling using an RC-network
- Safety mechanisms like thermal shutdown

---

## Quick Start

1. Open MATLAB and navigate to the desired model folder.
2. Open the Simulink `.slx` file (e.g. `static_model.slx`).
3. Run the simulation.
4. Modify input signals or parameters to explore different behaviors (e.g. VGS, VDS, ambient temperature).

---

## Learning Objectives

By following this tutorial, you will:
- Understand how to extract and use datasheet parameters in modeling
- Analyze the impact of parasitic elements on switching behavior
- Implement a temperature-dependent model using Simulink and MATLAB scripting
- Simulate thermal shutdown and temperature-sensitive electrical parameters (e.g., RDS(on), VTH)

---

## Example Results

Simulations include:
- **IDS-VDS curves** for various VGS levels
- **Switching waveforms** with and without parasitic capacitances
- **Thermal behavior over time**, showing RDS(on) and VTH variations
- **Power dissipation** and thermal shutdown conditions

---

## Requirements

- MATLAB R2021a or later
- Simulink
- Control System Toolbox (for some simulations)
- Basic knowledge of power electronics and MOSFET operation

---

## Contributions

This repository is open to contributions! Feel free to fork and extend the models for:
- Other power devices (e.g., SiC or GaN MOSFETs)
- Advanced thermal modeling (e.g., convection, aging)
- Improved visualizations and GUI for educational purposes


