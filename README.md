# 🛠️ Forked Ideal Gas Code – Based on ANSYS Original

## 📄 Copyright Notice

**Original Copyright © 1988–1998 ANSYS, Inc.**  
**Forked and Enhanced by Arash Jalil-Khabbazi**

All Rights Reserved.

This repository contains a forked and extended version of proprietary code originally developed by **ANSYS, Inc.** for simulating thermodynamic and transport properties of ideal gases. The foundational structure of this code is **unpublished proprietary software** and is **protected by U.S. copyright law**.  
It was furnished under a written license agreement and is considered **confidential** by ANSYS, Inc. Usage, reproduction, or disclosure is permitted only in accordance with the terms of that agreement.

---

## 📌 About This Fork

This fork has been modified and enhanced by **Arash Jalil-Khabbazi** to demonstrate the computational structure and implementation of an **ideal gas model** within the ANSYS FLUENT framework using **User-Defined Functions (UDFs)**.

### 🔍 Key Features

- ✅ **Ideal Gas Thermodynamics**: Computes density, enthalpy, entropy, specific heat, and speed of sound based on temperature, pressure, and species composition.
- ✅ **Transport Properties**: Includes empirical models for viscosity and thermal conductivity.
- ✅ **Species-Specific Modeling**: Currently supports a binary gas mixture (Hydrogen and Methane).
- ✅ **Modular Design**: Functions are organized for clarity and ease of integration with FLUENT’s real gas interface.
- ✅ **Reference State Handling**: Includes reference enthalpy and entropy offsets for consistency with thermodynamic tables.

> ⚠️ **Disclaimer:** This code is intended for academic and illustrative purposes only. It simplifies many real-gas behaviors and is not intended for production-level simulations. More advanced models (e.g., Peng-Robinson, Redlich-Kwong) are recommended for high-accuracy applications.

---

## 🧪 Supported Species

By default, this code is configured for a **binary gas mixture**:
- **Hydrogen (H₂)**
- **Methane (CH₄)**

These species are defined in the `MIXTURE_Setup()` function and their properties (molecular weight, critical temperature/pressure, specific heat coefficients) are hardcoded in the respective initialization functions.

---

## 🔧 How to Add More Species

To extend the model to support more gases:

1. **Update the species count**:
   ```c
   #define n_specs 3  // or more
   ```
