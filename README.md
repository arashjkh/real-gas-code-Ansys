# 🔬 Ideal Gas Mixture Property Model for ANSYS Fluent

## 📌 Overview

This repository contains a forked and extended version of a proprietary thermodynamic model originally developed by **ANSYS, Inc.** for simulating **ideal gas mixtures** in **ANSYS Fluent** using a **User-Defined Function (UDF)**.

The forked version by **Arash Jalil-Khabbazi** expands upon the original structure, demonstrating how to calculate thermophysical properties of multi-species gas mixtures based on polynomial fits for specific heat and critical property data. The sample configuration currently includes **hydrogen (H₂)** and **methane (CH₄)**.

---

## 🧠 What This Code Does

This UDF provides Fluent with custom implementations for the following gas properties of mixtures:

- ✅ Density  
- ✅ Specific heat (Cp)  
- ✅ Enthalpy and entropy  
- ✅ Speed of sound  
- ✅ Viscosity and thermal conductivity  
- ✅ Temperature and pressure derivatives of enthalpy and density  
- ✅ Mixture molecular weight and gas constant  

All properties are computed assuming **ideal gas behavior** using user-defined polynomial coefficients and static inputs.

---

## 🔧 How It Works

- Written in **C**, using Fluent's UDF API (`udf.h`)
- Compatible with **FLUENT’s compiled UDF libraries**
- Designed to work in **SI units**
- Currently supports **2 gases**:  
  - H₂ (Hydrogen)  
  - CH₄ (Methane)  
- Can be **manually extended** to support more species by editing:
  - `#define n_specs`  
  - `Mw()`, `Cp_Parameters()`, `Tcrit()`, `Pcrit()`, `Vcrit()`

---

## 📂 Application

This model is ideal for use in:
- Academic CFD research
- Educational demos of ideal gas modeling
- Sensitivity studies comparing Fluent’s built-in models vs. UDFs
- Simple multi-species non-reacting flows (e.g., H₂/CH₄ mixing, pipeline studies)

It is **not intended** for high-fidelity reacting flow simulations or real-gas effects beyond ideal assumptions.

---

## ⚠️ Limitations

- ❌ Only 2 species are supported by default (can be extended manually)
- ❌ Assumes **ideal gas behavior** — no non-ideal EOS support
- ❌ All inputs (Cp coefficients, critical properties) are hardcoded
- ❌ Requires recompilation when changing species or parameters

---

## 📜 Legal and Copyright Notice

**Original Copyright © 1988–1998 ANSYS, Inc.**  
**Forked and Extended by Arash Jalil-Khabbazi**  
All Rights Reserved.

This repository includes a modified version of proprietary source code originally developed by **ANSYS, Inc.** The base code is protected under U.S. copyright law as an **unpublished work** and is furnished under a written license agreement. It is considered **confidential** and may **not be used, copied, or disclosed** except in accordance with the terms of that agreement.

> This fork is intended **solely for academic and demonstrative purposes** and does not distribute or disclose any part of the ANSYS FLUENT binary or internal globals.

---

## 📬 Contact

For questions or additional information, please contact:  
**📧 arashjkh@gmail.com**

---

