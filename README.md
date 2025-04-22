# ⚡ Automated System-wide Strength Estimation Tool (ASSET)

**ASSET** is a free and open-source tool built with Python and [PSS®E](https://new.siemens.com/global/en/products/energy/services/transmission-distribution-smart-grid/pss-software/pss-e.html), designed to automatically assess and analyze grid strength across large-scale power systems. It supports multiple system operating conditions and contingencies, and outputs results in a tabular format for easy interpretation.

---

## 🚀 Features

- ✅ **Automated grid strength computation** for a list of given buses or Points of Interconnection (POIs).
- ⚠️ **Identify critical network contingencies** impacting each bus/POI.
- 🔍 **Analyze the impact of pre-determined contingencies** on system strength at specified buses.
- 📍 **Identify suitable locations for grid strength devices** _(coming soon)_.
- 📈 **Optimize location and size of synchronous condensers** and other grid strengthening devices for improving system strength _(coming soon)_.


---

## 🖼️ Workflow Overview

![ASSET Workflow](input/asset_flowchart.png)


---

## 🧪 Applications

ASSET has been used in the analysis of real-world power systems, including:

1. [Grid strength studies for U.S. Eastern Interconnect](https://www.nrel.gov/docs/fy24osti/88003.pdf) 🌐
2. [Puerto Rico grid resiliencey studies](https://www.nrel.gov/docs/fy24osti/88615.pdf) 🇵🇷
3. [Grid strength studies for U.S. Western Interconnect (WECC) system](https://www.osti.gov/servlets/purl/2500279/) 🌐
4. [Subnational strategies to improve grid quality and reduce energy costs in Argentina](https://www.nrel.gov/docs/fy25osti/91767.pdf) 🇦🇷
5. Grid oscillation studies in India _(ongoing work)_ 🇮🇳

---

## 📚 Citation


If you use ASSET in your research or publications, please cite the following papers:

```bash
P. Sharma, L. Rese, B. Wang, B. Vyakaranam and S. Shah, "Grid strength analysis for integrating 30 GW of offshore wind generation by 2030 in the U.S. Eastern Interconnection," 22nd Wind and Solar Integration Workshop (WIW 2023), Copenhagen, Denmark, 2023, pp. 36-43, doi: 10.1049/icp.2023.2716.
```
```bash
P. Sharma and S. Shah, "Application of the Extra Element Theorem for Grid Strength Analysis in IBR-Dominated Systems," 2025 IEEE Power & Energy Society General Meeting (PESGM), Austin, Texas, USA, 2025
```
---

## 📝 License

ASSET is released under the **BSD License**.  Copyright © National Renewable Energy Laboratory.

**NREL Software Record of Invention**  
Pranav Sharma, Shahil Shah, Bin Wang, Leonardo Rese –  “Automated System‑wide Strength Evaluation Tool (ASSET)”.

---

## ⚙️ Installation

**Clone the repository**

```bash
git clone https://github.com/NREL/ASSET_Public.git
cd asset

```

---

## ✉️ Contact

For questions, feedback, or other inquiries, please reach out to **shahil.shah@nrel.gov**.

For our related work in the domain of Power systems stability, please visit:  [**Grid Impedance Scan Tool**](https://www.nrel.gov/grid/impedance-measurement)
         
---

## 🤝 Contributing

We welcome contributions of all kinds—code, documentation, testing, or feature suggestions.

1. **Fork** the repository and create your branch from `main`.
2. **Commit** your changes with clear messages.
3. **Open a Pull Request** (PR) describing what you changed and why.
4. One of the maintainers will review your PR, suggest any revisions, and merge when ready.

If you’d like to discuss an idea before coding, please open an **Issue** or email us directly at **shahil.shah@nrel.gov**.

> **New to Git or GitHub?** Check out the [GitHub Docs “Fork a repo” guide](https://docs.github.com/en/get-started/quickstart/fork-a-repo).

Thank you for helping make ASSET better!
