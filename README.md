# 3-Omega Thermal Conductivity Calculator

A web-based tool for simulating and visualizing 3-omega voltage signals based on the Cahill integral formulation. The application computes the temperature oscillations and third harmonic voltage response (V₃ω) of a microfabricated wire heater used in thermal conductivity measurements.

---

## Features

- Simulates V₃ω (in-phase and out-of-phase) signals.
- Interactive D3.js plot of ln(2ω) vs. V₃ω and ΔT.
- Supports frequency domain filtering.
- Exports data as `.dat` files.
- Generates PDF reports with plots and parameters.

---

## Input Parameters

Each parameter is set via numeric input fields in the sidebar.

| Parameter                       | Symbol     | Units       | Input ID                |
|--------------------------------|------------|-------------|--------------------------|
| Total Wire Resistance          | R₀         | Ω           | `resistance`             |
| Wire Length                    | L          | mm          | `wire-length`            |
| Wire Width                     | 2b         | mm          | `wire-width`             |
| Thermal Coeff. of Resistance   | α          | 1/K         | `tcr`                    |
| Substrate Thermal Conductivity| κₛ         | W/mK        | `kappa-substrate`        |
| Substrate Thickness            | dₛ         | mm          | `thickness-substrate`    |
| Substrate Thermal Diffusivity | Dₛ         | mm²/s       | `diffusivity-substrate`  |
| Heating Current (peak)        | I₀         | A           | `current`                |

---

## Output

**Calculated Values:**

- Power per unit length: `P_rms/L`
- Linearity domain limits:
  - `f_min = (25·Dₛ)/(4·π·dₛ²)`
  - `f_max = Dₛ / (100·π·b²)`

**Plot:**

- X-axis (bottom): ln(2ω)
- X-axis (top): Frequency (Hz, log scale)
- Y-axis (left): ΔT (K)
- Y-axis (right): V₃ω (RMS) [V]

**Plot lines:**

- Blue: V₃ω in-phase (real part)
- Orange: V₃ω out-of-phase (imaginary part)
- Thick lines: linear domain only

---

## Export Options

- **Save .dat**: Download tabulated simulation data in plain text.
- **Save .pdf**: Generate a report with parameters and plot.

---

## Calculation Method

Implements Cahill's 3-omega integral model:
- Complex integration of a sinc²-weighted kernel.
- Frequency sweep from 0.02 Hz to 10⁵ Hz (log scale).
- Integration range up to η_max = 50/b using 1000-point Simpson’s rule.

Main computed quantities:
- Temperature oscillation: ΔT(ω)
- Third-harmonic voltage (RMS): V₃ω(ω) = α·ΔT·V₁/√2

---

## Dependencies

- [D3.js v7](https://d3js.org)
- [jsPDF](https://github.com/parallax/jsPDF)
- [html2canvas](https://html2canvas.hertzen.com)

---

## Author

**sim_3-omega, v0.9**  
NitaD — Univ Paris-Saclay

---

## License

This is a research tool. No warranty implied.

