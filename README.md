# radial_plot
This is a Python script for generating a radial plot from fission track data following the procedure of Galbraith (1990).

The input data can be both ages, eg. for LA-ICP-MS method, or Ns and Ni counts for the external detector method (EDM). The data have to be saved in csv format and should have headers in the following format:

For inputting the ages:

**Age**	 | **1s**	 | **Composition**

For EDM method:

**Ns**	|  **Ni**	|  **Area (cm2)**	 | **Composition**	|  **Rho_d**	|  **Nd**	 | **zeta**	|  **zeta_1s**

The user can input any chosen compositional data (e.g. DPar, Cl wt%) and it is recommended to label the header with the used compositional data type.

Example input files are provided and can be downloaded together with the script.

---
*References:*

Galbraith, R. F. (1990). The radial plot: graphical assessment of spread in ages. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17(3), 207-214.
