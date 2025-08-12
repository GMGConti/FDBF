# Fracture Displacement Basis Function Method

## Overview

This Code demonstrates the **Fracture Displacement Basis Function Method** for efficient modeling of shear displacement and tensile opening in 2D fractured domains. It computes normalized basis functions for individual fractures assuming elliptical shear slip and tensile opening profiles, and uses them to simulate displacement and stress fields of fracture pattern under far-field. Rock parameters (first and second Lame constants), fluid pressure, far field stress, friction coefficients, and the fracture pattern are taken as inputs.

## Authors
FDBF was implemented by Giulia Conti and Patrick Jenny at the [Institute of Fluid Dynamics](http://www.ifd.mavt.ethz.ch), ETH Zurich. The FDBF method was introduced by Giulia Conti, Stephan K. Matthai, and Patrick Jenny in the reference mentioned below.

## License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).
Further, any research making use of this software should appropriately cite the reference given below, in keeping with normal academic practice.

## Reference
```LaTeX
@article{Conti_2025,
    title = {Fracture Displacement Basis Function (FDBF) Method for Efficient Geomechanical Calculations of Fractured Rock},
    volume = {},
    doi = {},
    journal = {International Journal of Rock Mechanics and Mining Sciences},
    author = {Giulia Conti and Stephan K. Matth\"ai and Patrick Jenny},
    year = {2025},
}
```

## Theory
#### Normailzed Fracture Basis Functions
For the shear normalized basis function, solve the force balance

$$
    \frac{\partial}{\partial x_j} \hat{\sigma}_{ij}= 0
$$

for the displacement field $u$

$$
\textrm{on   }\mathbb{R}^2\textrm{ with } \\ \lim_{|\bf x|\rightarrow\infty}\hat{u}_i=0
\ \ \ \wedge\ \ \ 
\left(\hat{u}_1,\ \hat{u}_2\right)^T
\ =\ 
\left(\pm(0.5\sqrt{1-x_1^2/0.5^2}),\ 0\right)^T
\textrm{ for }
\left(x_1,x_2\right)^T\in\left([-0.5,0.5],\ \pm\epsilon\right)^T
\textrm{ with }0<\epsilon\ll 1
$$

assuming linear elasticity

$$
\hat{\sigma}_{ij}=\lambda \frac{\partial \hat{u}_k}{\partial x_k} \delta\_{ij}
+G(\frac{\partial\hat{u}_j}{\partial x_i}
+\frac{\partial\hat{u}_i}{\partial x_j})
\ \ \ 
$$

$$
\rightarrow \frac{\partial}{\partial x_j}(\lambda\frac{\partial \hat{u}_k}{\partial x_k} \delta\_{ij}+G(\frac{\partial \hat{u}_j}{\partial x_i}+\frac{\partial \hat{u}_i}{\partial x_j}))=0
\ \ \ 
$$

The tensile opening basis function is calculated analogously. 

#### Specific Basis Functions

$$
{\bf\hat{u}}^{f}({\bf x})
={\bf R}^{f}{\bf\hat{u}}({\bf G}^{f}({\bf x}-\bar{\bf x}^{f})/L^{f})
$$

$$
{\bf\hat{\sigma}}^{f}({\bf x})
={\bf R}^{f}{\bf\hat{\sigma}}({\bf G}^{f}({\bf x}-\bar{\bf x}^{f})/L^{f}){\bf G}^{f}
\ \ \ 
$$

with

$$
\ \ \ 
\bf R^f=\left(
\begin{array}{cc}
\cos(\alpha^f)&-\sin(\alpha^f)\\
\sin(\alpha^f)&\cos(\alpha^f)
\end{array}
\right)
\ \ \ 
$$

and

$$
\ \ \ 
\bf G^f=\left(
\begin{array}{cc}
\cos(\alpha^f)&\sin(\alpha^f)\\
-\sin(\alpha^f)&\cos(\alpha^f)
\end{array}
\right)
$$

#### Normalized traction- and normal stresses induced by fracture f on fracture g:

$$
{\bf\hat\sigma}\_{t,n}^{f\rightarrow g} = \int\_{\Omega^g}{\bf\hat\sigma}\_{t,n}^{f} dL
$$

and

$$
\ \ \ 
{\bf\sigma}\_{t,n}^{\infty\rightarrow g} = \int\_{\Omega^g} {\bf\sigma}\_{t,n}^{\infty} dL
$$

#### Constraints to obtain a linear system for the shear displacement dofs $s_t^f$:

$$
\Rightarrow \forall g\in\{1,..,n\}:\ \ \
$$

$$
\sum\_{f=1}^n (s_t^f (\hat{\sigma}_t\^{f\rightarrow g} + \mu^g \hat{\sigma}_n\^{f\rightarrow g}))
\le
\mu^g (\sigma_n\^{\infty\rightarrow g}\- p^g L^g ) \- \sigma_t\^{\infty\rightarrow g}
\ \ \ 
$$

(linear system for $$s_t^f$$)

#### Constraints to obtain a linear system for the dofs $s_t^f$ and $s_n^f$:
If tensile opening is present, there are two degrees of freedom per fracture: shear displacement $s_t$, and tensile opening $s_n$ . 

$$
\Rightarrow \forall g\in\{1,..,n\}:\ \ \ 
$$

$$
\sum\_{f=1}^n (s_t^{f} \hat{\sigma}\_{slip,t}\^{f\rightarrow g} + s_n^{f} \hat{\sigma}\_{open,t}\^{f\rightarrow g})
\=- \sigma_t\^{\infty\rightarrow g}
\ \ \ 
$$

(linear system for $$s_t^f$$)

$$
\sum_{f=1}^n (s_t^f \hat{\sigma}\_{slip,n}\^{f\rightarrow g}+ s_n^f \hat{\sigma}\_{open,n}\^{f\rightarrow g})
\le p^g L^g \-\sigma_n\^{\infty\rightarrow g}\
\ \ \ 
$$

(linear system for $$s_n^f$$)

#### Displacement- and stress fields induced by far field stress and N fractures

$$
{\bf u}\ =\ \sum\_{f=1}^N\{s_t^{f}{\hat{\bf u}}\_{slip}^{f}+s_n^{f}{\hat{\bf u}}\_{open}^{f}\}\ +\ {\bf u}^\infty
$$

and

$$
\sigma=\sum_{f=1}^N\{s_t^{f}{\hat{\sigma}}\_{slip}^{f}+s_n^{f}{\hat{\sigma}}\_{open}^{f}\}\ +\ \sigma^\infty
$$

with far field displacement obtained from far field stress:
- given: $$\sigma_{ij}^\infty$$
- ansatz: $$(u_1^\infty, u_2^\infty)^T\ =\ (a_{11}x_1\ +\ a_{12}x_2,\ a_{21}x_1\ +\ a_{22}x_2)^T$$
- stress-strain relation: $\Rightarrow\ \ \ \lambda a_{kk}\delta_{ij}
\ +\ 
G\left(a_{ji}+a_{ij})\right)
\ =\ \sigma^\infty_{ij}\ \ \ $(linear system for $$a_{ij}$$)

## Installation
The code is implemented in `FDBF_Code.ipynb`. For simpler visualization it is a jupyter notebook.

## Usage

1. **Modify parameters** at the top of the code for computation or reading of the basis function:
   - `solverType` (0 for Jacobi, 1 for sparse solver)
   - `readBF` (0 to compute basis functions, 1 to load from `Data/`)
   - Material and grid settings of basis function (`L`, `nhalf`, `Lame1`, `Lame2`)
   - Note that the basis function only has to be computed once for specific material and grid settings. They should be saved and used thereafter.

2. **Modify parameters** specific for test case:
   - `L_Domain`: domain size of test case
   - `n_Domain`: grid size for rendering 
   - `n_fractures`: number of fractures 
   - `readfromFile`: 1: n_fractures is adapted to number of fractures in the text file, 0: fracture need to be defined below.
   - `s11_inf`,`s22_inf`,`s12_inf`: far field stresses
   - `mu`: friction coefficient (values between 0 and 1)
   - `pressure`: fluid pressure in fractures
   - `showUFarField`: (0 without, 1 with far field displacement)


## Data Files

- `FDBF_Code.ipynb`: Main code, Jupyter Notebook
- `Data/*.csv`: Example of precomputed shear and opening basis functions for Lame1=6.3 GPa and Lame2=6.9 GPa and basis function sizes of L=8.
- `Data/Odling33.txt`: Example of fracture pattern testfile with format [Length, mid_x, mid_y, angle] (Example of paper).

## Dependencies

- Python 3.9.6
- numpy  1.25.2
- scipy 1.13.1
- matplotlib 3.7.2
- ipywidgets 8.1.0 (for interactive visualtization in Jupyter Notebook)
- IPython 8.14.0 (for displaying in Jupyter Notebook)









