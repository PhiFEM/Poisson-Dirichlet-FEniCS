# Phi-FEM for Poisson Dirichlet

Repository containing the codes used in [Duprez, Michel, and Alexei Lozinski. "ϕ-FEM: a finite element method on domains defined by level-sets." SIAM Journal on Numerical Analysis 58.2 (2020): 1008-1028.](https://hal.science/hal-03685445/document)

To run the python files you need to install [*FEniCS*](https://fenicsproject.org/). 
The easiest way is to use anaconda with the following commands: 

```
conda create -n fenics -c conda-forge fenics mshr
conda activate fenics 
pip install -U matplotlib sympy numpy 
```
