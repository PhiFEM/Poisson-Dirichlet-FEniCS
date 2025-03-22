# $\varphi$-FEM for Poisson Dirichlet

Repository containing the codes used in [Duprez, Michel, and Alexei Lozinski. "ϕ-FEM: a finite element method on domains defined by level-sets." SIAM Journal on Numerical Analysis 58.2 (2020): 1008-1028.](https://hal.science/hal-03685445/document)

## This repository is for reproducibility purposes only

It is "frozen in time" and not maintained.
To use our latest $\varphi$-FEM code please refer to the [phiFEM repository](https://github.com/PhiFEM/Poisson-Dirichlet-fenicsx).

## Usage

### Prerequisites

- [Git](https://git-scm.com/)
- [Docker](https://www.docker.com/)/[podman](https://podman.io/) (or any container engine with similar API as docker and podman.)

The image is based on the legacy FEniCS image: quay.io/fenicsproject/stable.
It contains also the `numpy`, `sympy` and `matplotlib` python libraries.

### Installation

Clone this repository
```bash
mkdir phifem-poisson-dirichlet
git clone https://github.com/PhiFEM/Poisson-Dirichlet-FEniCS.git phifem-poisson-dirichlet
cd phifem-poisson-dirichlet
```

Specify your container engine in the environment variable:
```bash
export CONTAINER_ENGINE=docker
```

Pull the container image from the `docker.io` registry:
```bash
sudo -E bash pull-image.sh
```

### Usage

Launch the image:
```bash
sudo -E bash run-image.sh
```

From the container, launch a test case, for example:
```bash
python3 phiFEM_test_case1.py
```

## Issues and support

Please use the issue tracker to report any issues.


## Authors (alphabetical)

[Michel Duprez](https://michelduprez.fr/), Inria Nancy Grand-Est  
[Alexei Lozinski](https://orcid.org/0000-0003-0745-0365), Université de Franche-Comté
