# PDDO-Upwind
PDDO-Upwind: On the Solution of Hyperbolic Equations Using the Peridynamic Differential Operator

Numerical solution of hyperbolic differential equations, such as the advection equation, poses challenges. Classically, this issue has been addressed by using a scheme known as the upwind scheme. It simply invokes more points from the upwind side of the flow stream when calculating derivatives. This study presents a generalized upwind scheme, referred to as directional nonlocality, for the numerical solution of linear and nonlinear hyperbolic Partial Differential Equations (PDEs) using the peridynamic differential operator (PDDO). The PDDO provides the nonlocal form of the differential equations by introducing an internal length parameter (horizon) and a weight function. The weight function controls the degree of interaction among the points within the horizon. A modification to the weight function, i.e., upwinded-weight function, accounts for directional nonlocality along which information travels. This modification results in a stable PDDO discretization of hyperbolic PDEs. Solutions are constructed in a consistent manner without special treatments through simple discretization. The capability of this approach is demonstrated by considering time dependent linear and nonlinear hyperbolic equations as well as the time invariant Eikonal equation. Numerical stability is ensured for the linear advection equation and the PD solutions compare well with the analytical/reference solutions.

## Citation
If you find our work useful in your research, please cite:
```
@article{BEKAR2022114574,
title = {On the solution of hyperbolic equations using the peridynamic differential operator},
journal = {Computer Methods in Applied Mechanics and Engineering},
volume = {391},
pages = {114574},
year = {2022},
issn = {0045-7825},
doi = {https://doi.org/10.1016/j.cma.2022.114574},
url = {https://www.sciencedirect.com/science/article/pii/S0045782522000032},
author = {Ali Can Bekar and Erdogan Madenci and Ehsan Haghighat},
keywords = {Peridynamics, Nonlocal, Hyperbolic, Advection, Eikonal},
}
```

## Contact
If you have any questions, please feel free to email <acbekar@email.arizona.edu>.
