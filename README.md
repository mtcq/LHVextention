# LHVextention 
## Code to accompany: [Nonlocality activation in a photonic quantum network](https://arxiv.org/abs/2309.06501)

#### Luis Villegas-Aguilar, Emanuele Polino, Farzad Ghafari, Marco Túlio Quintino, Kiarn Laverick, Ian R. Berkman, Sven Rogge, Lynden K. Shalm, Nora Tischler, Eric G. Cavalcanti, Sergei Slussarenko, Geoff J. Pryde

This is a repository for the article [Nonlocality activation in a photonic quantum network](https://arxiv.org/abs/2309.06501).

 This code requires:
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory.
- [cvx](http://cvxr.com/) - a free MATLAB toolbox for rapid prototyping of optimization problems.

The code consists in a single script considers that constructs a Local Hidden Variable (LHV) model (for dichotomic measurements) for the quantum state experimentally obtained, a state which is close to a qubit isotropic state with visibility alpha_target=0.6373.
The algorithm builds up on the existence of previous non-trivial LHV models. In particular, we make use of the existence of a LHV model for the two-qubit Werner state presented in [https://arxiv.org/abs/2302.04721](https://arxiv.org/abs/2302.04721).

This code may be easily adapted to ensure the existence of LHV models for various other states.
