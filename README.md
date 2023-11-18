# Polynomials of Random Matrices

Torben Krueger and Yuriy Nemish, 2023.

### Repository contains:
- `PolynomialsOfRandomMatrices.wl` - Mathematica package for numerical simulations of the eigenvalue distributions for polynomials of random matrices and related models
- `PolynomialsOfRandomMatrices_Examples.nb` - Mathematica notebook with typical simple examples of applications of functions from `PolynomialsOfRandomMatrices.wl`
- `Kronecker_Examples.nb` - Mathematica notebook with examples of applications of `PolynomialsOfRandomMatrices.wl` to Kronecker random matrices, in particular
  - comparing the empirical spectral density and the predicted (theoretical) self-adjoint pseudospecturum of non-Hermitian Kronecker random matrices
  - simulating the empirical spectral density of non-Hermitian Kronecker random matrices with nonnormal expectations
- `Polynomials_Examples.nb` - Mathematica notebook with examples of applications of `PolynomialsOfRandomMatrices.wl` to polynomials of random (Wigner) matrices, in particular
  - computing the limiting eigenvalue distribution of noncommutative polynomials in Wigner matrices by solving the Dyson equation for linearizations
  - numerically verifying the conjecture about the nonvanishing diagonal of the imaginary part of the solution to the Dyson equation for the minimal linearization in the bulk
  - numerically comparing the smallest eigenvalues of the imaginary part of the solution to the Dyson equation for the big and the minimal linearizations
- `QuantumDot_Examples.nb` - Mathematica notebook with examples of applications of `PolynomialsOfRandomMatrices.wl` to the study of the random matrix model of transport in quantum dots, in particular
  - computing the self-consistent density of the transport eigenvalues for various energy levels
  - computing the Fano factor for various energy levels

For each Mathematica notebook the repository also contains a .pdf file with the output of the corresponding .nb file.

Folders `ev_40K_interval` and `fano_phi_50_9k` contain auxiliary data files used in Mathematica notebooks.

### Usage
To use the package, load the package in your Mathematica notebook together with [NCAlgebra package](https://github.com/NCAlgebra).

### Literature
- `Kronecker_Examples.nb`: Alt, Erdos, Krueger and Nemish, 2019. Location of the spectrum of Kronecker random matrices [journal](https://doi.org/10.1214/18-aihp894), [arXiv](http://arxiv.org/abs/1706.08343)
- `Polynomials_Examples.nb`: Erdos, Krueger and Nemish, 2020. Local laws for polynomials of Wigner matrices [journal](https://doi.org/10.1016/j.jfa.2020.108507), [arXiv](http://arxiv.org/abs/1804.11340)
- `QuantumDot_Examples.nb`: Erdos, Krueger and Nemish, 2021. Scattering in quantum dots via noncommutative rational functions [journal](https://link.springer.com/article/10.1007/s00023-021-01085-6), [arXiv](https://arxiv.org/abs/1911.05112)
