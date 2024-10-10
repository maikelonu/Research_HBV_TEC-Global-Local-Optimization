# Research_HBV_TEC-Global-Local-Optimization
# HBV_OPTI

HBV_OPI is an R program designed to perform optimization and sensitivity analysis of hydrological models. To do so, it uses several local optimization methods, including:
(1) Nelder-Mead, (2) Broyden-Fletcher-Goldfarb-Shanno, (3) Hooke-Jeeves, (4) Variable Nonlinear Minimization, (5) Bound Optimization BY Quadratic Approximation, (6) Spectral Projected Gradient, (7) PORT Gradient Algorithm and (8) Levenberg-Marquardt Algorithm, as well as global optimization methods, including: (1) Generalized Simulated Annealing, (2) Differential Evolution Optimization, (3) Genetic Algorithms, (4) Shuffled Complex Evolution, (5) Enhanced Particle Swarm Optimization, (6) DIviding RECTangles for Global Optimization, (7) Controlled Random Search Local Mutation and (8) Augmented Lagrangian Minimization Algorithm.

## Installation

Use Git to clone and install the program

## Usage

Run HBV_OPTI.R script accordingly

## Contributing

Maikel Mendez Morales. Escuela de Ingeniería en Construcción, Instituto Tecnológico de Costa Rica. email: mamendez@itcr.ac.cr

Luis Alexander Calvo Valverde. Escuela de Ingeniería en Computación, Instituto Tecnológico de Costa Rica. email: lcalvo@itcr.ac.cr

## Publications

Mendez, M.; Calvo-Valverde, L.A. Comparison of global and local optimization methods for the calibration and sensitivity analysis of a conceptual hydrological model. TM 2019.
https://doi.org/10.18845/tm.v32i3.4477

![alt test](/edp01.png)

Graphical Abstract

![alt test](/FIG_TM_02.png)

Abstract: 

Eight global and eight local optimization methods were used to calibrate the HBV-TEC hydrological model on the upper Toro river catchment in Costa Rica for four different calibration periods (4, 8, 12 and 16 years). To evaluate their sensitivity to getting trapped in local minima, each method was tested against 50 sets of randomly-generated initial model parameters. All methods were then evaluated in terms of optimization performance and computational cost. Results show a comparable performance among various global and local methods as they highly
correlate to one another. Nonetheless, local methods are in general more sensitive to getting trapped in local minima, irrespective of the duration of the calibration period. Performance of the various methods seems to be independent to the total number of model calls, which may vary several orders of magnitude depending on the selected optimization method. The selection of an optimization method is largely influenced by its efficiency and the available computational resources regardless of global or local class.
Results

![alt test](/Opti02.png)

![alt test](/Opti03.png)

![alt test](/Opti04.png)

## License

[MIT](https://choosealicense.com/licenses/mit/)
