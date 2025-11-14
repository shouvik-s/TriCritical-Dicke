# Tri-Critical Dicke Model (TCDM)
_[Based on PRL 122, 193201 (2019) (1904.10576)]_
_[Reduces to Dicke Model at \epsilon = 0]_

This repository contains Julia codes to obtain 
1. The phase diagram of the TCDM.
2. The quadratures that get perfectly squeezed across the superradiant quantum phase transition (SRPT). 
3. Computations of two-parameter phase estimation. 

All physical parameters and scan ranges are controlled by a **Makefile**.


## Details of the three components
1. This part is parallelized by using a **distributed implementation of DMRG** with the help of the `ITensors`/`ITensorsMPS` and `Distributed` libraries.
2. Uses the output from **1** and applies a 3-parameter optimization to obtain the quadrature that is most squeezed; uses the `Optim` package.