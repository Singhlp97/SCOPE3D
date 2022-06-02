new_recon.f: main program

in.com: 
- input parameter file
- dt: time step
- nstep: number time steps 
- noutt: output sampling rate for spatially integrated quantities (ex. energy)
- noutf: output sampling rate for fields
- xl, yl, zl: grid parameters
- t0,gridalfa,gridbeta,gridgamma: grid parameters for not-equispaced grid along x
- de, beta_e, tau, eta, mm, pso, psoeq, phieq, eq_l, asym: physical parameters
- istart: restart parameter  
- omegax, omega, alfa1: filtering parameters for numerical noise control

par.inc: 
- setup number of points along x, y and z  (input)
- contains the definitions of the COMMON arrays

pstartup.f
- setup number of cores along z, nprocz (input)

init_mod.f 
- grid setup along x, y and z
- fft initialization
- setup filter parameters 
- setup tridiagonal matrices and LU factorization for filter, first and second derivative along x and for Poisson and Helmoltz operators   

condinit.f:  
- compute the initial conditions

rhs_xy.f:
- compute the terms on the right-hand side of the equations containing the spatial derivatives along x and y 

rhs_z.f:
- compute the terms on the right-hand side of the equations containing the spatial derivatives along z 

poisson_mod.f:
- solve the Poisson problem 
nabla^2 phi = U (U input, phi output)

helm_mod.f:
- solve the Helmholtz problem 
nabla^2 psi + (1/de^2) psi = F (F input, psi output)

out_field.f
- output file for the 3D fields (es. Psi.dat, Phi.dat,...) 

out_inv.f
- output file for energy components and invariant quantities (energuie.dat )


Directory PRO
- post process the simulation results. In particular,
  gamma_eta.pro 
    - tests the accuracy of the simulation resulting from the combination of the KH mode and resistive tearing mode
  gamma_de.pro
    - tests the accuracy of the simulation resulting from the combination of the KH mode and collisionless tearing mode




