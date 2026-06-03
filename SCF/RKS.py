"""
KS-RSCF Solver for DFT calculations using Psi4
"""

import time
import numpy as np
import psi4

def diag(diag, A, nel):
    """
    Diagonalizes the Fock matrix and builds the density matrix from N/2 lowest eigenvectors.
    """

    Fp = psi4.core.triplet(A, diag, A, True, False, True)
    nbf = A.shape[0]
    Cp = psi4.core.Matrix(nbf, nbf)
    eigvals = psi4.core.Vector(nbf)
    Fp.diagonalize(Cp, eigvals, psi4.core.DiagonalizeOrder.Ascending)

    C = psi4.core.doublet(A, Cp, False, False)
    Cocc = psi4.core.Matrix(nbf, 1)
    Cocc.np[:] = C.np[:, :nel//2]

    D = psi4.core.doublet(Cocc, Cocc, False, True)  
    
    return D, eigvals.np[-1]

def Vpot_init(build_superfunctional, wfn, alias, vname, restricted=True):
    """
    Initializes a Psi4 VBase potential object
    """

    sup = build_superfunctional(alias, restricted)[0]
    sup.set_deriv(1)
    sup.allocate()
    Vpot = psi4.core.VBase.build(wfn.basisset(), sup, vname)
    
    return Vpot

def Vpot_builder(Vpot, D, V, D_half):
    """
    Computes the potential on the grid for a given density
    """

    D_half.copy(D)
    D_half.scale(0.5)
    Vpot.set_D([ D_half ])
    Vpot.compute_V([ V ])
    e = Vpot.quadrature_values()['FUNCTIONAL']
    
    return e, V

def ks_solver(mol, EXC):
    """
    Main KS-RSCF Solver Loop
    """
    
    ## Convergence thresholds
    E_conv = 1.0e-8
    D_conv = 1.0e-8
    maxiter = 20
    
    ## Wavefunction & Basis Setup
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = mints.ao_overlap()
    nbf = wfn.nso()
    nel = wfn.nalpha() + wfn.nbeta()

    print(f'Number of basis functions:   {nbf}')

    ## Potential Initialization
    build_superfunctional = psi4.driver.dft.build_superfunctional
    D_half = psi4.core.Matrix(nbf, nbf)

    VXCpot = Vpot_init(build_superfunctional, wfn, EXC, "RV", restricted=True)
    VXCpot.initialize()
    VXC_null = psi4.core.Matrix(nbf, nbf)

    ## Calculate and store V, T, H_core, ERI (I), and diagonalization matrix (A)
    V = mints.ao_potential()
    T = mints.ao_kinetic()
    H = T.clone()
    H.add(V)
    I = np.asarray(mints.ao_eri())
    A = mints.ao_overlap()
    A.power(-0.5, 1.e-14)

    ## Initialize necessary matrices
    F = psi4.core.Matrix(nbf, nbf)
    J = psi4.core.Matrix(nbf, nbf)
    Vxc = psi4.core.Matrix(nbf, nbf)
    D_diff = psi4.core.Matrix(nbf, nbf)

    D, eigvals = diag(H, A, nel)
    
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    
    print('\nStarting SCF iterations:')
    print("\n    Iter               Energy         epsilon       Delta E         dRMS\n")
    t = time.time()

    for SCF_ITER in range(1, maxiter + 1):
        D_old = D
        
        ## Build J (Coulomb)
        J_np = np.einsum('pqrs,rs->pq', I, D.np, optimize=True)
        J.np[:] = J_np

        ## Build F = H + J
        F.copy(H)
        F.axpy(1.0, J)

        ## Build DFT potentials and calculate corresponding energies
        exc, Vxc = Vpot_builder(VXCpot, D, VXC_null, D_half)

        ## Add DFT potentials to Fock matrix
        F.axpy(1.0, Vxc)

        ## Energy calculation
        SCF_E = H.vector_dot(D)
        SCF_E += 0.5 * J.vector_dot(D)
        SCF_E += exc
        SCF_E += Enuc

        ## Diagonalize Fock matrix.
        D, mu = diag(F, A, nel)

        ## Density convergence check
        D_diff.copy(D)
        D_diff.subtract(D_old)
        dRMS = D_diff.rms()
        print('SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E'
            % (SCF_ITER, SCF_E, mu, (SCF_E - Eold), dRMS))
        
        if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
            break

        Eold = SCF_E
        
        if SCF_ITER == maxiter:
            SCF_D = D
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, SCF_D, mu, SCF_ITER
    
    print('\nTotal time for SCF iterations: %.3f seconds ' % (time.time() - t))

    return SCF_E, D, mu, SCF_ITER