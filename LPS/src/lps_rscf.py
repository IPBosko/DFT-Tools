"""
LPS-RSCF: Levi-Perdew-Sahni Restricted Self-Consistent Field Solver
Implementation of the LPS DFT method using Psi4 KS DFT code.
See LPS-RSCF.ipynb for details.
"""

import time
import numpy as np
import psi4

def diag_lps(diag, A, nel):
    """
    Diagonalizes the Fock matrix and builds the density matrix 
    from the lowest eigenvector only.
    """
    Fp = psi4.core.triplet(A, diag, A, True, False, True)
    nbf = A.shape[0]
    Cp = psi4.core.Matrix(nbf, nbf)
    eigvals = psi4.core.Vector(nbf)
    Fp.diagonalize(Cp, eigvals, psi4.core.DiagonalizeOrder.Ascending)

    C = psi4.core.doublet(A, Cp, False, False)
    Cocc = psi4.core.Matrix(nbf, 1)
    
    # Normalize lowest eigenvector to number of electrons
    Cocc.np[:] = np.sqrt(nel) * C.np[:, :1]

    D = psi4.core.doublet(Cocc, Cocc, False, True)  
    return D, eigvals.np[0]

def Vpot_init(build_superfunctional, wfn, alias, vname, restricted=True):
    """Initializes a Psi4 VBase potential object."""
    sup = build_superfunctional(alias, restricted)[0]
    sup.set_deriv(1)
    sup.allocate()
    Vpot = psi4.core.VBase.build(wfn.basisset(), sup, vname)
    return Vpot

def Vpot_builder(Vpot, D, V, D_half):
    """Computes the potential on the grid for a given density."""
    D_half.copy(D)
    D_half.scale(0.5)
    Vpot.set_D([ D_half ])
    Vpot.compute_V([ V ])
    e = Vpot.quadrature_values()['FUNCTIONAL']
    return e, V

def lps_solver(maxiter, TP, EXC, lam, mol, damp, FA, D_guess=None, DIIS=True, verbose=True):
    """
    Main LPS-RSCF Solver Loop.
    
    Args:
        maxiter (int):      Maximum SCF iterations.
        TP (list):          Pauli kinetic options [Functional_Name, Alpha].
        EXC (list):         XC options [X_Name, X_Alpha, C_Name, C_Alpha].
        lam (float):        Lambda parameter for T_vW subtraction.
        mol (Molecule):     Psi4 Molecule object.
        damp (list):        [damp_start, damp_end, cutoff].
        FA (list):          [Apply_FA (bool), Scaling_Factor (float)].
        D_guess (Matrix):   Optional density matrix for initial guess.
        DIIS (bool):        Enable DIIS extrapolation.
        verbose (bool):     Print iteration details.
        
    Returns:
        tuple: (SCF_E, SCF_D, SCF_ITER)
    """
    
    ## Convergence thresholds
    E_conv = 1.0e-5
    D_conv = 1.0e-5
    
    ## Wavefunction & Basis Setup
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = mints.ao_overlap()
    nbf = wfn.nso()
    nel = wfn.nalpha() + wfn.nbeta()

    if verbose:
        print(f'Number of basis functions:   {nbf}')

    ## Potential Initialization
    build_superfunctional = psi4.driver.dft.build_superfunctional
    D_half = psi4.core.Matrix(nbf, nbf)

    Pauli = {
    "name": "Pauli",
    "x_functionals": {"LDA_X": {"alpha": 0.00}},
    "c_functionals": {TP[0]: {"alpha": TP[1]}}
    }
    VPpot = Vpot_init(build_superfunctional, wfn, Pauli, "RV", restricted=True)
    VPpot.initialize()
    VP_null = psi4.core.Matrix(nbf, nbf)

    XC = {
    "name": "XC",
    "x_functionals": {EXC[0]: {"alpha": EXC[1]}},
    "c_functionals": {EXC[2]: {"alpha": EXC[3]}}
    }   
    VXCpot = Vpot_init(build_superfunctional, wfn, XC, "RV", restricted=True)
    VXCpot.initialize()
    VXC_null = psi4.core.Matrix(nbf, nbf)

    ## Initialize the von Weizsacker potential
    vW = {
        "name": "vW",
        "x_functionals": {"LDA_X": {"alpha": 0.00}},
        "c_functionals": {"GGA_K_VW": {"alpha": 1.00}}
    }
    VvWpot = Vpot_init(build_superfunctional, wfn, vW, "RV", restricted=True)
    VvWpot.initialize()
    VvW_null = psi4.core.Matrix(nbf, nbf)

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
    VG = psi4.core.Matrix(nbf, nbf)
    D_diff = psi4.core.Matrix(nbf, nbf)
    
    if D_guess is not None:
        ## Use D_GUESS as an initial guess
        D = D_guess.clone()
        mu = 0.0
    else:
        ## Calculate an initial Core guess
        D, mu = diag_lps(H, A, nel)
    
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    
    if verbose:
        print('\nStarting SCF iterations:')
        print("\n    Iter               Energy         ChemPot       Delta E         dRMS\n")
    t = time.time()

    if DIIS:
        diis_obj = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")

    for SCF_ITER in range(1, maxiter + 1):
        D_old = D
        
        ## Build J (Coulomb)
        J_np = np.einsum('pqrs,rs->pq', I, D.np, optimize=True)
        J.np[:] = J_np

        ## Build F = H + J
        F.copy(H)
        F.axpy(1.0, J)

        ## Add Fermi–Amaldi potential if requested
        if FA[0]:
            if nel == 0:
                F.axpy(0.0, J)
            else: 
                F.axpy(-FA[1]/nel, J)

        ## Build DFT potentials and calculate corresponding energies
        pau_e, VP = Vpot_builder(VPpot, D, VP_null, D_half)
        xc_e, VXC = Vpot_builder(VXCpot, D, VXC_null, D_half)
        vw_e, VvW = Vpot_builder(VvWpot, D, VvW_null, D_half)
        
        ## Caclulate G[n] = T_P[n] + E_xc[n] + (lam - 1)T_vW
        g_e = pau_e + xc_e + ( lam - 1.0 ) * vw_e 

        ## Add DFT potentials to Fock matrix
        VG.copy(VP)
        VG.axpy(1.0, VXC)
        VG.axpy((lam - 1.0), VvW)
        F.axpy(1.0, VG)

        if DIIS:
            diis_e = psi4.core.triplet(F, D, S, False, False, False)
            diis_e.subtract(psi4.core.triplet(S, D, F, False, False, False))
            diis_e = psi4.core.triplet(A, diis_e, A, False, False, False)
        
            diis_obj.add(F, diis_e)
            dRMS = diis_e.rms()

        ## Energy calculation
        SCF_E = H.vector_dot(D)
        SCF_E += 0.5 * J.vector_dot(D)
        if FA[0]:
            SCF_E += 0.5 * J.vector_dot(D) * ( - FA[1] / nel )
        SCF_E += g_e
        SCF_E += Enuc

        ##  DIIS convergence check
        if DIIS:
            if verbose:
                print('SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E'
                    % (SCF_ITER, SCF_E, mu, (SCF_E - Eold), dRMS))
            
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break
            
            Eold = SCF_E
            F = diis_obj.extrapolate()

        ## Diagonalize Fock matrix.
        D, mu = diag_lps(F, A, nel)

        ## non-DIIS convergence check
        if not DIIS:
            D_diff.copy(D)
            D_diff.subtract(D_old)
            dRMS = D_diff.rms()
            if verbose:
                print('SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E'
                    % (SCF_ITER, SCF_E, mu, (SCF_E - Eold), dRMS))
            
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break

            Eold = SCF_E
        
        ## Dynamic damping
        if (dRMS > damp[2]):
            current_damp = damp[0]
        else:
            current_damp = damp[1]
        D.scale(1.0 - current_damp)
        D.axpy(current_damp, D_old)
        
        if SCF_ITER == maxiter:
            SCF_D = D
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, SCF_D, mu, SCF_ITER
    
    if verbose:
        print('\nTotal time for SCF iterations: %.3f seconds ' % (time.time() - t))

    return SCF_E, D, mu, SCF_ITER