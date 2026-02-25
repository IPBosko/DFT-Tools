"""
LPS-USCF: Levi-Perdew-Sahni Unrestricted Self-Consistent Field Solver
Implementation of the LPS DFT method using Psi4 KS DFT code.
See LPS-USCF.ipynb for details.
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
    return D

def Vpot_init(build_superfunctional, wfn, alias, vname, restricted=True):
    """Initializes a Psi4 VBase potential object."""
    sup = build_superfunctional(alias, restricted)[0]
    sup.set_deriv(1)
    sup.allocate()
    Vpot = psi4.core.VBase.build(wfn.basisset(), sup, vname)
    return Vpot

def Vpot_builder(Vpot, Da, Db, Va, Vb):
    """Computes the potential on the grid for a given density."""
    Vpot.set_D([ Da, Db ])
    Vpot.compute_V([ Va, Vb ])
    e = Vpot.quadrature_values()['FUNCTIONAL']
    return e, Va, Vb

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
    D_conv = 1.0e-4
    
    ## Wavefunction & Basis Setup
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = mints.ao_overlap()
    nbf = wfn.nso()
    nalpha = wfn.nalpha()
    nbeta = wfn.nbeta()

    if verbose:
        print(f'Number of basis functions:   {nbf}')

    ## Potential Initialization
    build_superfunctional = psi4.driver.dft.build_superfunctional

    Pauli = {
    "name": "Pauli",
    "x_functionals": {"LDA_X": {"alpha": 0.00}},
    "c_functionals": {TP[0]: {"alpha": TP[1]}}
    }
    VPpot = Vpot_init(build_superfunctional, wfn, Pauli, "UV", restricted=False)
    VPpot.initialize()
    VPa_null = psi4.core.Matrix(nbf, nbf)
    VPb_null = psi4.core.Matrix(nbf, nbf)

    XC = {
    "name": "XC",
    "x_functionals": {EXC[0]: {"alpha": EXC[1]}},
    "c_functionals": {EXC[2]: {"alpha": EXC[3]}}
    }   
    VXCpot = Vpot_init(build_superfunctional, wfn, XC, "UV", restricted=False)
    VXCpot.initialize()
    VXCa_null = psi4.core.Matrix(nbf, nbf)
    VXCb_null = psi4.core.Matrix(nbf, nbf)

    ## Initialize the von Weizsacker potential
    vW = {
        "name": "vW",
        "x_functionals": {"LDA_X": {"alpha": 0.00}},
        "c_functionals": {"GGA_K_VW": {"alpha": 1.00}}
    }
    VvWpot = Vpot_init(build_superfunctional, wfn, vW, "UV", restricted=False)
    VvWpot.initialize()
    VvWa_null = psi4.core.Matrix(nbf, nbf)
    VvWb_null = psi4.core.Matrix(nbf, nbf)

    ## Calculate and store V, T, H_core, ERI (I), and diagonalization matrix (A)
    V = mints.ao_potential()
    T = mints.ao_kinetic()
    H = T.clone()
    H.add(V)
    I = np.asarray(mints.ao_eri())
    A = mints.ao_overlap()
    A.power(-0.5, 1.e-14)

    ## Initialize necessary matrices
    Fa = psi4.core.Matrix(nbf, nbf)
    Fb = psi4.core.Matrix(nbf, nbf)
    J = psi4.core.Matrix(nbf, nbf)
    FAa = psi4.core.Matrix(nbf, nbf)
    FAb = psi4.core.Matrix(nbf, nbf)
    VGa = psi4.core.Matrix(nbf, nbf)
    VGb = psi4.core.Matrix(nbf, nbf)
    D = psi4.core.Matrix(nbf, nbf)
    D_diff = psi4.core.Matrix(nbf, nbf)
    
    if D_guess is not None:
        ## Use D_GUESS as an initial guess
        Da = D_guess[0].clone()
        Db = D_guess[1].clone()
    else:
        ## Calculate an initial Core guess
        Da = diag_lps(H, A, nalpha)
        Db = diag_lps(H, A, nbeta)
    
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    
    if verbose:
        print('\nStarting SCF iterations:')
        print("\n    Iter            Energy            Delta E         dRMS\n")
    t = time.time()
    e_conv_list = []
    d_conv_list = []

    diis_objA = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")
    diis_objB = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")

    for SCF_ITER in range(1, maxiter + 1):
        Da_old = Da
        Db_old = Db
        
        ## Build J (Coulomb)
        Ja = np.einsum('pqrs,rs->pq', I, Da.np, optimize=True)
        Jb = np.einsum('pqrs,rs->pq', I, Db.np, optimize=True)
        J_np = Ja + Jb
        J.np[:] = J_np

        ## Build F = H + J
        Fa.copy(H)
        Fb.copy(H)
        Fa.axpy(1.0, J)
        Fb.axpy(1.0, J)

        ## Add Fermi–Amaldi potential if requested
        if FA[0]:
            FAa.np[:] = Ja   
            FAa.np[:] *= -FA[1]/nalpha
            FAb.np[:] = Jb
            if nbeta == 0:
                FAb.np[:] *= 0.0
            else:
                FAb.np[:] *= -FA[1]/nbeta
            Fa.axpy(1.0, FAa)
            Fb.axpy(1.0, FAb)

        ## Build DFT potentials and calculate corresponding energies
        pau_e, VPa, VPb = Vpot_builder(VPpot, Da, Db, VPa_null, VPb_null)
        xc_e, VXCa, VXCb = Vpot_builder(VXCpot, Da, Db, VXCa_null, VXCb_null)
        vw_e, VvWa, VvWb = Vpot_builder(VvWpot, Da, Db, VvWa_null, VvWb_null)
        
        ## Caclulate G[n] = T_P[n] + E_xc[n] + (lam - 1)T_vW
        g_e = pau_e + xc_e + ( lam - 1.0 ) * vw_e 

        ## Add DFT potentials to Fock matrix
        VGa.copy(VPa)
        VGb.copy(VPb)
        VGa.axpy(1.0, VXCa)
        VGb.axpy(1.0, VXCb)
        VGa.axpy((lam - 1.0), VvWa)
        VGb.axpy((lam - 1.0), VvWb)
        Fa.axpy(1.0, VGa)
        Fb.axpy(1.0, VGb)

        if DIIS:
            diis_eA = psi4.core.triplet(Fa, Da, S, False, False, False)
            diis_eA.subtract(psi4.core.triplet(S, Da, Fa, False, False, False))
            diis_eA = psi4.core.triplet(A, diis_eA, A, False, False, False)
            diis_objA.add(Fa, diis_eA)
            dRMS_a = diis_eA.rms()

            diis_eB = psi4.core.triplet(Fb, Db, S, False, False, False)
            diis_eB.subtract(psi4.core.triplet(S, Db, Fb, False, False, False))
            diis_eB = psi4.core.triplet(A, diis_eB, A, False, False, False)
            diis_objB.add(Fb, diis_eB)
            dRMS_b = diis_eB.rms()

            dRMS = dRMS_a + dRMS_b

            d_conv_list.append(np.log10(dRMS))

        ## Energy calculation
        SCF_Ea = H.vector_dot(Da)
        SCF_Eb = H.vector_dot(Db)
        SCF_E = SCF_Ea + SCF_Eb
        D.copy(Da)
        D.axpy(1.0, Db)
        SCF_E += 0.5 * J.vector_dot(D)
        if FA[0]:
            SCF_E += 0.5 * FAa.vector_dot(Da)
            SCF_E += 0.5 * FAb.vector_dot(Db)
        SCF_E += g_e
        SCF_E += Enuc

        e_conv_list.append(np.log10(abs(SCF_E - Eold)))

        ##  DIIS convergence check
        if DIIS:
            if verbose:
                print('SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E'
                    % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS))
            
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break
            
            Eold = SCF_E
            Fa = diis_objA.extrapolate()
            Fb = diis_objB.extrapolate()

        ## Diagonalize Fock matrix.
        Da = diag_lps(Fa, A, nalpha)
        Db = diag_lps(Fb, A, nbeta)

        ## non-DIIS convergence check
        if not DIIS:
            D_diff.copy(Da)
            D_diff.subtract(Da_old)
            dRMS_a = D_diff.rms()
            D_diff.copy(Db)
            D_diff.subtract(Db_old)
            dRMS_b = D_diff.rms()
            dRMS = dRMS_a + dRMS_b

            d_conv_list.append(np.log10(dRMS))

            if verbose:
                print('SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E'
                    % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS))
            
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break

            Eold = SCF_E
        
        ## Dynamic damping
        # if (dRMS > damp[2]):
        if (abs(SCF_E - Eold) > damp[2]):
            current_damp = damp[0]
            DIIS = True
        else:
            current_damp = damp[1]
            DIIS = True
        Da.scale(1.0 - current_damp)
        Da.axpy(current_damp, Da_old)
        Db.scale(1.0 - current_damp)
        Db.axpy(current_damp, Db_old)
        
        if SCF_ITER == maxiter:
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, Da, Db, SCF_ITER, e_conv_list, d_conv_list
    
    if verbose:
        print('\nTotal time for SCF iterations: %.3f seconds ' % (time.time() - t))

    return SCF_E, Da, Db, SCF_ITER, e_conv_list, d_conv_list