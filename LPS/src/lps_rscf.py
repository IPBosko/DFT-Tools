"""
LPS-RSCF: Levi-Perdew-Sahni Restricted Self-Consistent Field Solver
Implementation of the LPS DFT method using Psi4 KS DFT code.
See LPS-RSCF.ipynb for details.
"""

import time
import numpy as np
import psi4
import os
import sys
project_root = os.path.abspath('..')
if project_root not in sys.path:
    sys.path.append(project_root)
from src.build_density import diag_lps
from src.build_Vpot import Vpot_init, Vpot_builder

def lps_solver(mol, E_conv, D_conv, maxiter, TP, lam, EXC, FA, damp, DIIS, Guess=None, verbose=True):
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
        DIIS (list):        Enable DIIS extrapolation.
        verbose (bool):     Print iteration details.
        
    Returns:
        tuple: (SCF_E, SCF_D, SCF_ITER)
    """
    
    psi4.core.set_output_file('output.dat', False)
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    S = mints.ao_overlap()
    nbf = wfn.nso()
    nel = wfn.nalpha() + wfn.nbeta()

    if verbose:
        print(f'Number of basis functions:   {nbf}')

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

    vW = {
        "name": "vW",
        "x_functionals": {"LDA_X": {"alpha": 0.00}},
        "c_functionals": {"GGA_K_VW": {"alpha": 1.00}}
    }
    VvWpot = Vpot_init(build_superfunctional, wfn, vW, "RV", restricted=True)
    VvWpot.initialize()
    VvW_null = psi4.core.Matrix(nbf, nbf)

    V = mints.ao_potential()
    T = mints.ao_kinetic()
    H = T.clone()
    H.add(V)
    I = np.asarray(mints.ao_eri())
    A = mints.ao_overlap()
    A.power(-0.5, 1.e-14)

    F = psi4.core.Matrix(nbf, nbf)
    J = psi4.core.Matrix(nbf, nbf)
    VG = psi4.core.Matrix(nbf, nbf)
    
    if Guess is not None:
        D = Guess.clone()
        mu = 0.0
    else:
        D, mu = diag_lps(H, A, nel)
    
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    
    header = "\n    Iter               Energy         ChemPot       Delta E         dRMS\n"
    psi4.core.print_out(header)
    if verbose:
        print('\nStarting SCF iterations:')
        print(header)
    t = time.time()
    e_list = []
    e_conv_list = []
    d_conv_list = []

    diis_obj = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")

    for SCF_ITER in range(1, maxiter + 1):
        
        D_old = D
        J_np = np.einsum('pqrs,rs->pq', I, D.np, optimize=True)
        J.np[:] = J_np
        F.copy(H)
        F.axpy(1.0, J)
        if FA[0]:
            if nel == 0:
                F.axpy(0.0, J)
            else: 
                F.axpy(-FA[1]/nel, J)

        pau_e, VP = Vpot_builder(VPpot, D, VP_null, D_half)
        xc_e, VXC = Vpot_builder(VXCpot, D, VXC_null, D_half)
        vw_e, VvW = Vpot_builder(VvWpot, D, VvW_null, D_half)
        g_e = pau_e + xc_e + ( lam - 1.0 ) * vw_e 

        VG.copy(VP)
        VG.axpy(1.0, VXC)
        VG.axpy((lam - 1.0), VvW)
        F.axpy(1.0, VG)

        diis_e = psi4.core.triplet(F, D, S, False, False, False)
        diis_e.subtract(psi4.core.triplet(S, D, F, False, False, False))
        diis_e = psi4.core.triplet(A, diis_e, A, False, False, False)
        dRMS = diis_e.rms()
        d_conv_list.append(np.log10(dRMS))

        SCF_E = H.vector_dot(D)
        SCF_E += 0.5 * J.vector_dot(D)
        if FA[0]:
            SCF_E += 0.5 * J.vector_dot(D) * ( - FA[1] / nel )
        SCF_E += g_e
        SCF_E += Enuc
        e_conv_list.append(np.log10(abs(SCF_E - Eold)))
        e_list.append(SCF_E)

        output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, mu, (SCF_E - Eold), dRMS)
        psi4.core.print_out(output_str)
        if verbose:
            print(output_str.strip())
        if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
            break
        
        Eold = SCF_E

        if (dRMS > damp[2]):
            current_damp = damp[0]
            diis_active = DIIS[0]
        else:
            current_damp = damp[1]
            diis_active = DIIS[1]

        if diis_active:
            diis_obj.add(F, diis_e)
            F = diis_obj.extrapolate()
        
        D, mu = diag_lps(F, A, nel)
        D.scale(1.0 - current_damp)
        D.axpy(current_damp, D_old)
        
        if SCF_ITER == maxiter:
            SCF_D = D
            psi4.core.print_out("\nWARNING ! SCF did not converge. The final values are printed\n")
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, SCF_D, mu, SCF_ITER, e_list, e_conv_list, d_conv_list
    
    psi4.core.print_out('\nTotal time for SCF iterations: %.3f seconds \n' % (time.time() - t))
    if verbose:
        print('\nTotal time for SCF iterations: %.3f seconds ' % (time.time() - t))

    return SCF_E, D, mu, SCF_ITER, e_list, e_conv_list, d_conv_list

# D_diff.copy(D)
# D_diff.subtract(D_old)
# dRMS = D_diff.rms()
# d_conv_list.append(np.log10(dRMS))