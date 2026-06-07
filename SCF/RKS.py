"""
Spin-restricted KS SCF solver for DFT calculations using Psi4 objects
"""

import sys
import psi4
import time
import numpy as np
sys.path.append('/Users/ivanbosko/Documents/CODES/GIT/DFT-Tools/SCF')
import scf_helper 

def ks_solver(mol, EXC, damp, DIIS=True, verbose=False):
    """
    Spin-restricted KS SCF solver loop with Hartree damping
    """
    
    ## Convergence thresholds
    E_conv = 1.0e-8
    D_conv = 1.0e-8
    maxiter = 40
    
    current_damp = damp
    damping_switch_off = 1.0e1 * D_conv
    
    ## Wavefunction & Basis Setup
    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option("BASIS"))
    mints = psi4.core.MintsHelper(wfn.basisset())
    nbf = wfn.nso()
    nalpha = wfn.nalpha()

    ## Optional output
    basic_info_str = f'Number of basis functions: {nbf}\nNumber of electrons: {2*nalpha}\nNumber of occupied orbitals: {nalpha}\n'
    psi4.core.print_out(basic_info_str)
    if verbose:
        print(basic_info_str)

    ## Potential Initialization
    build_superfunctional = psi4.driver.dft.build_superfunctional
    D_half = psi4.core.Matrix(nbf, nbf)

    VXCpot = scf_helper.Vpot_init(build_superfunctional, wfn, EXC, "RV", restricted=True)
    VXCpot.initialize()
    VXC_null = psi4.core.Matrix(nbf, nbf)

    ## Calculate and store V, T, H (core), ERI (I), and diagonalization matrix (A)
    V = mints.ao_potential()
    T = mints.ao_kinetic()
    H = T.clone()
    H.add(V)
    I = np.asarray(mints.ao_eri())
    S = mints.ao_overlap()
    A = S.clone()
    A.power(-0.5, 1.e-14)

    ## Initialize necessary matrices
    F = psi4.core.Matrix(nbf, nbf)
    J = psi4.core.Matrix(nbf, nbf)
    K = psi4.core.Matrix(nbf, nbf)
    Vxc = psi4.core.Matrix(nbf, nbf)
    D_diff = psi4.core.Matrix(nbf, nbf)

    ## Initial guess (core Hamiltonian) density matrix
    D, homo = scf_helper.diag(H, A, nalpha)

    ## Initialize diis object
    if DIIS:
        diis_obj = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")
    
    Enuc = mol.nuclear_repulsion_energy()
    Eold = 0.0
    
    header = "\n    Iter               Energy         HOMO       Delta E         dRMS\n"
    psi4.core.print_out(header)
    if verbose:
        print('\nStarting SCF iterations:')
        print(header)
    t = time.time()

    for SCF_ITER in range(1, maxiter + 1):
        
        ## Saving the initial density matrix for SCF damping
        D_old = D
        
        ## Build J (Coulomb)
        J_np = np.einsum('pqrs,rs->pq', I, D.np, optimize=True)
        J.np[:] = J_np
        
        ## Build K (Fock exchange)
        if EXC["name"]=="EXX":
            K_np = np.einsum('prqs,rs->pq', I, D.np, optimize=True)
            K.np[:] = K_np

        ## Build F = H + J
        F.copy(H)
        F.axpy(1.0, J)
        if EXC["name"]=="EXX":
            F.axpy(-0.5, K)

        ## Build DFT potentials and calculate corresponding energies
        if EXC["name"]!="EXX":
            
            exc, Vxc = scf_helper.Vpot_builder(VXCpot, D, VXC_null, D_half)

            ## Add DFT potentials to Fock matrix
            F.axpy(1.0, Vxc)

        if DIIS:
            diis_e = psi4.core.triplet(F, D, S, False, False, False)
            diis_e.subtract(psi4.core.triplet(S, D, F, False, False, False))
            diis_e = psi4.core.triplet(A, diis_e, A, False, False, False)
        
            diis_obj.add(F, diis_e)
            dRMS = diis_e.rms()

        ## Energy calculation
        SCF_E = H.vector_dot(D)
        SCF_E += 0.5 * J.vector_dot(D)
        if EXC["name"]=="EXX":
            SCF_E -= 0.25 * K.vector_dot(D)
        else:
            SCF_E += exc
        SCF_E += Enuc

        ## DIIS convergence test and Fock extrapolation
        if DIIS:
            if verbose:
                print('SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E'
                    % (SCF_ITER, SCF_E, homo, (SCF_E - Eold), dRMS))
            
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break
            
            Eold = SCF_E
            F = diis_obj.extrapolate()
        
        ## Diagonalize Fock matrix.
        D, homo = scf_helper.diag(F, A, nalpha)

        ## Density convergence check
        if not DIIS:
            D_diff.copy(D)
            D_diff.subtract(D_old)
            dRMS = D_diff.rms()
        
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, homo, (SCF_E - Eold), dRMS)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break

            Eold = SCF_E

        ## Dynamic damping
        if dRMS < damping_switch_off:
            current_damp *= 0
        else:
            current_damp = damp 
        D.scale(1.0 - current_damp)
        D.axpy(current_damp, D_old)
        
        if SCF_ITER == maxiter:
            
            psi4.core.print_out("\nWARNING ! SCF did not converge. The final values are printed\n")
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, D, homo, SCF_ITER
        
    psi4.core.print_out('\nSCF converged in %d iterations and %.3f seconds \n' % (SCF_ITER, time.time() - t))
    print('\nSCF converged in %d iterations and %.3f seconds' % (SCF_ITER, time.time() - t))

    return SCF_E, D, homo, SCF_ITER