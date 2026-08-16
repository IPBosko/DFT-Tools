"""
Spin-unrestricted HF SCF solver using Psi4 objects for arbitrary spin fermions
"""

import sys
import psi4
import time
import numpy as np
sys.path.append('/Users/ivanbosko/Documents/CODES/GIT/DFT-Tools/SCF')
import scf_helper 

def scf_solver(mol, spin=0.5, damp=0.0, DIIS=True, verbose=False, frac=None, sym_break=False):
    """
    Spin-unrestricted HF SCF solver loop
    """
    
    ## Convergence thresholds
    E_conv = 1.0e-8
    D_conv = 1.0e-8
    maxiter = 220
    
    # Damping settings
    current_damp = damp
    if DIIS:
        damping_switch_off = 1.0e6 * D_conv
    else:
        damping_switch_off = 1.0e2 * D_conv
    
    ## Wavefunction & Basis Setup
    wfn, mints, nbf, nel, nalpha, nbeta = scf_helper.scf_main_objects(mol)

    # Multiplicity and spin channels (M = 2s + 1)
    M = int(2 * spin + 1)
    
    if M == 2:
        n_occ = [nalpha, nbeta]
    else:
        # Distribute fermions as evenly as possible across spin channels
        n_occ = [nel // M + (1 if i < nel % M else 0) for i in range(M)]

    ## Optional output
    basic_info_str = f'Number of basis functions: {nbf}\nNumber of electrons: {nel}\nSpin: {spin} (Spin channels: {M})\nOccupations: {n_occ}\n'
    psi4.core.print_out(basic_info_str)
    if verbose:
        print(basic_info_str)

    ## Calculate and store V, T, H (core), ERI (I), and diagonalization matrix (A)
    H, I, S, A = scf_helper.scf_building_blocks(mints)[2:]

    ## Initialize necessary matrices dynamically based on spin multiplicity
    D_spins = [psi4.core.Matrix(nbf, nbf) for _ in range(M)]
    D_old_spins = [psi4.core.Matrix(nbf, nbf) for _ in range(M)]
    F_spins = [psi4.core.Matrix(nbf, nbf) for _ in range(M)]
    J_spins = [psi4.core.Matrix(nbf, nbf) for _ in range(M)]
    K_spins = [psi4.core.Matrix(nbf, nbf) for _ in range(M)]
    
    D_tot = psi4.core.Matrix(nbf, nbf)
    J_tot = psi4.core.Matrix(nbf, nbf)
    D_diff = psi4.core.Matrix(nbf, nbf)

    ## Initial guess (core Hamiltonian) density matrices
    for i in range(M):
        D_spins[i] = scf_helper.diag(H, A, n_occ[i], restricted=False, frac=frac if i == 0 else None)

    ## Apply symmetry breaking if requested
    if sym_break and n_occ[0] > 0 and n_occ[0] < nbf:
        # Re-diagonalize Core Hamiltonian to get MO coefficients
        H_np = np.asarray(H)
        A_np = np.asarray(A)
        Fp = A_np.T @ H_np @ A_np
        eps, C_p = np.linalg.eigh(Fp)
        C = A_np @ C_p
        
        # Break symmetry for the first M // 2 channels 
        # (e.g., alpha for spin-1/2, alpha and beta for spin-3/2)
        for i in range(M // 2):
            if n_occ[i] > 0 and n_occ[i] < nbf:
                homo_idx = n_occ[i] - 1
                lumo_idx = n_occ[i]
                
                C_broken = C.copy()
                C_broken[:, homo_idx] = np.sqrt(0.75) * C[:, homo_idx] + np.sqrt(0.25) * C[:, lumo_idx]
                
                # Rebuild density matrix for channel i
                C_occ = C_broken[:, :n_occ[i]]
                D_broken_np = C_occ @ C_occ.T
                
                # Set the modified density matrix back to the Psi4 object
                D_spins[i].np[:] = D_broken_np

    if DIIS:
        ## Initialize diis object
        diis_objs = [psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest") for _ in range(M)]
    
    ## Nuclear energy
    Enuc = mol.nuclear_repulsion_energy()
    ## Initialize Eold for convegence check
    Eold = 0.0
    
    ## Optional output
    header = "\n    Iter               Energy            Delta E         dRMS\n"
    psi4.core.print_out(header)
    if verbose:
        print('\nStarting SCF iterations:')
        print(header)
    t = time.time()

    dRMS_spins = [0.0] * M

    ## SCF iterative procedure
    for SCF_ITER in range(1, maxiter + 1):
        
        ## Saving the initial density matrix for SCF damping
        for i in range(M):
            D_old_spins[i].copy(D_spins[i])
        
        ## Build J (Coulomb) and K (Fock exchange)
        for i in range(M):
            J_spins[i] = scf_helper.Jbuild(I, D_spins[i], J_spins[i])
            K_spins[i] = scf_helper.Kbuild(I, D_spins[i], K_spins[i])

        # Total J matrix (Coulomb interaction evaluates full density matrix)
        J_tot.copy(J_spins[0])
        for i in range(1, M):
            J_tot.axpy(1.0, J_spins[i])

        ## Build Fock matrices
        for i in range(M):
            F_spins[i].copy(H)
            F_spins[i].axpy(1.0, J_tot)
            F_spins[i].axpy(-1.0, K_spins[i])

        dRMS_tot = 0.0
        if DIIS:
            ## Build DIIS vector
            for i in range(M):
                diis_e, dRMS = scf_helper.diis_vector(F_spins[i], D_spins[i], S, A)
                diis_objs[i].add(F_spins[i], diis_e)
                dRMS_spins[i] = dRMS
                dRMS_tot += dRMS

        ## Energy calculation
        D_tot.copy(D_spins[0])
        for i in range(1, M):
            D_tot.axpy(1.0, D_spins[i])
            
        SCF_E = Enuc
        SCF_E += 0.5 * J_tot.vector_dot(D_tot)
        for i in range(M):
            SCF_E += H.vector_dot(D_spins[i])
            SCF_E -= 0.5 * K_spins[i].vector_dot(D_spins[i])

        ## DIIS convergence test and Fock extrapolation
        if DIIS:
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS_tot)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            if (abs(SCF_E - Eold) < E_conv and dRMS_tot < D_conv):
                break
            
            Eold = SCF_E
            for i in range(M):
                F_spins[i] = diis_objs[i].extrapolate()
        
        ## Diagonalize Fock matrix
        for i in range(M):
            D_spins[i] = scf_helper.diag(F_spins[i], A, n_occ[i], restricted=False, frac=frac if i == 0 else None)

        ## Density convergence check
        if not DIIS:
            dRMS_tot = 0.0
            for i in range(M):
                drms = scf_helper.density_RMS(D_diff, D_spins[i], D_old_spins[i])
                dRMS_spins[i] = drms
                dRMS_tot += drms
        
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS_tot)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            if (abs(SCF_E - Eold) < E_conv and dRMS_tot < D_conv):
                break

            Eold = SCF_E

        ## Dynamic damping
        for i in range(M):
            D_spins[i], current_damp = scf_helper.dynamic_damping(D_spins[i], D_old_spins[i], dRMS_spins[i], damp, damping_switch_off, current_damp)
        
        if SCF_ITER == maxiter:
            psi4.core.print_out("\nWARNING ! SCF did not converge. The final values are printed\n")
            print("\nWARNING ! SCF did not converge. The final values are printed")
            if M == 2:
                return SCF_E, D_spins[0], D_spins[1], SCF_ITER
            else:
                return SCF_E, D_spins, SCF_ITER
        
    psi4.core.print_out('\nSCF converged in %d iterations and %.3f seconds \n' % (SCF_ITER, time.time() - t))
    if verbose:
        print('\nSCF converged in %d iterations and %.3f seconds' % (SCF_ITER, time.time() - t))

    # Standard unpacking returns Da, Db for compatibility, or D_spins list for arbitrary spin logic
    if M == 2:
        return SCF_E, D_spins[0], D_spins[1], SCF_ITER
    else:
        return SCF_E, D_spins, SCF_ITER