"""
Spin-restricted KS SCF solver for DFT calculations using Psi4 objects
"""

import sys
import psi4
import time
import numpy as np
sys.path.append('/Users/ivanbosko/Documents/CODES/GIT/DFT-Tools/SCF')
import scf_helper 

def scf_solver(mol, EXC, damp=0.0, DIIS=True, verbose=False):
    """
    Spin-restricted KS SCF solver loop
    """
    
    ## Convergence thresholds
    E_conv = 1.0e-8
    D_conv = 1.0e-8
    maxiter = 40
    
    # Damping settings
    current_damp = damp
    if DIIS:
        damping_switch_off = 1.0e8 * D_conv
    else:
        damping_switch_off = 1.0e2 * D_conv
    
    ## Wavefunction & Basis Setup
    wfn,mints,nbf,nel,nalpha,nbeta = scf_helper.scf_main_objects(mol)

    ## Optional output
    basic_info_str = f'Number of basis functions: {nbf}\nNumber of electrons: {2*nalpha}\nNumber of occupied orbitals: {nalpha}\n'
    psi4.core.print_out(basic_info_str)
    if verbose:
        print(basic_info_str)

    ## Potential Initialization
    build_superfunctional = psi4.driver.dft.build_superfunctional
    VXCpot = scf_helper.Vpot_init(build_superfunctional, wfn, EXC, "RV", restricted=True)
    VXCpot.initialize()

    ## Calculate and store V, T, H (core), ERI (I), and diagonalization matrix (A)
    H, I, S, A = scf_helper.scf_building_blocks(mints)[2:]

    ## Initialize necessary matrices
    F, J, K, Vxc, VXC_null, D_diff, D_half = scf_helper.makeMatrices(nbf, 7)

    ## Initial guess (core Hamiltonian) density matrix
    D, homo = scf_helper.diag(H, A, nalpha)

    if DIIS:
        ## Initialize diis object
        diis_obj = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")
    
    ## Nuclear energy
    Enuc = mol.nuclear_repulsion_energy()
    ## Initialize Eold for convegence check
    Eold = 0.0
    
    ## Optional output
    header = "\n    Iter               Energy         HOMO       Delta E         dRMS\n"
    psi4.core.print_out(header)
    if verbose:
        print('\nStarting SCF iterations:')
        print(header)
    t = time.time()
    e_list = []
    e_conv_list = []
    d_conv_list = []

    ## SCF iterative procedure
    for SCF_ITER in range(1, maxiter + 1):
        
        ## Saving the initial density matrix for SCF damping
        D_old = D
        
        ## Build J (Coulomb) and K (Fock exchange)
        J = scf_helper.Jbuild(I, D, J)
        if EXC["name"]=="EXX":
            K = scf_helper.Kbuild(I, D, K)

        ## Build Fock matrix
        F.copy(H)
        F.axpy(1.0, J)
        if EXC["name"]=="EXX":
            F.axpy(-0.5, K)
        if EXC["name"]!="EXX":
            ## Build DFT potentials and calculate energies
            exc, Vxc = scf_helper.Vpot_builder(VXCpot, D, VXC_null, D_half)
            F.axpy(1.0, Vxc)

        if DIIS:
            ## Build DIIS vector
            diis_e, dRMS = scf_helper.diis_vector(F, D, S, A)
            diis_obj.add(F, diis_e)
            dRMS_val = dRMS if dRMS > 0 else 1e-16
            d_conv_list.append(np.log10(dRMS_val))

        ## Energy calculation
        SCF_E = H.vector_dot(D)
        SCF_E += 0.5 * J.vector_dot(D)
        if EXC["name"]=="EXX":
            SCF_E -= 0.25 * K.vector_dot(D)
        else:
            SCF_E += exc
        SCF_E += Enuc

        e_diff = abs(SCF_E - Eold)
        e_diff_val = e_diff if e_diff > 0 else 1e-16
        e_conv_list.append(np.log10(e_diff_val))
        e_list.append(SCF_E)

        ## DIIS convergence test and Fock extrapolation
        if DIIS:
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, homo, (SCF_E - Eold), dRMS)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break
            
            Eold = SCF_E
            F = diis_obj.extrapolate()
        
        ## Diagonalize Fock matrix
        D, homo = scf_helper.diag(F, A, nalpha)

        ## Density convergence check
        if not DIIS:
            
            dRMS = scf_helper.density_RMS(D_diff, D, D_old)
            dRMS_val = dRMS if dRMS > 0 else 1e-16
            d_conv_list.append(np.log10(dRMS_val))
        
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, homo, (SCF_E - Eold), dRMS)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break

            Eold = SCF_E

        ## Dynamic damping
        D, current_damp = scf_helper.dynamic_damping(D, D_old, dRMS, damp, damping_switch_off, current_damp)
        # psi4.core.print_out('\nCurrent damping %.1f\n' % current_damp)
        
        if SCF_ITER == maxiter:
            
            psi4.core.print_out("\nWARNING ! SCF did not converge. The final values are printed\n")
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, D, homo, SCF_ITER, e_list, e_conv_list, d_conv_list
        
    psi4.core.print_out('\nSCF converged in %d iterations and %.3f seconds \n' % (SCF_ITER, time.time() - t))
    if verbose:
        print('\nSCF converged in %d iterations and %.3f seconds' % (SCF_ITER, time.time() - t))

    return SCF_E, D, homo, SCF_ITER, e_list, e_conv_list, d_conv_list