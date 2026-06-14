"""
Spin-unrestricted KS SCF solver for DFT calculations using Psi4 objects
"""

import sys
import psi4
import time
import numpy as np
sys.path.append('/Users/ivanbosko/Documents/CODES/GIT/DFT-Tools/SCF')
import scf_helper 

def ks_solver(mol, EXC, damp=0.0, DIIS=True, verbose=False):
    """
    Spin-unrestricted KS SCF solver loop
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
    basic_info_str = f'Number of basis functions: {nbf}\nNumber of alpha electrons: {nalpha}\nNumber of beta electrons: {nbeta}\n'
    psi4.core.print_out(basic_info_str)
    if verbose:
        print(basic_info_str)

    ## Potential Initialization
    build_superfunctional = psi4.driver.dft.build_superfunctional
    VXCpot = scf_helper.Vpot_init(build_superfunctional, wfn, EXC, "UV", restricted=False)
    VXCpot.initialize()

    ## Calculate and store V, T, H (core), ERI (I), and diagonalization matrix (A)
    H, I, S, A = scf_helper.scf_building_blocks(mints)[2:]

    ## Initialize necessary matrices
    D,Da,Db,Fa,Fb,J,Ja,Jb,Ka,Kb,Vxca,Vxcb,VXCa_null,VXCb_null,D_diff = scf_helper.makeMatrices(nbf, 15)

    ## Initial guess (core Hamiltonian) density matrix
    Da = scf_helper.diag(H, A, nalpha, restricted=False)
    Db = scf_helper.diag(H, A, nbeta, restricted=False)

    if DIIS:
        ## Initialize diis object
        diis_obj_a = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")
        diis_obj_b = psi4.p4util.solvers.DIIS(max_vec=6, removal_policy="oldest")
    
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

    ## SCF iterative procedure
    for SCF_ITER in range(1, maxiter + 1):
        
        ## Saving the initial density matrix for SCF damping
        Da_old = Da
        Db_old = Db
        
        ## Build J (Coulomb) and K (Fock exchange)
        Ja = scf_helper.Jbuild(I, Da, Ja)
        Jb = scf_helper.Jbuild(I, Db, Jb)
        J.copy(Ja)
        J.axpy(1.0, Jb)
        if EXC["name"]=="EXX":
            Ka = scf_helper.Kbuild(I, Da, Ka)
            Kb = scf_helper.Kbuild(I, Db, Kb)

        ## Build Fock matrix
        Fa.copy(H)
        Fb.copy(H)
        Fa.axpy(1.0, Ja)
        Fa.axpy(1.0, Jb)
        Fb.axpy(1.0, Ja)
        Fb.axpy(1.0, Jb)
        if EXC["name"]=="EXX":
            Fa.axpy(-1.0, Ka)
            Fb.axpy(-1.0, Kb)
        if EXC["name"]!="EXX":
            ## Build DFT potentials and calculate energies
            exc,Vxc = scf_helper.Vpot_builder(VXCpot, Da, [VXCa_null,VXCb_null], Db, restricted=False)
            Vxca, Vxcb = Vxc[0], Vxc[1]
            Fa.axpy(1.0, Vxca)
            Fb.axpy(1.0, Vxcb)

        if DIIS:
            ## Build DIIS vector
            diis_e_a, dRMS_a = scf_helper.diis_vector(Fa, Da, S, A)
            diis_obj_a.add(Fa, diis_e_a)
            diis_e_b, dRMS_b = scf_helper.diis_vector(Fb, Db, S, A)
            diis_obj_b.add(Fb, diis_e_b)
            dRMS = dRMS_a + dRMS_b

        ## Energy calculation
        SCF_Ea = H.vector_dot(Da)
        SCF_Eb = H.vector_dot(Db)
        SCF_E = SCF_Ea + SCF_Eb
        D.copy(Da)
        D.axpy(1.0,Db)
        SCF_E += 0.5 * J.vector_dot(D)
        if EXC["name"]=="EXX":
            SCF_E -= 0.5 * Ka.vector_dot(Da)
            SCF_E -= 0.5 * Ka.vector_dot(Db)
        else:
            SCF_E += exc
        SCF_E += Enuc

        ## DIIS convergence test and Fock extrapolation
        if DIIS:
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break
            
            Eold = SCF_E
            Fa = diis_obj_a.extrapolate()
            Fb = diis_obj_b.extrapolate()
        
        ## Diagonalize Fock matrix
        Da = scf_helper.diag(Fa, A, nalpha, restricted=False)
        Db = scf_helper.diag(Fb, A, nbeta, restricted=False)

        ## Density convergence check

        if not DIIS:
            
            dRMS_a = scf_helper.density_RMS(D_diff, Da, Da_old)
            dRMS_b = scf_helper.density_RMS(D_diff, Db, Db_old)
            dRMS = dRMS_a + dRMS_b
        
            output_str = 'SCF Iter%3d: % 18.8f   % 1.5E   % 1.5E\n' % (SCF_ITER, SCF_E, (SCF_E - Eold), dRMS)
            psi4.core.print_out(output_str)
            if verbose:
                print(output_str.strip())
            if (abs(SCF_E - Eold) < E_conv and dRMS < D_conv):
                break

            Eold = SCF_E

        ## Dynamic damping
        Da, current_damp = scf_helper.dynamic_damping(Da, Da_old, dRMS_a, damp, damping_switch_off, current_damp)
        Db, current_damp = scf_helper.dynamic_damping(Db, Db_old, dRMS_b, damp, damping_switch_off, current_damp)
        
        if SCF_ITER == maxiter:
            
            psi4.core.print_out("\nWARNING ! SCF did not converge. The final values are printed\n")
            print("\nWARNING ! SCF did not converge. The final values are printed")
            return SCF_E, Da, Db, SCF_ITER
        
    psi4.core.print_out('\nSCF converged in %d iterations and %.3f seconds \n' % (SCF_ITER, time.time() - t))
    if verbose:
        print('\nSCF converged in %d iterations and %.3f seconds' % (SCF_ITER, time.time() - t))

    return SCF_E, Da, Db, SCF_ITER