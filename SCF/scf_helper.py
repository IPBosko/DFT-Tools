
import numpy as np
import psi4

def diag(diag, A, nalpha):
    """
    Diagonalizes the Fock matrix and builds the density matrix from N/2 lowest eigenvectors.
    """

    Fp = psi4.core.triplet(A, diag, A, True, False, True)
    nbf = A.shape[0]
    Cp = psi4.core.Matrix(nbf, nbf)
    eigvals = psi4.core.Vector(nbf)
    Fp.diagonalize(Cp, eigvals, psi4.core.DiagonalizeOrder.Ascending)

    C = psi4.core.doublet(A, Cp, False, False)
    Cocc = psi4.core.Matrix(nbf, nalpha)
    Cocc.np[:] = C.np[:, :nalpha]

    D = psi4.core.doublet(Cocc, Cocc, False, True) 
    D.scale(2.0) 
    
    return D, eigvals.np[nalpha-1]

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

def scf_building_blocks(mints):

    # Kinetic energy
    T = mints.ao_kinetic()
    # External potential
    V = mints.ao_potential()
    # H = T + V
    H = T.clone()
    H.add(V)
    # Electron-repulsion integrals (ERI) tensor
    I = np.asarray(mints.ao_eri())
    # Overlap matrix
    S = mints.ao_overlap()
    # Orthogonaliztion matrix A = S^{1/2}
    A = S.clone()
    A.power(-0.5, 1.e-14)