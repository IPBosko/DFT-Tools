import psi4
import numpy as np

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