import psi4

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