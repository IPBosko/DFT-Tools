def sic_bosko_d(N):
    """
    A closed-shell self-interaction correction
    to Dirac exchange functional derived analytically by Bosko
    """
    return 1 - (2 / N) ** (1 / 3)