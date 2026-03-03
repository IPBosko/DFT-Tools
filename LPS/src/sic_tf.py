def sic_bosko(N):
    """
    A closed-shell self-interaction correction
    to Thomas-Fermi functional derived analytically by Bosko
    """
    return 1 - (2 / N) ** (2 / 3)

def sic_absp1(N):
    """
    A closed-shell self-interaction correction to Thomas-Fermi
    functional derived in Acharya et al, 1980. Coefficient 
    C_ABSP1 = 1.412 is obtained by fitting to 55 HF neutral atoms
    """
    return 1 - 1.412 * N ** ( - 1 / 3 )

def sic_absp2(N):
    """
    A closed-shell self-interaction correction to Thomas-Fermi
    functional derived in Acharya et al, 1980. Coefficient 
    C_ABSP2 = 1.332 is obtained by fitting to 1200 HF atoms and ions
    """
    return 1 - 1.332 * N ** ( - 1 / 3 )

def sic_gr(N):
    """
    A closed-shell self-interaction correction to Thomas-Fermi
    functional derived by Gazquez and Robles in 1982 (Eq.54) 
    """
    return (1 - 2 / N) * (1 - 1.015 * N ** ( - 1 / 3 ) + 0.150 * N ** ( - 2 / 3 ))