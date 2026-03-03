def sic_p(N):
    return (1 - (2 / N) ** (2 / 3))

def sic_p_absp(N):
    return (1 - (2 / N) ** (1 / 3))

def sic_x(N):
    return (1 - (2 / N) ** (1 / 3))

def sic_p_test(N):
    return (1 - (1.412 * (N) ** (- 1 / 3)))