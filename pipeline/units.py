HARTREE_TO_KJ = 2625.499639      # kJ/mol per Eh
KJ_TO_KCAL = 1 / 4.184
ANGSTROM_TO_BOHR = 1 / 0.529177210903

# Standard atomic masses (amu) for the elements we expect in textbook mechanisms.
MASSES = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Na": 22.990, "P": 30.974, "S": 32.06, "Cl": 35.45, "K": 39.098,
    "Br": 79.904, "I": 126.904, "Li": 6.94, "B": 10.81, "Mg": 24.305,
}


def kj(eh: float) -> float:
    return eh * HARTREE_TO_KJ


def fmt_kj(x_kj: float) -> str:
    """'85.3 kJ/mol (20.4 kcal/mol)'"""
    return f"{x_kj:.1f} kJ/mol ({x_kj * KJ_TO_KCAL:.1f} kcal/mol)"
