"""XC functional dictionaries for use with RKS and tutorials.

Exported symbols:
- EXX, LDAx : dict
- FUNCTIONALS : dict mapping alias -> dict

Place additional functionals here and import them from notebooks or code.
"""

EXX = {
    "name": "EXX",
    "x_hf": {"alpha": 1.00},
    "x_functionals": {"LDA_X": {"alpha": 0.00}},
    "c_functionals": {"GGA_K_VW": {"alpha": 0.00}}
}

LDAx = {
    "name": "LDAx",
    "x_hf": {"alpha": 0.00},
    "x_functionals": {"LDA_X": {"alpha": 1.00}},
    "c_functionals": {"LDA_C_VWN_RPA": {"alpha": 0.00}}
}

PBEx = {
    "name": "PBEx",
    "x_hf": {"alpha": 0.00},
    "x_functionals": {"GGA_X_PBE": {"alpha": 1.00}},
    "c_functionals": {"GGA_C_PBE": {"alpha": 0.00}}
}

# Registry for lookup by name
FUNCTIONALS = {
    "EXX": EXX,
    "LDAx": LDAx,
    "PBEx": PBEx
}

__all__ = ["EXX", "LDAx", "PBEx", "FUNCTIONALS"]
