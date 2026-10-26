"""
metal_bonds.py

Identificacao de ligacoes de coordenacao metalica para representacao especial
(linhas pontilhadas) no EasyHybrid/vismol.

Uma ligacao e considerada "metalica" quando PELO MENOS UM dos seus dois atomos
e um metal. O conjunto de metais inclui metais de transicao + alcalinos +
alcalino-terrosos (configuravel via METAL_ELEMENTS).

Funcoes:
    is_metal(symbol) -> bool
    split_metal_bonds(elements, index_bonds) -> (covalent_flat, metal_flat)
        Recebe o array ACHATADO de pares [i0,j0, i1,j1, ...] e devolve dois
        arrays achatados: as ligacoes covalentes (sem metal) e as metalicas.
"""

import numpy as np


# Metais de transicao (blocos d e f) + alcalinos + alcalino-terrosos.
# Conforme escolhido: transicao + alcalinos/alcalino-terrosos.
METAL_ELEMENTS = {
    # alcalinos
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr',
    # alcalino-terrosos
    'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
    # metais de transicao (bloco d), periodos 4-7
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs',
    # lantanideos
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    # actinideos
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
    'Es', 'Fm', 'Md', 'No', 'Lr',
}


def _norm_symbol(sym):
    """ Normaliza o simbolo: 'FE' / 'fe' -> 'Fe'. Aceita simbolos com numero
        (ex.: 'Fe2') pegando so as letras iniciais. """
    s = ''.join(ch for ch in str(sym) if ch.isalpha())
    if not s:
        return s
    return s[0].upper() + s[1:].lower()


def is_metal(symbol):
    """ True se o simbolo quimico for um metal (segundo METAL_ELEMENTS). """
    return _norm_symbol(symbol) in METAL_ELEMENTS


def split_metal_bonds(elements, index_bonds):
    """ Separa as ligacoes em covalentes (sem metal) e metalicas (>=1 metal).

        elements    : lista de simbolos quimicos (1 por atomo)
        index_bonds : array ACHATADO de pares [i0,j0, i1,j1, ...]

        Retorna (covalent_flat, metal_flat): dois arrays numpy uint32
        achatados no MESMO formato de index_bonds. Qualquer um pode ser vazio.
    """
    ib = np.asarray(index_bonds).ravel()
    n_pairs = len(ib) // 2
    metal_mask = np.zeros(len(elements), dtype=bool)
    for idx, sym in enumerate(elements):
        if is_metal(sym):
            metal_mask[idx] = True

    cov = []
    met = []
    for k in range(n_pairs):
        i = int(ib[2 * k]); j = int(ib[2 * k + 1])
        is_metallic = False
        if i < len(metal_mask) and metal_mask[i]:
            is_metallic = True
        if j < len(metal_mask) and metal_mask[j]:
            is_metallic = True
        if is_metallic:
            met.extend((i, j))
        else:
            cov.extend((i, j))

    return (np.array(cov, dtype=np.uint32),
            np.array(met, dtype=np.uint32))
