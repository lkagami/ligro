#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  bond.py
#  
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
#        
#        ═══════════════════════════════════════════════════════════════════════════════
#          PIPELINE DE LIGAÇÕES — VismolObject (EasyHybrid3 / vismol)
#        ═══════════════════════════════════════════════════════════════════════════════
#        
#        ┌─────────────────────────────────────────────────────────────────────────────┐
#        │ ENTRADA: coordenadas (frames[0], shape (N,3)) + elementos (atoms[i].symbol)  │
#        └─────────────────────────────────────────────────────────────────────────────┘
#                                            │
#                                            ▼
#                ┌───────────────────────────────────────────────────┐
#                │ 1. CONECTIVIDADE (quais átomos estão ligados)       │
#                │    Critério por raios covalentes:                   │
#                │      d_ij < (r_i + r_j) * tolerância                │
#                │    → produz self.index_bonds  [i0,j0, i1,j1, ...]   │
#                │      (array achatado de pares)                      │
#                └───────────────────────────────────────────────────┘
#                                            │
#                                            ▼
#                ┌───────────────────────────────────────────────────┐
#                │ 2. CRIA OS OBJETOS Bond                             │
#                │    _bonds_from_pair_of_indexes_list()               │
#                │      • aplica exclude_list (ex.: H-H)               │
#                │      • monta self.bonds = {(i,j): Bond}  (dict,     │
#                │        chave normalizada menor→maior)               │
#                │      • atom.bonds.append(bond) nos dois átomos      │
#                └───────────────────────────────────────────────────┘
#                                            │
#                                            ▼
#                            ╔═══════════════════════════╗
#                            ║  Há ordens externas?      ║
#                            ║  (external_orders)         ║
#                            ╚═══════════════════════════╝
#                               │ SIM              │ NÃO
#                               ▼                  ▼
#                ┌──────────────────────┐   ┌──────────────────────────────────────┐
#                │ 3a. USA AS ORDENS    │   │ 3b. PERCEPÇÃO DE ORDEM               │
#                │     EXTERNAS         │   │     _perceive_bond_orders()          │
#                │  (MOL2/SDF/pDynamo   │   │     [flag multiple_bonds == True]    │
#                │   já trazem ordens)  │   │                                      │
#                │  → respeita o arquivo│   │  módulo bond_order_perception.py     │
#                └──────────────────────┘   │  (estilo Antechamber / Wang & Case): │
#                               │           │   • valência-alvo por elemento (APS) │
#                               │           │   • modelo de oxiânion (nitro,       │
#                               │           │     carboxilato, sulfato, fosfato)   │
#                               │           │   • busca por menor penalidade (tps) │
#                               │           │  → self.bond_order_list [1,2,1,3,...]│
#                               │           │     (1 ordem por ligação, alinhada   │
#                               │           │      a index_bonds)                  │
#                               │           └──────────────────────────────────────┘
#                               │                          │
#                               └──────────┬───────────────┘
#                                          ▼
#                ┌───────────────────────────────────────────────────┐
#                │ 4. ORDEM POR ÁTOMO (para o VBO do shader)         │
#                │    _build_bond_order_per_atom()                   │
#                │      cada átomo recebe a MAIOR ordem entre suas   │
#                │      ligações  → self.bond_order_per_atom         │
#                │    (obs.: a representação recalcula isso na hora  │
#                │     do desenho, alinhado ao frame — ver passo 7)  │
#                └───────────────────────────────────────────────────┘
#                                            │
#                                            ▼
#                ┌───────────────────────────────────────────────────┐
#                │ 5. CLASSIFICAÇÃO METAL / COVALENTE                │
#                │    _detect_metal_bonds()                          │
#                │    [flag metal_dashed_bonds == True]              │
#                │                                                   │
#                │    módulo metal_bonds.py:                         │
#                │      ligação é METÁLICA se ≥1 átomo é metal       │
#                │      (transição + alcalinos + alcalino-terrosos)  │
#                │                                                   │
#                │    → self.covalent_bonds  (sem metal)             │
#                │    → self.metal_bonds     (com metal)             │
#                │    (index_bonds e bonds NÃO mudam — só roteamento)│
#                └───────────────────────────────────────────────────┘
#                                            │
#                                            ▼
#        ═══════════════════════════════════════════════════════════════════════════════
#          DESENHO  (create_representation / draw)
#        ═══════════════════════════════════════════════════════════════════════════════
#                                            │
#                    ┌───────────────────────┼───────────────────────┐
#                    ▼                       ▼                       ▼
#          ┌──────────────────┐   ┌───────────────────┐   ┌──────────────────────┐
#          │ STICKS / LINES   │   │ METAL_DASH        │   │ (demais reps:        │
#          │ índices =        │   │ índices =         │   │  spheres, dots, ...) │
#          │  _covalent_      │   │  self.metal_bonds │   └──────────────────────┘
#          │  indexes()       │   │                   │
#          │                  │   │ DashedLines-      │
#          │ (covalent_bonds  │   │ Representation    │
#          │  ou index_bonds  │   │ (linha pontilhada)│
#          │  se flag off)    │   │ _ensure_metal_    │
#          └──────────────────┘   │  dash()           │
#                    │            └───────────────────┘
#                    ▼
#          ┌─────────────────────────────────────────────────────────┐
#          │ 7. VBO DE ORDEM (sticks)                                  │
#          │    _compute_bond_order_per_vertex(n_atoms)               │
#          │    [flag multiple_bonds: se off → tudo 1 (simples)]      │
#          │      lê index_bonds + bond_order_list                    │
#          │      dimensiona pelo FRAME (n_atoms = frames[0].shape[0])│
#          │      → atributo inteiro vert_bond_order (1/2/3 por átomo)│
#          └─────────────────────────────────────────────────────────┘
#                    │
#                    ▼
#          ┌─────────────────────────────────────────────────────────┐
#          │ 8. GEOMETRY SHADER (sticks.py)                           │
#          │    lê vert_bond_order e emite, por passada (u_pass):     │
#          │      ordem 1 → 1 cilindro                                │
#          │      ordem 2 → 2 cilindros paralelos (offset lateral)    │
#          │      ordem 3 → 3 cilindros em base triangular            │
#          │    (múltiplas usam raio reduzido eff_rad)                │
#          └─────────────────────────────────────────────────────────┘
#        



from vismol.utils.debug import dprint
import numpy as np
import math

class Bond:
    """ Class doc """
    
    def __init__(self, atom_i, atom_j, atom_index_i=None, atom_index_j=None, bond_order=1):
        """ Class initialiser """
        self.atom_i = atom_i
        self.atom_j = atom_j
        # - Remember that the "index" attribute refers to the numbering of atoms 
        # (it is not a zero base, it starts at 1 for the first atom)
        # these indices are zero base numbering (it starts at 0 for the first atom)
        # [EN] 2026-09-28 BUG FIX (user report: KeyError deleting an atom in
        # the Builder, then a white screen of picking-ID colours): this
        # used to be `if atom_index_i:`, so a legitimate index 0 was
        # treated as "not given" and replaced by atom_i.index-1. Every
        # caller passes the real 0-based index explicitly; whenever the
        # 1-based atom.index had drifted from it (Builder deletions
        # renumber atom_id), every bond of atom 0 pointed at the wrong --
        # or a non-existent -- atom. Only a MISSING index (None) falls
        # back to atom.index-1 now.
        if atom_index_i is not None:
            self.atom_index_i = atom_index_i
        else:
            self.atom_index_i = atom_i.index-1
        
        if atom_index_j is not None:
            self.atom_index_j = atom_index_j
        else:
            self.atom_index_j = atom_j.index-1
        
        self.bond_order = bond_order
        self.line_active = True
        self.stick_active = False
    
        self.bond_reference = None


    def get_bond_reference(self):
        """
        Compute reference bond lengths for different bond orders using the
        UFF (Universal Force Field) bond length model.

        The bond length r_ij is estimated as:

            r_ij = r_i + r_j + rBO + rEN

        where:
            r_i, r_j : covalent radii of atoms i and j
            rBO      : bond-order correction term
            rEN      : electronegativity correction term

        Results are stored in:
            self.bond_reference

        Dictionary format:
            {
                1   : single bond reference length,
                1.5 : aromatic bond reference length,
                2   : double bond reference length,
                3   : triple bond reference length
            }
        """

        # Skip recalculation if references already exist
        if self.bond_reference is not None:
            return

        # Empirical constant from UFF
        _lambda = 0.1332

        # Atomic parameters for atom i
        ri = self.atom_i.cov_rad
        Xi = self.atom_i.electronegativity

        # Atomic parameters for atom j
        rj = self.atom_j.cov_rad
        Xj = self.atom_j.electronegativity

        self.bond_reference = {}

        # Bond orders considered in this model
        # 1.5 corresponds to aromatic bonds
        bond_orders = [1, 1.5, 2]#, 3]

        for order in bond_orders:

            # Bond-order correction:
            # Higher bond order -> shorter bond
            rBO = -_lambda * (ri + rj) * math.log(order)

            # Electronegativity correction:
            # Accounts for asymmetric electron density distribution
            rEN = (
                ri * rj * (math.sqrt(Xi) - math.sqrt(Xj)) ** 2
            ) / (Xi * ri + Xj * rj)

            # Final reference bond length
            r_ij = ri + rj + rBO + rEN

            self.bond_reference[order] = r_ij
            

    def get_bond_order(self):
        """
        Determine bond order using nearest-neighbor classification.

        The observed interatomic distance is compared against the ideal
        bond lengths stored in self.bond_reference.

        The bond order whose reference distance is closest to the observed
        distance is selected.

        Example:
            observed distance = 1.39 Å

            reference distances:
                single   -> 1.54
                aromatic -> 1.40
                double   -> 1.34
                triple   -> 1.20

            Result:
                bond order = 1.5 (aromatic)
        """
        
        return 1

        # Ensure reference bond lengths are available
        if self.bond_reference is None:
            self.get_bond_reference()

        # Current measured distance between atoms
        distance = self.distance()

        # Find the bond order whose ideal bond length is closest
        # to the measured distance.
        #
        # min(..., key=...) iterates over all bond orders and returns
        # the one minimizing:
        #
        #     |observed_distance - reference_distance|
        #
        self.bond_order = min(
            self.bond_reference,
           
            key=lambda order: abs(
                distance - self.bond_reference[order]
            )
        )
        dprint(self.atom_i.name,self.atom_j.name ,distance, self.bond_reference, self.bond_order)
        return self.bond_order


    def get_bond_order_old(self):
        """Determina ordem de ligação com base na distância."""

        if self.bond_reference is None:
            self.get_bond_reference()

        dist = self.distance()

        refs = sorted(
            self.bond_reference.items(),
            key=lambda item: item[1]
        )

        self.bond_order = 0   # default: sem ligação

        for order, ref_dist in refs:
            if dist <= ref_dist:
                self.bond_order = order
                break

        return self.bond_order

    def get_bond_reference_old(self):
        """
        Calcula comprimentos de ligação r_ij para diferentes ordens
        usando o modelo do UFF.
        """
        if self.bond_reference is None:
            _lambda = 0.1332

            ri = self.atom_i.cov_rad
            rj = self.atom_j.cov_rad

            Xi = self.atom_i.electronegativity
            Xj = self.atom_j.electronegativity

            self.bond_reference = {}

            bond_orders = [1,  2, 3, 1.5]

            for n in bond_orders:
                rBO = -_lambda * (ri + rj) * math.log(n)

                rEN = (ri * rj * (math.sqrt(Xi) - math.sqrt(Xj))**2) / (Xi * ri + Xj * rj)

                r_ij = ri + rj + rBO + rEN

                self.bond_reference[n] = r_ij

    def distance (self, frame=0):
        """ Function doc """
        vec = self.atom_i.coords(frame) - self.atom_j.coords(frame)
        return np.linalg.norm(vec)
