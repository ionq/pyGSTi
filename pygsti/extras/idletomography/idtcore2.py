#***************************************************************************************************
# Copyright 2015, 2019, 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#***************************************************************************************************
from collections import defaultdict
from collections.abc import Iterable
import itertools
from dataclasses import dataclass
import time
from typing import Literal
from abc import ABC, abstractmethod

#import scipy
import scipy.linalg as _spl
from torch import mode
import stim as _stim

import numpy as _np

import pygsti.baseobjs as _baseobjs
from pygsti.baseobjs.nicelyserializable import NicelySerializable as _NicelySerializable
from ...circuits.circuit import Circuit as _Circuit
from ...baseobjs.label import Label as _Lbl
from pygsti.baseobjs.verbosityprinter import VerbosityPrinter
from pygsti.errorgenpropagation.localstimerrorgen import LocalStimErrorgenLabel as _LSE
from pygsti.baseobjs.errorgenlabel import LocalElementaryErrorgenLabel
from pygsti.circuits.cloudcircuitconstruction import create_kcoverage_template, _check_kcoverage_template
from pygsti.tools import errgenproptools as eprop
from pygsti.tools.listtools import remove_duplicates_in_place
from pygsti.protocols.protocol import CircuitListsDesign, CombinedExperimentDesign
from pygsti.data.dataset import DataSet
import pygsti.protocols as _proto
from pygsti.algorithms.germselection import compact_EVD as _compact_EVD

#type PauliBasisMap = dict[str, tuple[str, ...]] # if using Python 3.12+
PauliBasisMap = dict[str, tuple[str, ...]]
LEEL = LocalElementaryErrorgenLabel

# from pygsti.errorgenpropagation.errorpropagator import ErrorGeneratorPropagator


#***************************************************************************************************
# Helper objects
# MOVE to pauliobjs.py and replace what's there if we want to keep this.
#***************************************************************************************************

class NQPauliState(object):
    """
    A N-qubit state that is the tensor product of N
    1-qubit Pauli eigenstates.  These can be represented as
    a string of Xs, Ys and Zz (but not Is) *each* with a +/-
    sign indicating which of the two eigenstates is meant.

    A NQPauliState object can also be used to represent a POVM
    whose effects are the projections onto the 2^N tensor products
    of (the given) Pauli eigenstates.  The +/- sign in this case
    indicates which eigenstate is equated with the "0" (vs "1") outcome.
    """

    def __init__(self, string_rep, signs=None):
        """
        Create a NQPauliState

        Parameters
        ----------
        string_rep : str
            A string with letters in {X,Y,Z} (note: I is not allowed!),
            specifying the Pauli basis for each qubit.

        signs : tuple, optional
            A tuple of 0s and/or 1s.  A zero means the "+" eigenvector is
            either prepared or corresponds to the "0" outcome (if this
            NQPauliState is used to describe a measurment basis).  A one
            means the opposite: the "-" eigenvector is prepared and it
            corresponds to a "0" outcome.  The default is all zeros.
        """
        assert("I" not in string_rep), "'I' cannot be in a NQPauliState"
        self.rep = string_rep
        if signs is None:
            signs = (0,) * len(self.rep)
        self.signs = signs

    def __len__(self):
        return len(self.rep)

    def __str__(self):
        sgn = {1: '+', -1: '-'}
        return "".join(["%s%s" % (sgn[s], let)
                        for s, let in zip(self.signs, self.rep)])

    def __repr__(self):
        return "State[" + str(self) + "]"

    def __eq__(self, other):
        return (self.rep == other.rep) and (self.signs == other.signs)

    def __hash__(self):
        return hash(str(self))

    def to_fiducial_circuit(self, pauli_basis_map: PauliBasisMap):
        """
        Convert this Pauli basis state to a fiducial operation sequence.

        When the returned operation sequence follows a preparation in the `|0...0>`
        Z-basis state or is followed by a Z-basis measurement (with all "+"
        signs), then the Pauli state preparation or measurement described by
        this object will be performed.

        Parameters
        ----------
        pauli_basis_map : dict
            A dictionary w/keys like `"+X"` or `"-Y"` and values that
            are tuples of gate *names* (not labels, which include qubit or
            other state-space designations), e.g. `("Gx","Gx")`.  This
            dictionary describes how to prepare or measure in Pauli bases.

        Returns
        -------
        Circuit
        """
        opstr = []
        sgn = {1: '+', -1: '-'}
        nQubits = len(self.signs)
        for i, (s, let) in enumerate(zip(self.signs, self.rep)):
            key = sgn[s] + let  # e.g. "+X", "-Y", etc
            if key not in pauli_basis_map:
                raise ValueError("'%s' is not in `pauli_basis_map` (keys = %s)"
                                 % (key, str(list(pauli_basis_map.keys()))))
            opstr.extend([_Lbl(opname, i) for opname in pauli_basis_map[key]])
            # pauli_basis_map just has 1Q gate *names* -- need to make into labels
        return _Circuit(opstr, num_lines=nQubits).parallelize()
    
    def restrict_to_qubits(self, qubits):
        new_rep = ''.join(self.rep[i] for i in qubits)
        new_signs = tuple(self.signs[i] for i in qubits)
        return NQPauliState(new_rep, new_signs)


class NQPauliOp(object):
    """
    A N-qubit pauli operator, consisting of
    a 1-qubit Pauli operation on each qubits.
    """

    def __init__(self, string_rep):
        """
        Create a NQPauliOp.

        Parameters
        ----------
        string_rep : str
            A string with letters in {I,X,Y,Z}, specifying the Pauli operator
            for each qubit.
        """
        self.rep = string_rep

    def __len__(self):
        return len(self.rep)

    def __str__(self):
        return self.rep

    def __repr__(self): 
        return "NQPauliOp[%s]" % self.rep

    def __eq__(self, other):
        return (self.rep == other.rep)
    
    def __hash__(self):
        return hash(self.rep)
    
    def to_fiducial_circuit(self, pauli_basis_map: PauliBasisMap):
        """
        Convert this Pauli operator to a fiducial operation sequence.

        When the returned operation sequence follows a preparation in the `|0...0>`
        Z-basis state or is followed by a Z-basis measurement (with all "+"
        signs), then the Pauli state preparation or measurement described by
        this object will be performed.

        Parameters
        ----------
        pauli_basis_map : dict
            A dictionary w/keys like `"+X"` or `"-Y"` and values that
            are tuples of gate *names* (not labels, which include qubit or
            other state-space designations), e.g. `("Gx","Gx")`.  This
            dictionary describes how to prepare or measure in Pauli bases.

        Returns
        -------
        Circuit
        """
        opstr = []
        nQubits = len(self.rep)
        for i, p in enumerate(self.rep):
            if p not in pauli_basis_map:
                raise ValueError("'%s' is not in `pauli_basis_map` (keys = %s)"
                                 % (p, str(list(pauli_basis_map.keys()))))
            opstr.extend([_Lbl(opname, i) for opname in pauli_basis_map[p]])
            # pauli_basis_map just has 1Q gate *names* -- need to make into labels
        return _Circuit(opstr, num_lines=nQubits).parallelize()
    
    def restrict_to_qubits(self, qubits):
        new_rep = ''.join(self.rep[i] for i in qubits)
        return NQPauliOp(new_rep)
    
    def identity_except_on_qubits(self, qubits):
        new_rep = ''.join(self.rep[i] if i in qubits else 'I' for i in range(len(self.rep)))
        return NQPauliOp(new_rep)
    
    def support_indices(self):
        return tuple(i for i, p in enumerate(self.rep) if p != 'I')

#***************************************************************************************************
# Utility functions
#***************************************************************************************************

def add_independent_rows(rows_to_add: Iterable, mx_start: Iterable | None = None):
    if mx_start is None:
        mx = []    
        cur_rank = 0
    else:        
        mx = list(mx_start)
        cur_rank = _np.linalg.matrix_rank(mx)
    
    kept_rows = []
    for i, row in enumerate(rows_to_add):
        test = mx + [row]
        test_r = _np.linalg.matrix_rank(_np.array(test, 'd'))
        if test_r > cur_rank:
            mx = test
            kept_rows.append(i)
            cur_rank = test_r    
    return _np.array(mx), kept_rows, cur_rank


def meas_fidpair_to_meas_paulis(full_weight_meas, max_meas_wt: int | None = None):
    nq = len(full_weight_meas.rep)
    maxwt = max_meas_wt if max_meas_wt is not None else nq
    meas_paulis = []
    assert all(p in 'XYZ' for p in full_weight_meas.rep), f"Measurement fiducial {full_weight_meas} must be full weight (no I's)!"
    for meas_wt in reversed(range(1, maxwt + 1)):
        for support in itertools.combinations(range(nq), meas_wt):
            meas_paulis.append(
                NQPauliOp(''.join([(full_weight_meas.rep[i] if i in support else 'I') for i in range(nq)]))
            )
    return meas_paulis

#***************************************************************************************************
# Jacobian construction routines
#***************************************************************************************************

def _compute_prep_tableau(prep_pauli: NQPauliState, prep_pauli_basis_map: PauliBasisMap):
    prepc = prep_pauli.to_fiducial_circuit(prep_pauli_basis_map)
    prep_tableau = prepc.convert_to_stim_tableau()
    if prep_tableau is None:  # then identity
        nqubits = len(prepc.line_labels)
        prep_tableau = _stim.Tableau(nqubits)
    return prep_tableau

def get_alpha(eeg: LocalElementaryErrorgenLabel, prep_pauli: NQPauliState, meas_pauli: NQPauliOp, prep_pauli_basis_map: PauliBasisMap,
              prep_tableau_cache: dict[NQPauliState, _stim.Tableau] | None = None):
    if prep_tableau_cache is not None:
        if prep_pauli not in prep_tableau_cache:
            prep_tableau_cache[prep_pauli] = _compute_prep_tableau(prep_pauli, prep_pauli_basis_map)
        prep_tableau = prep_tableau_cache[prep_pauli]
    else:
        prep_tableau = _compute_prep_tableau(prep_pauli, prep_pauli_basis_map)
    meas_paulistr = _stim.PauliString(meas_pauli.rep)
    return eprop.alpha_pauli(eeg, prep_tableau, meas_paulistr)

def construct_jacobian(eegs: list[LocalElementaryErrorgenLabel], prep_meas_pairs: list[tuple[NQPauliState, NQPauliOp]],
                       prep_pauli_basis_map: PauliBasisMap,
                       prep_tableau_cache: dict[NQPauliState, _stim.Tableau] | None = None):
    if prep_tableau_cache is None:
        prep_tableau_cache = {}
    num_prepmeas = len(prep_meas_pairs)  
    num_eegs = len(eegs)
    jac = _np.zeros((num_prepmeas, num_eegs), 'd')          
    for j, eeg in enumerate(eegs):
        for i, (prep_pauli, meas_pauli) in enumerate(prep_meas_pairs): 
            a = get_alpha(eeg, prep_pauli, meas_pauli, prep_pauli_basis_map=prep_pauli_basis_map, prep_tableau_cache=prep_tableau_cache)
            if a != 0:
                jac[i, j] = a
    return jac

#***************************************************************************************************
# Circuit construction routines
#***************************************************************************************************

def tile_pauli_fidpairs(base_fidpairs: list[tuple[NQPauliState, NQPauliState]], nqubits: int, weight: int, drop_ident: bool):
    """
    Tiles a set of base fiducial pairs on `weight` qubits to a
    set of fiducial pairs on `nqubits` qubits such that every set
    of `weight` qubits takes on the values in each base pair in
    at least one of the returned pairs.

    Parameters
    ----------
    base_fidpairs : list
        A list of 2-tuples of :class:`NQPauliState` objects (on `maxweight`
        qubits).

    nqubits : int
        The number of qubits.

    maxweight : int
        The maximum weight errors the base qubits are meant to
        detect.  Equal to the number of qubits in the base pairs.

    Returns
    -------
    list
        A list of 2-tuples of :class:`NQPauliState` objects (on `nqubits`
        qubits).
    """
    nqubit_fidpairs = []
    tmpl = create_kcoverage_template(nqubits, weight)
    _check_kcoverage_template(tmpl, nqubits, weight)
    for base_prep, base_meas in base_fidpairs:
        if drop_ident and 'I' in base_meas.rep: continue
        for tmpl_row in tmpl:
            #Replace 0...weight-1 integers in tmpl_row with Pauli basis
            # designations (e.g. +X) to construct NQPauliState objects.
            prep = NQPauliState(''.join(base_prep.rep[i] for i in tmpl_row),
                                signs=tuple(base_prep.signs[i] for i in tmpl_row))
            meas = NQPauliOp(''.join(base_meas.rep[i] for i in tmpl_row))
            nqubit_fidpairs.append((prep, meas))

    remove_duplicates_in_place(nqubit_fidpairs)
    return nqubit_fidpairs

def embed_eeg(eeg, targets, eeg_nq):
    nq_bels = []
    for bel in eeg.basis_element_labels:
        x = ['I'] * eeg_nq
        for i, p in zip(targets, bel):
            x[i] = p
        nq_bels.append(''.join(x))
    return LEEL(eeg.errorgen_type, tuple(nq_bels))

def build_eeg_list_uptoweight(eg_types: tuple[str, ...], uptoweigth: int):
    non_ident_bels = [''.join(p) for p in itertools.product(['I', 'X', 'Y', 'Z'], repeat=uptoweigth)][1:]
    eegs = []
    if 'H' in eg_types:
        eegs.extend([LEEL('H', (bel,)) for bel in non_ident_bels])
    if 'S' in eg_types:
        eegs.extend([LEEL('S', (bel,)) for bel in non_ident_bels])
    if 'C' in eg_types:
        eegs.extend([LEEL('C', (bel1,bel2)) for i, bel1 in enumerate(non_ident_bels) for bel2 in non_ident_bels[i+1:]])
    if 'A' in eg_types:
        eegs.extend([LEEL('A', (bel1,bel2)) for i, bel1 in enumerate(non_ident_bels) for bel2 in non_ident_bels[i+1:]])
    return eegs

def tile_eegs(base_eegs: list[LocalElementaryErrorgenLabel], nqubits: int, weight: int):
    nqubit_eegs = []
    for eeg in base_eegs:
        for inds in itertools.combinations(range(nqubits), weight):
            nqubit_eegs.append(embed_eeg(eeg, inds, nqubits))
    remove_duplicates_in_place(nqubit_eegs)
    return nqubit_eegs

def get_eeg_types_num_qubits_and_max_weight(eegs_to_probe: list[LocalElementaryErrorgenLabel]) -> tuple[set[str], int, int]:
    eeg_nq = set(); eeg_types = set(); eeg_weights = set()
    for eeg in eegs_to_probe:
        eeg_types.add(eeg.errorgen_type)
        bel_lens = set(len(bel) for bel in eeg.basis_element_labels)
        assert len(bel_lens) == 1, f"All basis elements of an elementary error generator must have the same number of qubits.  Found: {bel_lens} for {eeg}!"
        nq = next(iter(bel_lens))
        eeg_nq.add(nq)
        eeg_weights.add(sum(any(bel[i] != 'I' for bel in eeg.basis_element_labels) for i in range(nq)))
    assert len(eeg_nq) == 1, f"All elementary error generators in `eegs_to_probe` must have the same # of qubits.  Found: {eeg_nq}!"
    assert eeg_types.issubset({'H', 'S', 'C', 'A'}), f"Unknown elementary error generator types in {eeg_types}!"
    eeg_num_qubits = next(iter(eeg_nq))
    max_eeg_weight = max(eeg_weights)
    return eeg_types, eeg_num_qubits, max_eeg_weight


class IdleTomographyExperimentDesign(CircuitListsDesign):
    def __init__(self, circuit_lists: list[list[_Circuit]], 
                 fiducial_pairs: list[tuple[NQPauliState, NQPauliOp]],
                 fidpair_to_circuits_map: dict[tuple[NQPauliState, NQPauliOp], list[_Circuit]],
                 extrinsic_rate_circuits: dict[tuple[NQPauliState, NQPauliOp], list[_Circuit]],
                 Ls: list[int],
                 eegs_to_probe: list[LocalElementaryErrorgenLabel],
                 prep_fiducial_map: PauliBasisMap,
                 meas_fiducial_map: PauliBasisMap,
                 jacobian: _np.ndarray | None = None):
        super().__init__(circuit_lists)
        self.eeg_types, self.eeg_num_qubits, self.max_eeg_weight = get_eeg_types_num_qubits_and_max_weight(eegs_to_probe)
        self.fiducial_pairs = fiducial_pairs
        self.fidpair_to_circuits_map = fidpair_to_circuits_map
        self.Ls = Ls
        self.eegs_to_probe = eegs_to_probe
        self.prep_fiducial_map = prep_fiducial_map
        self.meas_fiducial_map = meas_fiducial_map
        self.extrinsic_rate_circuits = extrinsic_rate_circuits
        self.jacobian = jacobian

    @property
    def num_qubits(self):
        return len(self.fiducial_pairs[0][0]) # number of qubits is length of prep/meas Paulis in fidpairs
    
    def get_jacobian(self, method: str = "brute", verbosity: int | VerbosityPrinter = 0):
        if self.jacobian is not None:
            return self.jacobian
        elif method == "brute":
            return self._compute_jacobian_brute(verbosity=verbosity)
        else:
            raise ValueError(f"Unknown method '{method}' for computing Jacobian!")
        
    def _compute_jacobian_brute(self, verbosity: int | VerbosityPrinter = 0):
        printer = VerbosityPrinter.create_printer(verbosity)
        jac_full = _np.zeros((len(self.extrinsic_rate_circuits), len(self.eegs_to_probe)), 'd')
        printer.log(f"Brute-force constructing full jacobian with {len(self.extrinsic_rate_circuits)} rows and {len(self.eegs_to_probe)} columns... ({jac_full.size * 8. / 1e9:.2f} GB)")
        t0 = time.time()
        cache = {}
                            
        for j, nq_eeg in enumerate(self.eegs_to_probe):
            for i, (nq_prep, nq_meas, qubits_to_measure) in enumerate(self.extrinsic_rate_circuits.keys()):   
                nq_meas_on_targets = nq_meas.identity_except_on_qubits(qubits_to_measure)
                jac_full[i, j] = get_alpha(nq_eeg, nq_prep, nq_meas_on_targets,
                                           prep_pauli_basis_map=self.prep_fiducial_map, prep_tableau_cache=cache)

        nnz = (jac_full != 0).sum()      
        printer.log(f"Done brute-force constructing jacobian in {time.time() - t0:.2f} seconds.")
        printer.log(f"Nonzero elements: {nnz} / {jac_full.size} ({100*nnz/jac_full.size:.2f}%)")
        return jac_full


class IdleTomographyCircuitBuilder(_NicelySerializable):
    @classmethod
    def create(cls, eegs_to_probe: list[LocalElementaryErrorgenLabel],
                prep_fiducial_map: PauliBasisMap,
                meas_fiducial_map: PauliBasisMap, 
                idle_germ: str = "Gidle",
                seed: int = 42,
                mode: Literal['bruteRank', 'compactEVD'] = 'compactEVD',
                verbosity: int = 0):
        (fiducial_pairs, pauli_pair_to_fidpair_map, jacobian) = \
            cls._find_fiducials_random(eegs_to_probe, prep_fiducial_map, seed, mode, verbosity=verbosity)
        return cls(eegs_to_probe, prep_fiducial_map, meas_fiducial_map, idle_germ,
                    fiducial_pairs, pauli_pair_to_fidpair_map, jacobian)
    
    @classmethod
    def create_separate_hamiltonian_stochastic(cls, eegs_to_probe: list[LocalElementaryErrorgenLabel],
                prep_fiducial_map: PauliBasisMap,
                meas_fiducial_map: PauliBasisMap, 
                idle_germ: str = "Gidle",
                seed: int = 42,
                mode: Literal['bruteRank', 'compactEVD'] = 'compactEVD',
                verbosity: int = 0):
        hamiltonian_eegs = [eeg for eeg in eegs_to_probe if eeg.errorgen_type == 'H']
        stochastic_eegs = [eeg for eeg in eegs_to_probe if eeg.errorgen_type != 'H']

        (fiducial_pairs, pauli_pair_to_fidpair_map, jacobian) = \
            cls._find_fiducials_random(hamiltonian_eegs, prep_fiducial_map, seed, mode, verbosity=verbosity)
        ham_builder = cls(hamiltonian_eegs, prep_fiducial_map, meas_fiducial_map, idle_germ,
                            fiducial_pairs, pauli_pair_to_fidpair_map, jacobian)
        
        (fiducial_pairs, pauli_pair_to_fidpair_map, jacobian) = \
            cls._find_fiducials_random(stochastic_eegs, prep_fiducial_map, seed, mode, verbosity=verbosity)
        stoch_builder = cls(stochastic_eegs, prep_fiducial_map, meas_fiducial_map, idle_germ,
                            fiducial_pairs, pauli_pair_to_fidpair_map, jacobian)
        return ham_builder, stoch_builder


    def __init__(self, eegs_to_probe: list[LocalElementaryErrorgenLabel],
                prep_fiducial_map: PauliBasisMap,
                meas_fiducial_map: PauliBasisMap, 
                idle_germ: str,
                fiducial_pairs: list[tuple[NQPauliState, NQPauliOp]],
                pauli_pair_to_fidpair_map: dict[tuple[NQPauliState, NQPauliOp], tuple[NQPauliState, NQPauliOp, tuple[int, ...]]],
                jacobian: _np.ndarray
                ):
        
        super().__init__()
        self._dbcoordinates = None
        self.eeg_types, self.eeg_num_qubits, self.max_eeg_weight = get_eeg_types_num_qubits_and_max_weight(eegs_to_probe)
        self.eegs_to_probe = eegs_to_probe
        self.prep_fiducial_map = prep_fiducial_map
        self.meas_fiducial_map = meas_fiducial_map
        self.idle_germ = idle_germ
        self.fiducial_pairs = fiducial_pairs
        self.pauli_pair_to_fidpair_map = pauli_pair_to_fidpair_map
        self.jacobian = jacobian

    def _to_nice_serialization(self):
        state = super()._to_nice_serialization()
        state['eegs_to_probe'] = [str(eeg) for eeg in self.eegs_to_probe]
        state['prep_fiducial_map'] = {k: list(v) for k, v in self.prep_fiducial_map.items()}
        state['meas_fiducial_map'] = {k: list(v) for k, v in self.meas_fiducial_map.items()}
        state['idle_germ'] = self.idle_germ

        state['fiducial_pairs'] = [
            {'prep_rep': prep.rep, 'prep_signs': list(prep.signs), 'meas_rep': meas.rep}
            for prep, meas in self.fiducial_pairs
        ]

        # pauli_pair_to_fidpair_map: keys are (NQPauliState, NQPauliOp), values are (NQPauliState, NQPauliOp, tuple[int,...])
        state['pauli_pair_to_fidpair_map'] = [
            {
                'key_prep_rep': key_prep.rep, 'key_prep_signs': list(key_prep.signs),
                'key_meas_rep': key_meas.rep,
                'val_prep_rep': val_prep.rep, 'val_prep_signs': list(val_prep.signs),
                'val_meas_rep': val_meas.rep,
                'val_support': list(val_support),
            }
            for (key_prep, key_meas), (val_prep, val_meas, val_support) in self.pauli_pair_to_fidpair_map.items()
        ]

        state['jacobian'] = self._encodemx(self.jacobian)
        return state

    @classmethod
    def _from_nice_serialization(cls, state):
        eegs_to_probe = [LocalElementaryErrorgenLabel.cast(s) for s in state['eegs_to_probe']]
        prep_fiducial_map = {k: tuple(v) for k, v in state['prep_fiducial_map'].items()}
        meas_fiducial_map = {k: tuple(v) for k, v in state['meas_fiducial_map'].items()}
        idle_germ = state['idle_germ']

        fiducial_pairs = [
            (NQPauliState(fp['prep_rep'], tuple(fp['prep_signs'])), NQPauliOp(fp['meas_rep']))
            for fp in state['fiducial_pairs']
        ]

        pauli_pair_to_fidpair_map = {}
        for entry in state['pauli_pair_to_fidpair_map']:
            key = (NQPauliState(entry['key_prep_rep'], tuple(entry['key_prep_signs'])),
                   NQPauliOp(entry['key_meas_rep']))
            val = (NQPauliState(entry['val_prep_rep'], tuple(entry['val_prep_signs'])),
                   NQPauliOp(entry['val_meas_rep']),
                   tuple(entry['val_support']))
            pauli_pair_to_fidpair_map[key] = val

        jacobian = cls._decodemx(state['jacobian'])

        return cls(eegs_to_probe, prep_fiducial_map, meas_fiducial_map, idle_germ,
                   fiducial_pairs, pauli_pair_to_fidpair_map, jacobian)

    @property
    def pauli_pairs(self):
        return list(self.pauli_pair_to_fidpair_map.keys())

    @staticmethod
    def _find_fiducials_random(eegs_to_probe: list[LocalElementaryErrorgenLabel], prep_fiducial_map: PauliBasisMap,
                               seed: int = 42, mode: Literal['bruteRank', 'compactEVD'] = 'compactEVD', verbosity: int | VerbosityPrinter = 0):
        printer = VerbosityPrinter.create_printer(verbosity)
        _, eeg_num_qubits, max_eeg_weight = get_eeg_types_num_qubits_and_max_weight(eegs_to_probe)
        num_eegs = len(eegs_to_probe)

        # get EEGs by support
        eegs_by_support = defaultdict(set)
        for eeg in eegs_to_probe:
            s = eeg.support_indices()
            eegs_by_support[s].add(eeg.restrict_to(s))

        # Get base fiducial pairs and jacobian for each support
        base_fidpairs = {}; base_pauli_pairs = {}; base_jacs = {}
        cache = {}
        for support,eegs in eegs_by_support.items():

            # get local fiducials that give full rank on this support
            eegs = tuple(eegs)
            if eegs not in cache:
                fidpairs, pauli_pairs, jac, jac_rank = IdleTomographyCircuitBuilder._find_fiducials_brute(
                    eegs, prep_fiducial_map, tile=False, mode=mode, verbosity=verbosity-2)
                cache[eegs] = (fidpairs, pauli_pairs, jac)
            else:
                fidpairs, pauli_pairs, jac = cache[eegs]
            base_fidpairs[support] = fidpairs
            base_pauli_pairs[support] = pauli_pairs
            base_jacs[support] = jac

            rank_str = f"full-rank={jac_rank} jacobian." if jac_rank == jac.shape[1] \
                else f"RANK-DEFFICIENT JACOBIAN: {jac_rank} < {jac.shape[1]}"
            printer.log(f"Support {support}: {len(eegs)} EEGs: Found {len(fidpairs)} fiducial pairs => {rank_str}")
            #printer.log(f"Jacobian shape: {jac.shape}, rank: {_np.linalg.matrix_rank(jac)}")
            
        # Now try to construct full n-qubit circuits and jacobian by ~"tiling" these support-specific fiducials and EEGs
        # Simple heuristic: 
        #  construct 1Q prep/meas pair per qubit based on base fidpairs that overlap that qubit, then randomly
        #  (or systematically?) combine these into n-qubit fidpairs and attempt to construct full rank jacobian.
        
        fidpairs1Q_by_qubit = []
        for i in range(eeg_num_qubits):
            fidpairs1Q = defaultdict(lambda: 0)
            for support, fidpairs in base_fidpairs.items():
                if i in support:
                    for prep, meas in fidpairs:
                        for k in range(len(support)):
                            p = prep.restrict_to_qubits((k,))
                            m = meas.restrict_to_qubits((k,))
                            fidpairs1Q[(p, m)] += 1
            sorted_fidpairs1Q = dict(sorted(fidpairs1Q.items(), key=lambda x: x[1], reverse=True))
            fidpairs1Q_by_qubit.append(sorted_fidpairs1Q)

        print(f"1Q fidpair counts by qubit:")
        for i, fidpairs1Q in enumerate(fidpairs1Q_by_qubit):
            cnt_dict = {f"{p}/{m}": cnt for (p, m), cnt in fidpairs1Q.items()}
            print(f"Qubit {i}: {cnt_dict}")

        # Greedy addition of fiducial pairs
        fidpairs = []; pauli_pairs = []
        pauli_pair_to_fidpair_map = {}
        attempted_fidpairs = set()
        if mode == 'bruteRank':
            jac = _np.zeros((0, num_eegs))  # Initialize an empty Jacobian
        elif mode == 'compactEVD':
            jtj = _np.zeros((num_eegs, num_eegs))  # Initialize an empty J^T J matrix for compact EVD
            projU = None
            jac_rows_to_concat = []  # Keep track of the Jacobian rows corresponding to the added fiducial pairs, to construct the full Jacobian at the end
        else:
            raise ValueError(f"Unknown mode '{mode}' for finding fiducials!")
        jac_rank = 0
        rng = _np.random.default_rng(seed)
        prep_tableau_cache = {}

        printer.log(f"Attempting to tile support-specific fiducials to achieve Jacobian of rank {num_eegs} on {eeg_num_qubits} qubits...")
        while jac_rank < num_eegs and len(attempted_fidpairs) < 1000:
            # construct candidate fiducial pair by combining 1Q fidpairs for each qubit, prioritizing those that appear most frequently in the base fidpairs
            
            candidate_fidpair = None; ntries = 0
            while not candidate_fidpair:        
                component_fidpairs = [rng.choice(list(fidpairs1Q_by_qubit[i].keys())) for i in range(eeg_num_qubits)]
                candidate_prep = NQPauliState(''.join(p.rep[0] for p, _ in component_fidpairs), tuple(p.signs[0] for p, _ in component_fidpairs))
                candidate_meas = NQPauliOp(''.join(m.rep[0] for _, m in component_fidpairs))
                candidate_fidpair = (candidate_prep, candidate_meas)
                if candidate_fidpair in attempted_fidpairs:
                    candidate_fidpair = None  # already have this one, try again
                    ntries += 1
                    if ntries > 100:  # prevent infinite loop
                        raise RuntimeError("Unable to find a new candidate fiducial pair after 100 tries.")
                else:
                    attempted_fidpairs.add(candidate_fidpair)

            (prep, full_wt_meas) = candidate_fidpair
            prep_meas_pairs = [(prep, meas_pauli) for meas_pauli in meas_fidpair_to_meas_paulis(full_wt_meas, max_meas_wt=max_eeg_weight)
                                    if (prep, meas_pauli) not in pauli_pair_to_fidpair_map] # only add pairs we haven't already added to the jacobian
            jac_rows_for_fidpair = construct_jacobian(eegs_to_probe, prep_meas_pairs, prep_fiducial_map, prep_tableau_cache)
            if mode == 'bruteRank':
                test_jac = _np.concatenate((jac, jac_rows_for_fidpair), axis=0)
                rank = _np.linalg.matrix_rank(test_jac)
            elif mode == 'compactEVD':
                low_rank_rep = IdleTomographyCircuitBuilder._low_rank_rep(jac_rows_for_fidpair)
                if projU is not None:                    
                    rank = jac_rank + IdleTomographyCircuitBuilder._rank_increase_from_update(low_rank_rep, projU, evd_tol=1e-10)
                else:
                    assert jac_rank == 0, "If projU is None, then jac_rank should be 0!"
                    rank = low_rank_rep.shape[1]  # if jac_rank is 0, then the rank increase is just the rank of the new low-rank representation
            else:
                raise ValueError(f"Unknown mode '{mode}' for finding fiducials!")

            if rank > jac_rank:
                printer.log(f"Adding pair (attempt {len(attempted_fidpairs)}): {prep}/{full_wt_meas}  rank: {jac_rank} => {rank} <? {num_eegs}", 2)
                fidpairs.append(candidate_fidpair)
                pauli_pairs.extend(prep_meas_pairs)
                pauli_pair_to_fidpair_map.update({pp: (prep, full_wt_meas, pp[1].support_indices()) for pp in prep_meas_pairs})
                jac_rank = rank

                if mode == 'bruteRank':
                    jac = test_jac
                elif mode == 'compactEVD':
                    jac_rows_to_concat.append(jac_rows_for_fidpair)
                    jtj += low_rank_rep @ low_rank_rep.conj().T  # update J^T J with the low-rank contribution from the new fiducial pair                    
                    rank_check, projU = IdleTomographyCircuitBuilder._compact_rank_and_projector(jtj)
                    assert rank_check == jac_rank, f"Rank check failed!  Expected rank {jac_rank} but got {rank_check} from compact EVD of J^T J!"
                else:
                    raise ValueError(f"Unknown mode '{mode}' for finding fiducials!")
  
        if jac_rank < num_eegs:
            printer.log(f"Could not find full-rank Jacobian (rank {jac_rank} < {num_eegs}) in {len(attempted_fidpairs)} attempted fiducial pairs.")
        else:
            printer.log(f"Full rank ({num_eegs}) achieved with {len(fidpairs)} fiducial pairs after {len(attempted_fidpairs)} attempts.")
        
        if mode == 'compactEVD':
            # construct the full Jacobian from the rows corresponding to the added fiducial pairs
            jac = _np.concatenate(jac_rows_to_concat, axis=0)  

        assert jac.shape == (len(pauli_pairs), num_eegs), "Jacobian shape should be (#prep-meas pairs, #EEGs)"
        assert len(pauli_pairs) == len(set(pauli_pairs)) == len(pauli_pair_to_fidpair_map), "Pauli pairs should be unique"
        assert pauli_pairs == list(pauli_pair_to_fidpair_map.keys()), "Pauli pairs should match keys of pauli_pair_to_fidpair_map"
        return fidpairs, pauli_pair_to_fidpair_map, jac


    @staticmethod
    def _low_rank_rep(jac: _np.ndarray, evd_tol: float = 1e-10) -> _np.ndarray:
        """Compute low rank representation of jacobian `jac`.

        Typically `jac` is a part of a larger jacobian, namely the rows corresponding
        to one or several fiducial pairs.
    

        Parameters
        ----------
        jac
            Jacobian matrix, with shape = (number of prep-meas pairs, num_eegs)
        evd_tol, optional
            Eigenvalue tolerance for compact EVD, by default 1e-10

        Returns
        -------
            Compact representation of "JtJ" jacobian, shape=(num_eegs, rank)
        """
        jjt = jac @ (jac.conj().T)  # (num_pauli_pairs,num_eegs) @ (num_eegs,num_pauli_pairs) => (num_pauli_pairs,num_pauli_pairs)
        e, U = _compact_EVD(jjt, evd_tol)  # jjt hermetian => U unitary, e real
        # U.shape == (num_pauli_pairs, rank), e.shape == (rank,)

        # convert to compact EVD of jac.dag @ jac via left singular vectors:        
        # Multiply U by jac.conj().T and rescale the columns
        # by the corresponding singular value, i.e. the sqrt of the
        # eigenvalue. Use some broadcasting for fast rescaling.
        U_prime = ((jac.conj().T) @ U) / _np.sqrt(e.reshape((1,len(e))))
        ret = U_prime @ _np.diag(_np.sqrt(e)) # ret.shape = (num_eegs, rank)
        assert _np.allclose(ret, jac.conj().T @ U) # check that this is the same!

        # Move this to a unit test:
        assert _np.allclose(ret @ ret.conj().T, jac.conj().T @ jac) # check that this is the compact EVD of jac.dag @ jac
        return ret # a compact representation of "JtJ" jacobian, shape=(num_eegs, rank)


    @staticmethod
    def _compact_rank_and_projector(jtj: _np.ndarray, evd_tol: float = 1e-10) -> tuple[_np.ndarray, _np.ndarray, _np.ndarray]:
        """Compute rank of a hermitian matrix and projecctor onto complement of column space.

        Used for efficient updates of jacobian rank when adding fiducial pairs one at a time,
        without needing to store the full jacobian in memory.

        Parameters
        ----------
        jtj
            square matrix to compute decomposition of with shape = (n, n).  Must be hermitian, 
            e.g. jac.dag @ jac for some jacobian matrix jac.
        evd_tol, optional
            Eigenvalue tolerance for compact EVD, by default 1e-10

        Returns
        -------
            tuple[_np.ndarray, _np.ndarray, _np.ndarray]
                A tuple containing the eigenvalues, eigenvectors, and projector onto the complement of the column space.
        """
        # jtj is assumed to be hermitian, so U should be real and inv(U) = U.conj().T
        e, U = _compact_EVD(jtj, evd_tol) # U.shape = (n, rank), e.shape = (rank,)
        assert _np.isreal(e).all() and _np.isreal(U).all(), "EVD of hermitian matrix should be real!"
    
        #construct the projector, shape = (n, n)
        projU = _np.eye(jtj.shape[0]) - U@U.T # .conj() not needed bc real.
        rank = len(e)
        return rank, projU
    
    @staticmethod
    def _rank_increase_from_update(update_low_rank_rep: _np.ndarray, proj_U: _np.ndarray, evd_tol: float = 1e-10) -> int:
        """Compute additional rank from updating a hermitian matrix `mx = J.dag @ J` with a low-rank update given by `update_low_rank_rep`.

        Parameters
        ----------
        update_low_rank_rep
            The low-rank representation of the update, e.g., from the _low_rank_rep method.  Shape = (n, update_rank)
        proj_U
            The projector onto the complement of the column space of `mx` before the update, shape = (n, n).
        evd_tol, optional
            Eigenvalue tolerance for compact EVD, by default 1e-10

        Returns
        -------
            int
                The increase in rank of `mx` after the low-rank update (>=0).
        """
        # Form projector matrix, whose column space forms an orthonormal basis for the
        #  component of update that is in the complement of U.
        proj_update= proj_U @ update_low_rank_rep
        
        #Next take the RRQR decomposition of this matrix:
        q_update, r_update, _ = _spl.qr(proj_update, mode='economic', pivoting=True)
        
        #Construct P by taking the columns of q_update corresponding to non-zero values of r_A on the diagonal.
        nonzero_indices_update= _np.nonzero(_np.abs(_np.diag(r_update)) > _np.sqrt(evd_tol))
    
        rank_increase = len(nonzero_indices_update[0])
        return rank_increase

    @staticmethod
    def _find_fiducials_brute(eegs_to_probe: list[LocalElementaryErrorgenLabel], prep_fiducial_map: PauliBasisMap,
                              tile: bool = False, mode: Literal['bruteRank', 'compactEVD'] = 'compactEVD',
                              verbosity: int = 0) -> tuple[list[tuple[NQPauliState, NQPauliOp]], list[tuple[NQPauliState, NQPauliOp]], _np.ndarray, int]:
        printer = VerbosityPrinter.create_printer(verbosity)
        eeg_types, eeg_num_qubits, max_eeg_weight = get_eeg_types_num_qubits_and_max_weight(eegs_to_probe)
        num_eegs = len(eegs_to_probe)

        # Construct jacobian - for all (or a group of, with easy mod) errorgens of *exactly* weight=nq
        all_weight_nq_bels = list(itertools.product(['X', 'Y', 'Z'], repeat=eeg_num_qubits))
        printer.log(f"{eeg_num_qubits}-qubit analysis with {eeg_types} eegs.")
        
        if 'C' in eeg_types or 'A' in eeg_types:
            possible_signs = tuple(itertools.product([+1,-1], repeat=eeg_num_qubits))
        else:
            possible_signs = (tuple([+1] * eeg_num_qubits),)

        # start with full-rank measure ops
        possible_preps = [NQPauliState(''.join(p), signs) for p in all_weight_nq_bels for signs in possible_signs]
        possible_meas = [NQPauliOp(''.join(p)) for p in all_weight_nq_bels]
        possible_fidpairs = [(prep_pauli, meas_pauli) for prep_pauli in possible_preps for meas_pauli in possible_meas]
        
        printer.log(f"{num_eegs} EEGs to estimate.")
        printer.log(f"{len(possible_fidpairs)} potential circuits.")
        
        # Greedy addition of fiducial pairs
        fidpairs = []; pauli_pairs = []
        if mode == 'bruteRank':
            jac = _np.zeros((0, num_eegs))  # Initialize an empty Jacobian
            jac_rank = 0
        elif mode == 'compactEVD':
            jtj = _np.zeros((num_eegs, num_eegs), 'd')
            jac_rank = 0

            # If initial fidpairs, construct initial jacobian like this?
            # for j, idx in enumerate(nonzero_weight_indices):
            #     if j==0:
            #         temp_DDD = derivDaggerDeriv[idx] @ derivDaggerDeriv[idx].T
            #     else:
            #         temp_DDD += derivDaggerDeriv[idx] @ derivDaggerDeriv[idx].T
            # currentDDDList.append(temp_DDD)

        else:
            raise ValueError(f"Unknown mode '{mode}' for finding fiducials!")

        # Add pair that gives largest rank increase
        candidate_fidpairs = possible_fidpairs.copy()
        
        jac_rows_for_fidpair = {}
        low_rank_rep_for_fidpair = {}
        prep_tableau_cache = {}
        for i, (prep, full_wt_meas) in enumerate(candidate_fidpairs):
            prep_meas_pairs = [(prep, meas_pauli) for meas_pauli in meas_fidpair_to_meas_paulis(full_wt_meas)]
            printer.log(f"Computing jacobian for candidate {i}/{len(candidate_fidpairs)}: {prep}/{full_wt_meas} (shape = ({num_eegs},{len(prep_meas_pairs)}))", 3)
            jac_rows = construct_jacobian(eegs_to_probe, prep_meas_pairs, prep_fiducial_map, prep_tableau_cache)
            jac_rows_for_fidpair[(prep, full_wt_meas)] = jac_rows                        
            if mode == 'compactEVD':                
                low_rank_rep_for_fidpair[(prep, full_wt_meas)] = IdleTomographyCircuitBuilder._low_rank_rep(jac_rows) # compact low rank update

        while jac_rank < num_eegs and len(fidpairs) < len(possible_fidpairs):

            # Outer loop
            if mode == 'compactEVD':
                rank_check, projU = IdleTomographyCircuitBuilder._compact_rank_and_projector(jtj, evd_tol=1e-10)
                assert rank_check == jac_rank, f"Rank check failed!  Expected rank {jac_rank} but got {rank_check} from compact EVD of J^T J!"

            # Inner loop: compute rank increase for each candidate fiducial pair, and find best one to add
            test_ranks = []
            if mode == 'bruteRank':            
                for i, (prep, full_wt_meas) in enumerate(candidate_fidpairs):                
                    jac_addon = jac_rows_for_fidpair[(prep, full_wt_meas)]
                    test_jac = _np.concatenate((jac, jac_addon), axis=0)
                    rank = _np.linalg.matrix_rank(test_jac)
                    test_ranks.append(rank)

            elif mode == 'compactEVD':
                for i, (prep, full_wt_meas) in enumerate(candidate_fidpairs):                
                    jac_addon_low_rank = low_rank_rep_for_fidpair[(prep, full_wt_meas)]
                    rank = jac_rank + IdleTomographyCircuitBuilder._rank_increase_from_update(jac_addon_low_rank, projU, evd_tol=1e-10)
                    test_ranks.append(rank)

            else:
                raise ValueError(f"Unknown mode '{mode}' for finding fiducials!")
            
            # Candidates that don't increase the rank can be removed from consideration going forward,
            #  since there's no chance of them increasing the rank later (they add the same rows to the
            #  jacobian as they would now, and those rows are linearly dependent on existing rows).
            candidates_to_delete = [i for i, r in enumerate(test_ranks) if r <= jac_rank]

            i = _np.argmax(test_ranks)
            rank = test_ranks[i]

            if rank > jac_rank:
                prep, full_wt_meas = candidate_fidpairs[i]
                printer.log(f"Adding pair (of {len(candidate_fidpairs)} candidates): {prep}/{full_wt_meas}  rank: {jac_rank} => {rank}", 2)
                pauli_pairs_addon = [(prep, meas_pauli) for meas_pauli in meas_fidpair_to_meas_paulis(full_wt_meas)]
                fidpairs.append((prep, full_wt_meas))
                pauli_pairs.extend(pauli_pairs_addon)
                jac_rank = rank
                candidates_to_delete.append(i)  # also remove the candidate we just added

                if mode == 'bruteRank':                    
                    jac_addon = jac_rows_for_fidpair[(prep, full_wt_meas)]
                    jac = _np.concatenate((jac, jac_addon), axis=0)
                    #assert _np.linalg.matrix_rank(jac) == jac_rank # Sanity check (removed for speed)

                elif mode == 'compactEVD':
                    # Update current J.dag @ J matrix (`jtj``):
                    jac_addon_low_rank = low_rank_rep_for_fidpair[(prep, full_wt_meas)]
                    jtj += jac_addon_low_rank @ jac_addon_low_rank.conj().T
                
                else:
                    raise ValueError(f"Unknown mode '{mode}' for finding fiducials!")

                if jac_rank == num_eegs:
                    printer.log(f"Full rank ({num_eegs}) achieved with {len(fidpairs)} fiducial pairs.")
                    break

            # Remove candidates that don't increase rank
            for i in sorted(candidates_to_delete, reverse=True):
                del candidate_fidpairs[i]
        else:
            printer.log(f"Exhausted all {len(possible_fidpairs)} fiducial pairs and Jacobian rank {jac_rank} < {num_eegs}.")
        

        if mode == 'compactEVD':
            # rebuild `jac` since we haven't been updating it and need it to for next step
            jac = _np.concatenate([jac_rows_for_fidpair[fidpair] for fidpair in fidpairs], axis=0)

        orig_shape = jac.shape
        jac, lin_indep_row_inds, rank = add_independent_rows(jac)
        pauli_pairs = [pauli_pairs[i] for i in lin_indep_row_inds]        
        assert rank == jac_rank, f"Rank after restricting to independent rows ({rank}) should match rank beforehand ({jac_rank})."
        if jac.shape[0] < orig_shape[0]:
            printer.log(f"Restricted Jacobian to {jac.shape[0]} independent rows (from {orig_shape[0]}).")
        else:
            printer.log(f"All rows of Jacobian are independent (shape {jac.shape}).")

        if tile:
            # Tile fiducial pairs to get more fiducial pairs to use in the final set of circuits
            fidpairs = tile_pauli_fidpairs(fidpairs, eeg_num_qubits, max_eeg_weight, drop_ident=True)
            pauli_pairs = tile_pauli_fidpairs(pauli_pairs, eeg_num_qubits, max_eeg_weight, drop_ident=False)
            printer.log(f"Tiled fiducial pairs to {len(fidpairs)} pairs, {len(pauli_pairs)} Pauli pairs.")
            jac = construct_jacobian(eegs_to_probe, pauli_pairs, prep_fiducial_map, prep_tableau_cache)
            rank = _np.linalg.matrix_rank(jac)
            assert rank == jac_rank, f"Rank after tiling ({rank}) should match rank beforehand ({jac_rank})."

        return fidpairs, pauli_pairs, jac, jac_rank

    def _create_idt_circuit(self, prep_pauli_basis: NQPauliState, meas_pauli_op: NQPauliOp, L: int):
        assert len(prep_pauli_basis) == len(meas_pauli_op), "prep and meas Paulis must have same number of qubits"
        
        # idle_germ
        nq = len(prep_pauli_basis)
        if isinstance(self.idle_germ, str):
            idle_germ_circuit = _Circuit(self.idle_germ, num_lines=nq)
        else:
            raise ValueError("`idle_germ` must be a string.")
        
        return (prep_pauli_basis.to_fiducial_circuit(self.prep_fiducial_map) 
                + idle_germ_circuit.repeat(L)
                + meas_pauli_op.to_fiducial_circuit(self.meas_fiducial_map))


    def _create_circuits_for_fidpairs(self, fidpairs: list[tuple[NQPauliState, NQPauliOp]], Ls: list[int]):
        fidpair_to_circuits_map = defaultdict(list)
        circuits_by_L = []
        for L in Ls:
            unique_circuits = []
            for prep, meas in fidpairs:
                c = self._create_idt_circuit(prep, meas, L)
                fidpair_to_circuits_map[(prep, meas)].append(c)

                indx = next((i for i, x in enumerate(unique_circuits) if x == c), None)
                if indx is None:                
                    unique_circuits.append(c)
            circuits_by_L.append(unique_circuits)

        assert all(len(fidpair_to_circuits_map[fidpair]) == len(Ls) for fidpair in fidpairs), "Each fidpair should have a circuit for each L!"
        return circuits_by_L, fidpair_to_circuits_map
    

    def create_idt_experiment_design(self, Ls: list[int], verbosity: int | VerbosityPrinter = 0) -> IdleTomographyExperimentDesign:
        printer = VerbosityPrinter.create_printer(verbosity)

        # Create circuits for each fidpair and L
        circuits_by_L, fidpair_to_circuits_map = self._create_circuits_for_fidpairs(self.fiducial_pairs, Ls)
        printer.log(f"Created {sum([len(lst) for lst in circuits_by_L])} unique circuits for {len(self.fiducial_pairs)} fidpairs and {len(Ls)} L values.")

        extrinsic_rate_circuits = {}
        for pauli_pair in self.pauli_pairs:
            nq_prep, nq_meas, measure_on = self.pauli_pair_to_fidpair_map[pauli_pair]
            extrinsic_rate_circuits[(nq_prep, nq_meas, measure_on)] = fidpair_to_circuits_map[(nq_prep, nq_meas)]

        return IdleTomographyExperimentDesign(circuits_by_L, self.fiducial_pairs, fidpair_to_circuits_map,
                                            extrinsic_rate_circuits, Ls, self.eegs_to_probe, 
                                            self.prep_fiducial_map, self.meas_fiducial_map, self.jacobian)

        
def combine_builders(hamiltonian_builder: IdleTomographyCircuitBuilder, stochastic_builder: IdleTomographyCircuitBuilder) -> IdleTomographyCircuitBuilder:
    assert hamiltonian_builder.prep_fiducial_map == stochastic_builder.prep_fiducial_map, "Prep fiducial maps must match to combine builders!"
    assert hamiltonian_builder.meas_fiducial_map == stochastic_builder.meas_fiducial_map, "Meas fiducial maps must match to combine builders!"
    assert hamiltonian_builder.idle_germ == stochastic_builder.idle_germ, "Idle germs must match to combine builders!"

    # check that there are no overlaps in the fiducial pairs used for the hamiltonian vs stochastic errorgens
    if set(hamiltonian_builder.fiducial_pairs).intersection(set(stochastic_builder.fiducial_pairs)):
        raise ValueError("The same fiducial pair cannot be used for both hamiltonian and stochastic errorgens!")
    if set(hamiltonian_builder.pauli_pairs).intersection(set(stochastic_builder.pauli_pairs)):
        raise ValueError("The same pauli pair cannot be used for both hamiltonian and stochastic errorgens!")
        
    eegs_to_probe = hamiltonian_builder.eegs_to_probe + stochastic_builder.eegs_to_probe
    fiducial_pairs = hamiltonian_builder.fiducial_pairs + stochastic_builder.fiducial_pairs
    pauli_pair_to_fidpair_map ={**hamiltonian_builder.pauli_pair_to_fidpair_map, **stochastic_builder.pauli_pair_to_fidpair_map}
    
    npairs = len(hamiltonian_builder.pauli_pairs) + len(stochastic_builder.pauli_pairs)
    jacobian = _np.zeros((npairs,  len(eegs_to_probe)), 'd')
    jacobian[:len(hamiltonian_builder.pauli_pairs), :len(hamiltonian_builder.eegs_to_probe)] = hamiltonian_builder.jacobian
    jacobian[len(hamiltonian_builder.pauli_pairs):, len(hamiltonian_builder.eegs_to_probe):] = stochastic_builder.jacobian

    return IdleTomographyCircuitBuilder(eegs_to_probe, hamiltonian_builder.prep_fiducial_map,
                                        hamiltonian_builder.meas_fiducial_map, hamiltonian_builder.idle_germ,
                                        fiducial_pairs, pauli_pair_to_fidpair_map, jacobian)


#***************************************************************************************************
# Core Idle tomography routines
#***************************************************************************************************

def parity(outcome: str, restrict_to: list[int]):
    return [outcome[i] for i in restrict_to].count('1') % 2

def noise_free_parity(prep: NQPauliState, meas: NQPauliOp):
    if prep.rep == meas.rep:
        if prep.signs.count('-') % 2:
            return -1.0
        else:
            return +1.0
    else:
        return 0.0

class IdleTomography(_proto.Protocol):
    """
    The idle tomography protocol.

    Parameters
    ----------

    verbosity : int, optional
        The 'verbosity' option is an integer specifying the level of
        detail printed to stdout during the calculation.

    name : str, optional
        The name of this protocol, also used to (by default) name the
        results produced by this protocol.  If None, the class name will
        be used.
    """

    def __init__(self, prep_fid_map: PauliBasisMap, meas_fid_map: PauliBasisMap, compute_jacobian_method: str = "brute", name: str = None, verbosity: int = 0):
        super().__init__(name)
        self.prep_fid_map = prep_fid_map
        self.meas_fid_map = meas_fid_map
        self.compute_jacobian_method = compute_jacobian_method
        self.verbosity = verbosity

    
    def run(self, data, memlimit=None, comm=None):
        """
        Run this protocol on `data`.

        Parameters
        ----------
        data : ProtocolData
            The input data.

        memlimit : int, optional
            A rough per-processor memory limit in bytes.

        comm : mpi4py.MPI.Comm, optional
            When not ``None``, an MPI communicator used to run this protocol
            in parallel.

        Returns
        -------
        ModelEstimateResults
        """
        top_edesign = data.edesign

        if isinstance(top_edesign, IdleTomographyExperimentDesign):
            idt_edesigns = [top_edesign]
        elif isinstance(top_edesign, CombinedExperimentDesign):
            idt_edesigns = [e for _, e in top_edesign.items() if isinstance(e, IdleTomographyExperimentDesign)]
            if not idt_edesigns:
                raise ValueError("CombinedExperimentDesign must contain at least one IdleTomographyExperimentDesign!")
        else:
            raise ValueError("Experiment design must be an IdleTomographyExperimentDesign or a CombinedExperimentDesign containing at least one IdleTomographyExperimentDesign!  Found: %s" % type(edesign))
        
        ds = data.dataset

        printer = VerbosityPrinter.create_printer(self.verbosity, comm)
        # if self.record_output and not printer.is_recording():
        #     printer.start_recording()

        extrinsic_rates = {}; intrinsic_rates = {}; jacobians = []
        for i, edesign in enumerate(idt_edesigns):
            # Note: we could also iterate over sub-datasets corresponding to each edesign
            #  but this shouldn't be necessary.

            printer.log(f"Computing extrinsic rates for edesign {i+1}/{len(idt_edesigns)}...")
            local_extrinsic_rates = self.compute_extrinsic_rates(edesign, ds, printer - 1)

            # REMOVE (DEBUG)
            # if set(local_extrinsic_rates.keys()).intersection(set(extrinsic_rates.keys())):
            #     print(list(local_extrinsic_rates.keys()))
            #     print("----")
            #     print(list(extrinsic_rates.keys()))
            #     import pdb; pdb.set_trace()
            #     print("HERE")


            assert not set(local_extrinsic_rates.keys()).intersection(set(extrinsic_rates.keys())), \
                "Extrinsic rates should be disjoint across different edesigns (since they should be probing different sets of errorgens)!"
            assert not set(edesign.eegs_to_probe).intersection(set(intrinsic_rates.keys())), \
                "Intrinsic rates should be disjoint across different edesigns (since they should be probing different sets of errorgens)!"
        
            printer.log(f"Computed {len(local_extrinsic_rates)} extrinsic rates (for {len(edesign.eegs_to_probe)} intrinsic rates).")
            if len(local_extrinsic_rates) < len(edesign.eegs_to_probe):
                printer.log(f"WARNING: This should be more than the {len(edesign.eegs_to_probe)} intrinsic rates!!!") 

            # Compute and use jacobians to convert from extrinsic rates to intrinsic rates
            #  must consider nonzero alpha (jac element) whenever errorgen has overlap with targets
            jac = edesign.get_jacobian(printer)
            
            # Compute intrinsic rates using jacobian
            t1 = time.time()
            extrinsic_vec = _np.array(list(local_extrinsic_rates.values()))
            intrinsic_vec, residuals, solve_rank, solve_sing_vals = _np.linalg.lstsq(jac, extrinsic_vec, rcond=None)
            local_intrinsic_rates = dict(zip(edesign.eegs_to_probe, intrinsic_vec))
            printer.log(f"Solved for intrinsic rates in {time.time() - t1:.2f} seconds.  Solve_rank = {solve_rank}, smallest singular value = {min(solve_sing_vals)}")

            extrinsic_rates.update(local_extrinsic_rates)
            intrinsic_rates.update(local_intrinsic_rates)
            jacobians.append(jac)
        
        # Block diagonal jacobian since we assert different errorgens are probed by different circuits
        jac_full = _np.zeros((len(extrinsic_rates), len(intrinsic_rates)))
        row = 0
        for jac in jacobians:
            nrows, ncols = jac.shape
            jac_full[row:row+nrows, :ncols] = jac
            row += nrows

        return IdleTomographyResults(top_edesign, ds, extrinsic_rates, intrinsic_rates, jac_full)

    def compute_extrinsic_rates(self, edesign: IdleTomographyExperimentDesign, ds: DataSet, verbosity: int | VerbosityPrinter = 0):
        printer = VerbosityPrinter.create_printer(verbosity)
        extrinsic_rates = {}

        for (nq_prep, nq_meas, measure_on), circuits_vs_L in edesign.extrinsic_rate_circuits.items():
            xs = []; ys = []  # for fitting a line to
            for L, circuit in zip(edesign.Ls, circuits_vs_L):
                fs = ds[circuit].fractions
                avg_parity = sum([(-1)**parity(outcome[0], measure_on) * f for outcome, f in fs.items()])
                xs.append(L)
                ys.append(avg_parity) #- nf_parity)
            a, b = _np.polyfit(xs, ys, deg=1)

            printer.log(f"    nqubit fidpair: {(nq_prep, nq_meas, measure_on)}: rate = {a}", 2)
            printer.log(f"     using data => {xs=}, {ys=} ==> {a}*x + {b}", 2)
            extrinsic_rates[(nq_prep, nq_meas, measure_on)] = a

        return extrinsic_rates


@dataclass
class IdleTomographyResults:
    """
    A container for idle tomography results: intrinsic and observable errors,
    along with supporting information.
    """
    edesign: IdleTomographyExperimentDesign
    dataset: DataSet
    #fit_order: int
    extrinsic_rates: dict
    intrinsic_rates: dict
    jacobian: _np.ndarray
