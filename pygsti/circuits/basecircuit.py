"""
Refactored Circuit classes: BaseCircuit, EditableCircuit, and StaticCircuit
"""
#***************************************************************************************************
# Copyright 2015, 2019, 2025 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#***************************************************************************************************

from __future__ import annotations
from typing import Dict, Tuple, Union, Optional, List, TYPE_CHECKING
if TYPE_CHECKING:
    try:
        import qiskit
        import stim
    except:
        pass

import collections as _collections
import itertools as _itertools
import warnings as _warnings
from abc import ABC, abstractmethod

import numpy as _np
from pygsti.baseobjs.label import Label as _Label, CircuitLabel as _CircuitLabel
from pygsti.baseobjs import outcomelabeldict as _ld, _compatibility as _compat
from pygsti.tools import internalgates as _itgs
from pygsti.tools import slicetools as _slct
from pygsti.tools.legacytools import deprecate as _deprecate_fn


# Import helper functions from the original circuit module
from pygsti.circuits.circuit import (
    _label_to_nested_lists_of_simple_labels,
    _sslbls_of_nested_lists_of_simple_labels,
    _accumulate_explicit_sslbls,
    _op_seq_str_suffix,
    _op_seq_to_str,
    to_label
)


class BaseCircuit(ABC):
    """
    Abstract base class for quantum circuits.

    This class contains all common logic shared between EditableCircuit and StaticCircuit.
    It provides read-only operations and properties that work for both mutable and immutable
    circuit representations.

    Attributes
    ----------
    default_expand_subcircuits : bool
        By default, expand sub-circuit labels.

    line_labels : tuple
        The line labels (often qubit labels) of this circuit.

    layertup : tuple
        This Circuit's layers as a standard Python tuple of layer Labels.

    tup : tuple
        This Circuit as a standard Python tuple of layer Labels and line labels.

    str : str
        The Python string representation of this Circuit.
    """
    default_expand_subcircuits = True

    @abstractmethod
    def _bare_init(self, labels, line_labels, editable, name='', stringrep=None, occurrence=None,
                   compilable_layer_indices_tup=()):
        """Internal bare initialization - subclasses must implement this."""
        pass

    def __init__(self):
        """Base initialization - subclasses should call this."""
        self._labels = None  # Will be set by subclasses
        self._line_labels = None  # Will be set by subclasses
        self._occurrence_id = None
        self._compilable_layer_indices_tup = ()
        self._name = None
        self.auxinfo = {}

    @classmethod
    def cast(cls, obj):
        """
        Convert `obj` into a circuit of this class type.

        Parameters
        ----------
        obj : object
            Object to convert

        Returns
        -------
        BaseCircuit subclass
        """
        if isinstance(obj, cls): 
            return obj
        if isinstance(obj, BaseCircuit):
            # Convert between subclasses
            if isinstance(cls, StaticCircuit):
                return obj.to_static()
            else:
                return obj.to_editable()
        if isinstance(obj, (tuple, list)): 
            return cls.from_tuple(obj)
        if isinstance(obj, str): 
            return cls(obj)
        raise ValueError("Cannot create a %s object from '%s'" % (cls.__name__, str(type(obj))))

    @classmethod
    @abstractmethod
    def from_tuple(cls, tup):
        """
        Creates a circuit from a tuple.

        Parameters
        ----------
        tup : tuple
            The tuple to convert.

        Returns
        -------
        BaseCircuit subclass
        """
        pass

    @property
    def line_labels(self):
        """The line labels (often qubit labels) of this circuit."""
        return self._line_labels

    @property
    def name(self):
        """
        The name of this circuit.

        Note: the name is *not* a part of the hashed value.
        The name is used to name the CircuitLabel returned from to_label().
        """
        return self._name

    @property
    def occurrence(self):
        """The occurrence id of this circuit."""
        return self._occurrence_id

    @property
    def compilable_layer_indices(self):
        """Tuple of the layer indices corresponding to "compilable" layers."""
        return self._compilable_layer_indices_tup

    @property
    def compilable_by_layer(self):
        """Boolean array indicating whether each layer is "compilable" or not."""
        ret = _np.zeros(self.depth, dtype=bool)
        ret[list(self._compilable_layer_indices_tup)] = True
        return ret

    @property
    @abstractmethod
    def layertup(self):
        """
        This Circuit's layers as a standard Python tuple of layer Labels.

        Returns
        -------
        tuple
        """
        pass

    @property
    def tup(self):
        """
        This Circuit as a standard Python tuple of layer Labels and line labels.

        Returns
        -------
        tuple
        """
        comp_lbl_flag = ('__CMPLBL__',) if self._compilable_layer_indices_tup else ()
        layertup = self.layertup

        if self._occurrence_id is None:
            if self._line_labels in (('*',), ()):
                return layertup + comp_lbl_flag
            else:
                return layertup + ('@',) + self._line_labels + comp_lbl_flag
        else:
            if self._line_labels in (('*',), ()):
                return layertup + ('@', '(', ')') + ('@', self._occurrence_id) + comp_lbl_flag
            else:
                return layertup + ('@',) + self._line_labels + ('@', self._occurrence_id) + comp_lbl_flag

    @property
    @abstractmethod
    def str(self):
        """
        The Python string representation of this Circuit.

        Returns
        -------
        str
        """
        pass

    @property
    def layerstr(self):
        """Just the string representation of the circuit layers (no '@<line_labels>' suffix)"""
        return self._labels_lines_str()[0]

    @property
    def linesstr(self):
        """Just the string representation of the circuit's line labels (the '@<line_labels>' suffix)"""
        return self._labels_lines_str()[1]

    def _labels_lines_str(self):
        """Split the string representation up into layer-labels & line-labels parts"""
        if '@' in self.str:
            return self.str.split('@')
        else:
            return self.str, ''

    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        return self._labels.__iter__()

    def __contains__(self, x):
        """Note: this is not covered by __iter__ for case of contained CircuitLabels"""
        return any([(x == layer or x in layer) for layer in self._labels])

    def __eq__(self, x):
        if isinstance(x, BaseCircuit):
            if len(self) != len(x):
                return False
            else:
                return self.tup == x.tup
        elif x is None:
            return False
        else:
            tup_x = tuple(x)
            if len(self.layertup) != len(tup_x):
                return False
            else:
                return self.layertup == tup_x

    def __lt__(self, x):
        if isinstance(x, BaseCircuit):
            return self.tup < x.tup
        else:
            return self.layertup < tuple(x)

    def __gt__(self, x):
        if isinstance(x, BaseCircuit):
            return self.tup > x.tup
        else:
            return self.layertup > tuple(x)

    @property
    def num_lines(self):
        """
        The number of lines in this circuit.

        Returns
        -------
        int
        """
        return len(self._line_labels)

    @property
    def num_layers(self):
        """
        The number of layers in the circuit.

        Returns
        -------
        int
        """
        return len(self._labels)

    @property
    def depth(self):
        """
        The circuit depth.

        This is the number of layers in simple circuits. For circuits containing
        sub-circuit blocks, this includes the full depth of these blocks.

        Returns
        -------
        int
        """
        return sum([lbl.depth for lbl in self.layertup])

    @property
    def width(self):
        """
        The circuit width.

        This is the number of qubits on which the circuit acts.

        Returns
        -------
        int
        """
        return len(self._line_labels)

    @property
    def size(self):
        """
        Returns the circuit size.

        This is the sum of the sizes of all the gates in the circuit.

        Returns
        -------
        int
        """
        def label_size(lbl):
            if lbl.IS_SIMPLE:
                return len(lbl.sslbls) if (lbl.sslbls is not None) else len(self._line_labels)
            else:
                return sum([label_size(sublbl) for sublbl in lbl.components])
        return sum([label_size(lbl) for lbl in self.layertup])

    def to_label(self, nreps=1):
        """
        Construct and return this entire circuit as a CircuitLabel.

        Parameters
        ----------
        nreps : int, optional
            The number of times this circuit will be repeated.

        Returns
        -------
        CircuitLabel
        """
        eff_line_labels = None if self._line_labels == ('*',) else self._line_labels
        return _CircuitLabel(self._name, self.layertup, eff_line_labels, nreps)

    def _proc_layers_arg(self, layers):
        """Pre-process the layers argument used by many methods"""
        if layers is None:
            layers = list(range(len(self._labels)))
        elif isinstance(layers, slice):
            if layers.start is None and layers.stop is None:
                layers = list(range(len(self._labels)))
            else:
                layers = _slct.indices(layers, len(self._labels))
        elif not isinstance(layers, (list, tuple)):
            layers = (layers,)
        return layers

    def _proc_lines_arg(self, lines):
        """Pre-process the lines argument used by many methods"""
        if lines is None:
            lines = self._line_labels
        elif isinstance(lines, slice):
            if lines.start is None and lines.stop is None:
                lines = self._line_labels
            else:
                lines = _slct.indices(lines)
        elif not isinstance(lines, (list, tuple)):
            lines = (lines,)
        return lines

    def _proc_key_arg(self, key):
        """Pre-process the key argument used by many methods"""
        if isinstance(key, tuple):
            if len(key) != 2:
                return IndexError("Index must be of the form <layerIndex>,<lineIndex>")
            else:
                return key[0], key[1]
        else:
            return key, None

    @abstractmethod
    def _layer_components(self, ilayer):
        """Get the components of the `ilayer`-th layer as a list/tuple."""
        pass

    def __getitem__(self, key):
        layers, lines = self._proc_key_arg(key)
        return self.extract_labels(layers, lines, strict=True)

    @abstractmethod
    def copy(self, editable='auto'):
        """
        Returns a copy of the circuit.

        Parameters
        ----------
        editable : {True, False, "auto"}
            Whether returned copy is editable.

        Returns
        -------
        BaseCircuit subclass
        """
        pass

    @abstractmethod
    def to_editable(self):
        """
        Convert to an EditableCircuit.

        Returns
        -------
        EditableCircuit
        """
        pass

    @abstractmethod
    def to_static(self):
        """
        Convert to a StaticCircuit.

        Returns
        -------
        StaticCircuit
        """
        pass

    def extract_labels(self, layers=None, lines=None, strict=True):
        """
        Get a subregion - a "rectangle" - of this Circuit.

        Parameters
        ----------
        layers : int, slice, or list/tuple of ints
            Which layers to select.

        lines : str/int, slice, or list/tuple of strs/ints
            Which lines to select.

        strict : bool, optional
            When True, only gates lying completely within the selected
            region are included.

        Returns
        -------
        Label or Circuit subclass
            The requested portion of this circuit.
        """
        nonint_layers = not isinstance(layers, int)
        layers = self._proc_layers_arg(layers)
        lines = self._proc_lines_arg(lines)

        if len(layers) == 0 or len(lines) == 0:
            return self.__class__._fastinit((), tuple(lines), False)

        ret = []
        for i in layers:
            ret_layer = []
            for l in self._layer_components(i):
                sslbls = l.sslbls if isinstance(l, _Label) else \
                         _sslbls_of_nested_lists_of_simple_labels(l)
                if sslbls is None:
                    sslbls = set(self._line_labels)
                else:
                    sslbls = set(sslbls)
                if (strict and sslbls.issubset(lines)) or \
                   (not strict and len(sslbls.intersection(lines)) >= 0):
                    ret_layer.append(l)
            ret.append(_Label(ret_layer) if len(ret_layer) != 1 else ret_layer[0])

        if nonint_layers:
            if not strict: 
                lines = "auto" # since we may have included lbls on other lines
            # don't worry about string rep for now...
            return self.__class__._fastinit(tuple(ret), tuple(lines) if strict else lines, False)
        else:
            return _Label(ret[0])

    @abstractmethod
    def _is_line_idling(self, line_label, idle_layer_labels=None):
        """Whether the line is idling in every circuit layer."""
        pass

    def idling_lines(self, idle_layer_labels=None):
        """
        Returns the line labels corresponding to idling lines.

        Parameters
        ----------
        idle_layer_labels : iterable, optional
            Layer-labels that should be treated as idle operations.

        Returns
        -------
        tuple
        """
        return tuple([x for x in self._line_labels 
                     if self._is_line_idling(x, idle_layer_labels)])

