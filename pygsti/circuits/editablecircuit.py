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
                lines = "auto"
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


class EditableCircuit(BaseCircuit):
    """
    A mutable quantum circuit.

    EditableCircuit represents a quantum circuit that can be modified in place.
    It cannot be hashed and should not be used as dictionary keys. Call to_static()
    or done_editing() to convert it to an immutable StaticCircuit.

    When editable=True, a rich set of operations may be used to construct the
    circuit in place. Layer labels are stored as nested lists of simple labels.
    """

    def __init__(self, layer_labels=(), line_labels='auto', num_lines=None,
                 stringrep=None, name='', check=True, expand_subcircuits="default",
                 occurrence=None, compilable_layer_indices=None):
        """
        Creates a new EditableCircuit object.

        Parameters
        ----------
        layer_labels : iterable of Labels or str
            List of layer labels specifying the gates for the circuit.

        line_labels : iterable, optional
            The label for each circuit line.

        num_lines : int, optional
            Specify instead of line_labels to use integers 0 to num_lines-1.

        stringrep : string, optional
            A string representation for the circuit.

        name : str, optional
            A name for this circuit.

        check : bool, optional
            Whether to check consistency.

        expand_subcircuits : bool or "default"
            Whether to expand sub-circuits.

        occurrence : hashable, optional
            An occurrence id for this circuit.

        compilable_layer_indices : tuple, optional
            Circuit-layer indices that may be compiled.
        """
        super().__init__()
        from pygsti.circuits.circuitparser import CircuitParser as _CircuitParser

        # Handle various input types for layer_labels
        layer_labels_objs = None
        if isinstance(layer_labels, str):
            cparser = _CircuitParser()
            layer_labels_objs, line_labels_from_parser, occurrence_from_parser, compilable_from_parser = \
                cparser.parse(layer_labels)
            if line_labels == 'auto' and line_labels_from_parser is not None:
                line_labels = line_labels_from_parser
            if occurrence is None and occurrence_from_parser is not None:
                occurrence = occurrence_from_parser
            if compilable_layer_indices is None and compilable_from_parser is not None:
                compilable_layer_indices = compilable_from_parser
            layer_labels = layer_labels_objs

        if expand_subcircuits == "default":
            expand_subcircuits = self.default_expand_subcircuits
        if expand_subcircuits and layer_labels is not None:
            layer_labels = _expand_subcircuits_in_layer_labels(layer_labels)

        # Parse stringrep if needed
        if stringrep is not None and (layer_labels is None or check):
            cparser = _CircuitParser()
            chk_labels, chk_line_labels, chk_occurrence, chk_compilable = cparser.parse(stringrep)
            if layer_labels is None:
                layer_labels = chk_labels

        if layer_labels is None:
            layer_labels = ()

        # Set line_labels
        if line_labels == 'auto':
            explicit_sslbls = _accumulate_explicit_sslbls(layer_labels)
            my_line_labels = tuple(sorted(explicit_sslbls)) if len(explicit_sslbls) > 0 else ('*',)
        else:
            my_line_labels = tuple(line_labels) if line_labels is not None else ('*',)

        if (num_lines is not None) and (num_lines != len(my_line_labels)):
            if len(my_line_labels) == 1 and my_line_labels[0] == '*':
                my_line_labels = tuple(range(num_lines))
            else:
                raise ValueError("Conflicting num_lines and line_labels arguments!")

        # Convert layer_labels to nested lists of simple labels
        labels = []
        for layer_lbl in layer_labels:
            labels.append(_label_to_nested_lists_of_simple_labels(layer_lbl, my_line_labels))

        # Set compilable indices
        compilable_layer_indices_tup = tuple(compilable_layer_indices) if compilable_layer_indices else ()

        # Initialize
        self._labels = labels
        self._line_labels = my_line_labels
        self._occurrence_id = occurrence
        self._compilable_layer_indices_tup = compilable_layer_indices_tup
        self._str = None  # Lazy generation
        self._name = name
        self.auxinfo = {}

    @classmethod
    def from_tuple(cls, tup):
        """
        Creates an EditableCircuit from a tuple.

        Parameters
        ----------
        tup : tuple
            The tuple to convert.

        Returns
        -------
        EditableCircuit
        """
        if '@' in tup:
            k = tup.index('@')
            return cls(tup[0:k], tup[k + 1:])
        else:
            return cls(tup)

    @classmethod
    def _fastinit(cls, labels, line_labels, editable, name='', stringrep=None, occurrence=None,
                  compilable_layer_indices_tup=()):
        """Fast initialization for internal use."""
        ret = cls.__new__(cls)
        ret._labels = list(labels) if not isinstance(labels, list) else labels
        ret._line_labels = tuple(line_labels)
        ret._occurrence_id = occurrence
        ret._compilable_layer_indices_tup = compilable_layer_indices_tup
        ret._str = stringrep
        ret._name = name
        ret.auxinfo = {}
        return ret

    @property
    def layertup(self):
        """This Circuit's layers as a standard Python tuple of layer Labels."""
        return tuple([layer_lbl if isinstance(layer_lbl, _Label) else _Label(layer_lbl) 
                     for layer_lbl in self._labels])

    @property
    def str(self):
        """The Python string representation of this Circuit."""
        if self._str is None:
            self._str = _op_seq_to_str(self._labels, self._line_labels, 
                                      self._occurrence_id, self._compilable_layer_indices_tup)
        return self._str

    @str.setter
    def str(self, value):
        """Set the Python string representation of this Circuit."""
        from pygsti.circuits.circuitparser import CircuitParser as _CircuitParser
        cparser = _CircuitParser()
        chk, chk_labels, chk_occurrence, chk_compilable_inds = cparser.parse(value)

        if not all([my_layer in (chk_lbl, [chk_lbl]) for chk_lbl, my_layer in zip(chk, self._labels)]):
            raise ValueError("String representation is inconsistent with circuit layers!")
        
        self._str = value

    @property
    def line_labels(self):
        """The line labels (often qubit labels) of this circuit."""
        return self._line_labels

    @line_labels.setter
    def line_labels(self, value):
        """Set the line labels."""
        if value == self._line_labels:
            return
        removed_line_labels = set(self._line_labels) - set(value)
        if removed_line_labels:
            raise ValueError("Cannot remove line labels %s - delete lines first!" % str(removed_line_labels))
        self._line_labels = tuple(value)
        self._str = None

    @property
    def occurrence(self):
        """The occurrence id of this circuit."""
        return self._occurrence_id

    @occurrence.setter
    def occurrence(self, value):
        """Set the occurrence id."""
        self._occurrence_id = value
        self._str = None

    @property
    def compilable_layer_indices(self):
        """Tuple of the layer indices corresponding to "compilable" layers."""
        return self._compilable_layer_indices_tup

    @compilable_layer_indices.setter
    def compilable_layer_indices(self, val):
        """Set the compilable layer indices."""
        self._compilable_layer_indices_tup = tuple(val) if (val is not None) else ()

    def __hash__(self):
        raise TypeError("EditableCircuit is not hashable. Call done_editing() or to_static() first.")

    def _layer_components(self, ilayer):
        """Get the components of the `ilayer`-th layer as a list/tuple."""
        return self._labels[ilayer] if isinstance(self._labels[ilayer], list) else [self._labels[ilayer]]

    def _remove_layer_component(self, ilayer, indx):
        """Removes the `indx`-th component from the `ilayer`-th layer"""
        if isinstance(self._labels[ilayer], list):
            del self._labels[ilayer][indx]
        else:
            assert(indx == 0), "Only index 0 exists for a single-simple-Label level"
            self._labels[ilayer] = []

    def _append_layer_component(self, ilayer, val):
        """Add `val` to the `ilayer`-th layer"""
        if isinstance(self._labels[ilayer], list):
            self._labels[ilayer].append(val)
        else:
            self._labels[ilayer] = [self._labels[ilayer], val]

    def _replace_layer_component(self, ilayer, indx, val):
        """Replace `indx`-th component of `ilayer`-th layer with `val`"""
        if isinstance(self._labels[ilayer], list):
            self._labels[ilayer][indx] = val
        else:
            assert(indx == 0), "Only index 0 exists for a single-simple-Label level"
            self._labels[ilayer] = val

    def copy(self, editable='auto'):
        """
        Returns a copy of the circuit.

        Parameters
        ----------
        editable : {True, False, "auto"}
            Whether returned copy is editable. If "auto", copy is editable.

        Returns
        -------
        EditableCircuit or StaticCircuit
        """
        if editable == "auto":
            editable = True

        if editable:
            ret = EditableCircuit.__new__(EditableCircuit)
            labels_copy = [layer[:] if isinstance(layer, list) else layer for layer in self._labels]
            ret._labels = labels_copy
            ret._line_labels = self._line_labels
            ret._occurrence_id = self._occurrence_id
            ret._compilable_layer_indices_tup = self._compilable_layer_indices_tup
            ret._str = None
            ret._name = self._name
            ret.auxinfo = self.auxinfo.copy()
            return ret
        else:
            return self.to_static()

    def to_editable(self):
        """Convert to an EditableCircuit (returns self since already editable)."""
        return self

    def to_static(self):
        """
        Convert to a StaticCircuit (immutable version).

        Returns
        -------
        StaticCircuit
        """
        static_labels = tuple([layer_lbl if isinstance(layer_lbl, _Label) else _Label(layer_lbl) 
                              for layer_lbl in self._labels])
        return StaticCircuit._fastinit(static_labels, self._line_labels, False,
                                      self._name, self._str, self._occurrence_id,
                                      self._compilable_layer_indices_tup)

    def done_editing(self):
        """
        Convert this EditableCircuit to a StaticCircuit.

        Returns
        -------
        StaticCircuit
        """
        return self.to_static()

    def clear(self):
        """Removes all the gates in the circuit (preserving the number of lines)."""
        self._labels = []

    def __setitem__(self, key, val):
        layers, lines = self._proc_key_arg(key)
        return self.set_labels(val, layers, lines)

    def __delitem__(self, key):
        layers, lines = self._proc_key_arg(key)
        if layers is None:
            self.delete_lines(lines)
        elif lines is None:
            self.delete_layers(layers)
        else:
            self.clear_labels(layers, lines)

    def set_labels(self, lbls, layers=None, lines=None):
        """
        Write `lbls` to the block defined by the `layers` and `lines` arguments.

        Parameters
        ----------
        lbls : Label, list/tuple of Labels, or Circuit
            Labels to write to the specified region.

        layers : int, slice, or list/tuple of ints
            Which layers to set.

        lines : str/int, slice, or list/tuple of strs/ints
            Which lines to set.

        Returns
        -------
        None
        """
        all_layers = bool(layers is None)
        int_layers = isinstance(layers, int)
        layers = self._proc_layers_arg(layers)

        all_lines = bool(lines is None)
        lines = self._proc_lines_arg(lines)

        # Convert lbls to appropriate format
        if int_layers:
            if isinstance(lbls, BaseCircuit):
                lbls = lbls.to_label()
            lbls = to_label(lbls)
            lbls_sslbls = None if (lbls.sslbls is None) else set(lbls.sslbls)
        else:
            if isinstance(lbls, BaseCircuit):
                lbls = lbls.layertup
            lbls = tuple(map(to_label, lbls))
            lbls_sslbls = None if any([l.sslbls is None for l in lbls]) \
                else set(_itertools.chain(*[l.sslbls for l in lbls]))

        if len(layers) == 0 or len(lines) == 0: 
            return

        # Handle layer expansion if needed
        if all_layers:
            while len(lbls) > len(self._labels):
                self._labels.append([])
        elif len(layers) > 1:
            assert(len(layers) == len(lbls)), \
                "Block width mismatch: assigning %d layers to %d layers" % (len(lbls), len(layers))

        # Handle new line labels
        if lbls_sslbls is not None:
            new_line_labels = set(lbls_sslbls) - set(self._line_labels)
            if all_lines and len(new_line_labels) > 0:
                self._line_labels = tuple(sorted(set(self._line_labels) | new_line_labels))

        # Remove all labels in block to be assigned
        self._clear_labels(layers, lines)

        def_sslbls = None if all_lines else lines
        if not int_layers:
            for i, lbls_comp in zip(layers, lbls):
                self._labels[i].extend(_label_to_nested_lists_of_simple_labels(lbls_comp, def_sslbls))
        else:
            self._labels[layers[0]].extend(_label_to_nested_lists_of_simple_labels(lbls, def_sslbls))

    def _clear_labels(self, layers, lines, clear_straddlers=False):
        """Remove all labels in a block given by layers and lines."""
        for i in layers:
            inds_to_delete = []
            for k, lbl in enumerate(self._layer_components(i)):
                sslbls = _sslbls_of_nested_lists_of_simple_labels(lbl)
                if sslbls is None:
                    sslbls = set(self._line_labels)
                else:
                    sslbls = set(sslbls)
                if clear_straddlers and len(sslbls.intersection(lines)) > 0:
                    inds_to_delete.append(k)
                elif sslbls.issubset(lines):
                    inds_to_delete.append(k)
            for k in reversed(inds_to_delete):
                self._remove_layer_component(i, k)

    def clear_labels(self, layers=None, lines=None, clear_straddlers=False):
        """Removes all the gates within the given circuit region."""
        layers = self._proc_layers_arg(layers)
        lines = self._proc_lines_arg(lines)
        self._clear_labels(layers, lines, clear_straddlers)

    def delete_layers(self, layers=None):
        """Deletes one or more layers from the circuit."""
        layers = self._proc_layers_arg(layers)
        for i in reversed(sorted(layers)):
            del self._labels[i]

    def delete_lines(self, lines, delete_straddlers=False):
        """Deletes one or more lines from the circuit."""
        lines = self._proc_lines_arg(lines)
        for i in range(len(self._labels)):
            inds_to_delete = []
            for k, lbl in enumerate(self._layer_components(i)):
                sslbls = _sslbls_of_nested_lists_of_simple_labels(lbl)
                if sslbls is None:
                    sslbls = set(self._line_labels)
                else:
                    sslbls = set(sslbls)
                if delete_straddlers and len(sslbls.intersection(lines)) > 0:
                    inds_to_delete.append(k)
                elif sslbls.issubset(lines):
                    inds_to_delete.append(k)
            for k in reversed(inds_to_delete):
                self._remove_layer_component(i, k)
        self._line_labels = tuple([x for x in self._line_labels if x not in lines])

    def _is_line_idling(self, line_label, idle_layer_labels=None):
        """Whether the line is idling in every circuit layer."""
        all_sslbls = _sslbls_of_nested_lists_of_simple_labels(self._labels, idle_layer_labels)
        if all_sslbls is None:
            return False
        return bool(line_label not in all_sslbls)

    # Additional methods for circuit manipulation
    def insert_layer_inplace(self, circuit_layer, j):
        """Insert a layer at position j."""
        if j is None: 
            j = len(self._labels)
        elif j < 0: 
            j = len(self._labels) + j
        layer_as_nested_lists = _label_to_nested_lists_of_simple_labels(circuit_layer, self._line_labels)
        self._labels.insert(j, layer_as_nested_lists)
        self._str = None

    def append_circuit_inplace(self, circuit):
        """Append another circuit to this one."""
        if isinstance(circuit, BaseCircuit):
            for layer in circuit._labels:
                self._labels.append(layer if isinstance(layer, list) else [layer])
        else:
            raise ValueError("Can only append circuits")
        self._str = None


class StaticCircuit(BaseCircuit):
    """
    An immutable quantum circuit.

    StaticCircuit represents a quantum circuit that cannot be modified.
    It can be hashed and used as dictionary keys. Layer labels are stored
    as a tuple of Label objects.

    To modify a StaticCircuit, you must create a new modified copy using
    methods that return new circuits, or convert to EditableCircuit first.
    """

    def __init__(self, layer_labels=(), line_labels='auto', num_lines=None,
                 stringrep=None, name='', check=True, expand_subcircuits="default",
                 occurrence=None, compilable_layer_indices=None):
        """
        Creates a new StaticCircuit object.

        Parameters are the same as EditableCircuit, but the resulting circuit
        is immutable and hashable.
        """
        super().__init__()
        from pygsti.circuits.circuitparser import CircuitParser as _CircuitParser

        # Similar initialization as EditableCircuit but create static labels
        layer_labels_objs = None
        if isinstance(layer_labels, str):
            cparser = _CircuitParser()
            layer_labels_objs, line_labels_from_parser, occurrence_from_parser, compilable_from_parser = \
                cparser.parse(layer_labels)
            if line_labels == 'auto' and line_labels_from_parser is not None:
                line_labels = line_labels_from_parser
            if occurrence is None and occurrence_from_parser is not None:
                occurrence = occurrence_from_parser
            if compilable_layer_indices is None and compilable_from_parser is not None:
                compilable_layer_indices = compilable_from_parser
            layer_labels = layer_labels_objs

        if expand_subcircuits == "default":
            expand_subcircuits = self.default_expand_subcircuits
        if expand_subcircuits and layer_labels is not None:
            layer_labels = _expand_subcircuits_in_layer_labels(layer_labels)

        if layer_labels is None:
            layer_labels = ()

        # Set line_labels
        if line_labels == 'auto':
            explicit_sslbls = _accumulate_explicit_sslbls(layer_labels)
            my_line_labels = tuple(sorted(explicit_sslbls)) if len(explicit_sslbls) > 0 else ('*',)
        else:
            my_line_labels = tuple(line_labels) if line_labels is not None else ('*',)

        if (num_lines is not None) and (num_lines != len(my_line_labels)):
            if len(my_line_labels) == 1 and my_line_labels[0] == '*':
                my_line_labels = tuple(range(num_lines))
            else:
                raise ValueError("Conflicting num_lines and line_labels arguments!")

        # Convert to static tuple of Labels
        labels = tuple([to_label(lbl) for lbl in layer_labels])

        compilable_layer_indices_tup = tuple(compilable_layer_indices) if compilable_layer_indices else ()

        # Initialize
        self._labels = labels
        self._line_labels = my_line_labels
        self._occurrence_id = occurrence
        self._compilable_layer_indices_tup = compilable_layer_indices_tup
        self._str = stringrep
        self._name = name
        self.auxinfo = {}
        
        # Compute hash
        self._hashable_tup = self.tup
        self._hash = hash(self._hashable_tup)

    @classmethod
    def from_tuple(cls, tup):
        """Creates a StaticCircuit from a tuple."""
        if '@' in tup:
            k = tup.index('@')
            return cls(tup[0:k], tup[k + 1:])
        else:
            return cls(tup)

    @classmethod
    def _fastinit(cls, labels, line_labels, editable, name='', stringrep=None, occurrence=None,
                  compilable_layer_indices_tup=()):
        """Fast initialization for internal use."""
        ret = cls.__new__(cls)
        ret._labels = tuple(labels) if not isinstance(labels, tuple) else labels
        ret._line_labels = tuple(line_labels)
        ret._occurrence_id = occurrence
        ret._compilable_layer_indices_tup = compilable_layer_indices_tup
        ret._str = stringrep
        ret._name = name
        ret.auxinfo = {}
        ret._hashable_tup = ret.tup
        ret._hash = hash(ret._hashable_tup)
        return ret

    @property
    def layertup(self):
        """This Circuit's layers as a standard Python tuple of layer Labels."""
        return self._labels

    @property
    def str(self):
        """The Python string representation of this Circuit."""
        if self._str is None:
            self._str = _op_seq_to_str(self._labels, self._line_labels, 
                                      self._occurrence_id, self._compilable_layer_indices_tup)
        return self._str

    def __hash__(self):
        return self._hash

    def __getstate__(self):
        """For pickling."""
        return self.__dict__

    def __setstate__(self, state_dict):
        """For unpickling."""
        for k, v in state_dict.items():
            self.__dict__[k] = v
        if not hasattr(self, '_hash') or self._hash is None:
            if hasattr(self, '_hashable_tup') and self._hashable_tup is not None:
                self._hash = hash(self._hashable_tup)
            else:
                self._hashable_tup = self.tup
                self._hash = hash(self._hashable_tup)

    def _layer_components(self, ilayer):
        """Get the components of the `ilayer`-th layer as a list/tuple."""
        if self._labels[ilayer].IS_SIMPLE: 
            return [self._labels[ilayer]]
        else: 
            return self._labels[ilayer].components

    def copy(self, editable='auto'):
        """
        Returns a copy of the circuit.

        Parameters
        ----------
        editable : {True, False, "auto"}
            Whether returned copy is editable. If "auto", copy is static.

        Returns
        -------
        StaticCircuit or EditableCircuit
        """
        if editable == "auto":
            editable = False

        if editable:
            return self.to_editable()
        else:
            # Static copy - can share immutable data
            ret = StaticCircuit.__new__(StaticCircuit)
            ret._labels = self._labels
            ret._line_labels = self._line_labels
            ret._occurrence_id = self._occurrence_id
            ret._compilable_layer_indices_tup = self._compilable_layer_indices_tup
            ret._str = self._str
            ret._name = self._name
            ret.auxinfo = self.auxinfo.copy()
            ret._hashable_tup = self._hashable_tup
            ret._hash = self._hash
            return ret

    def to_editable(self):
        """
        Convert to an EditableCircuit (mutable version).

        Returns
        -------
        EditableCircuit
        """
        labels_as_lists = [_label_to_nested_lists_of_simple_labels(lbl, self._line_labels) 
                          for lbl in self._labels]
        return EditableCircuit._fastinit(labels_as_lists, self._line_labels, True,
                                        self._name, None, self._occurrence_id,
                                        self._compilable_layer_indices_tup)

    def to_static(self):
        """Convert to a StaticCircuit (returns self since already static)."""
        return self

    def _is_line_idling(self, line_label, idle_layer_labels=None):
        """Whether the line is idling in every circuit layer."""
        layers = [x for x in self._labels if x not in idle_layer_labels] if idle_layer_labels else self._labels
        all_sslbls = None if any([layer.sslbls is None for layer in layers]) \
            else set([sslbl for layer in layers for sslbl in layer.sslbls])
        if all_sslbls is None:
            return False
            return bool(line_label not in all_sslbls)

    # Methods that return modified copies (no in-place modification allowed)
    def insert_layer(self, circuit_layer, j):
        """Insert a layer at position j, returning a new StaticCircuit."""
        editable = self.to_editable()
        editable.insert_layer_inplace(circuit_layer, j)
        return editable.to_static()

    def append_circuit(self, circuit):
        """Append another circuit, returning a new StaticCircuit."""
        editable = self.to_editable()
        editable.append_circuit_inplace(circuit)
        return editable.to_static()

    def delete_layers(self, layers=None):
        """Delete layers, returning a new StaticCircuit."""
        editable = self.to_editable()
        editable.delete_layers(layers)
        return editable.to_static()

    def delete_lines(self, lines, delete_straddlers=False):
        """Delete lines, returning a new StaticCircuit."""
        editable = self.to_editable()
        editable.delete_lines(lines, delete_straddlers)
        return editable.to_static()

    def clear_labels(self, layers=None, lines=None, clear_straddlers=False):
        """Clear labels, returning a new StaticCircuit."""
        editable = self.to_editable()
        editable.clear_labels(layers, lines, clear_straddlers)
        return editable.to_static()


# Helper function for expanding subcircuits
def _expand_subcircuits_in_layer_labels(layer_labels):
    """Expand any subcircuits in layer_labels."""
    # This is a placeholder - the actual implementation would need to
    # handle CircuitLabel expansion logic from the original Circuit class
    return layer_labels


# Maintain backward compatibility by aliasing Circuit to EditableCircuit
# (or you could make Circuit a factory that returns the appropriate type)
Circuit = EditableCircuit
