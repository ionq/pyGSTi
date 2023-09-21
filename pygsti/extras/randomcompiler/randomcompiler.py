"""
Base random compilation routines
"""
#***************************************************************************************************
# Copyright 2015, 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#***************************************************************************************************

import numpy as _np

from pygsti.baseobjs import StateSpace as _StateSpace
from pygsti.baseobjs.nicelyserializable import NicelySerializable as _NicelySerializable
from pygsti.circuits import Circuit as _Circuit


class RandomCompiler(_NicelySerializable):
    """A prescription for randomly compiling the gates of one circuit into other gates.

    Parameters
    ----------
    idle_lbl: Label
        A label to insert for idle operations during compilation.
        Value of sslbls is a dummy and will be replaced during padding.

    seed: int, RandomState, or None
        A seed or random state used for sampling random frames
    """
    def __init__(self, idle_lbl, seed=None):
        self.idle_lbl = idle_lbl

        if not isinstance(seed, _np.random.RandomState):
            self.rand_state = _np.random.RandomState(seed)
        else:
            self.rand_state = seed
    
    def random_compile_circuit(self, circuit, random_first_frame=False,
                               random_last_frame=False, return_frames=False):
        """Recompile a circuit.

        Parameters
        ----------
        circuit: Circuit
            Initial circuit to randomly compile
        
        random_first_frame: bool
            Whether to sample the first frame randomly (True) or use an all-zeros
            initial frame (False, the default).
        
        random_last_frame: bool
            Whether to sample the last frame randomly (True) or use an all-zeros
            final frame (False, the default).

        return_frames: bool
            Whether to return the first and last frames

        Returns
        -------
        rc_circuit: Circuit
            Randomly compiled circuit
        
        first_frame: Any
            First random frame. Only returned when return_frames is True

        last_frame: Any
            Final random frame. Only returned when return_frames is True
        """
        line_labels = circuit.line_labels
        rc_circuit = _Circuit(line_labels=line_labels, editable=True)

        # Initialize first frame
        padded_layer = self._pad_layer(circuit.layer_label(0).components, line_labels)
        first_frame = self.get_default_frame(padded_layer, line_labels)
        if random_first_frame:
            first_frame = self.get_next_frame(padded_layer, line_labels, first_frame)
        
        # Run through every layer, pulling a new postframe and performing RC
        preframe = first_frame
        for i in range(circuit.depth):
            padded_layer = self._pad_layer(circuit.layer_label(i).components, line_labels)

            postframe = self.get_next_frame(padded_layer, line_labels, preframe)
            if i == circuit.depth-1 and not random_last_frame:
                postframe = self.get_default_frame(padded_layer, line_labels)
            
            new_layer = self.layer_compilation(padded_layer, line_labels, preframe, postframe)
            
            rc_circuit.insert_layer_inplace(new_layer, i)

            preframe = postframe
        last_frame = postframe

        rc_circuit.done_editing()

        if return_frames:
            return rc_circuit, first_frame, last_frame
        
        return rc_circuit
    
    def get_default_frame(self, layer, line_labels):
        """Definition of a default frame to use given a padded layer.

        Parameters
        ----------
        layer: list
            A Circuit layer to recompile. Should be padded with idles such
            that every line has an operation.
                    
        line_labels: list
            Line labels in the Circuit.

        Returns
        -------
        frame: Any
            The default frame for a given layer.
            Format is specified by the derived class.
        """
        raise NotImplementedError("Derived classes should implement this!")

    def get_next_frame(self, layer, line_labels, preframe):
        """Compute or sample the next frame given a layer and preframe.

        Parameters
        ----------
        layer: list
            A Circuit layer to recompile.
        
        line_labels: list
            Line labels in the Circuit.
        
        preframe: Any
            The random frame preceding the given layer.
            Format is specified by the derived class.
        
        Returns
        -------
        postframe: Any
            The random frame following the given layer.
            Format is specified by derived class.
        """
        raise NotImplementedError("Derived classes should implement this!")
    
    def layer_compilation(self, layer, line_labels, preframe, postframe):
        """Recompile a layer based on random frames sandwiching the layer.

        Parameters
        ----------
        layer: list
            A Circuit layer to recompile.
        
        line_labels: list
            Line labels in the Circuit.

        preframe: Any
            The random frame preceding the given layer.
            Format is specified by the derived class.
        
        postframe: Any
            The random frame following the given layer.
            Format is specified by derived class.
        
        Returns
        -------
        rc_layer: list
            The recompiled layer.
        """
        raise NotImplementedError("Derived classes should implement this!")

    def _pad_layer(self, layer, line_labels):
        used_qubits = set()
        for comp in layer.components:
            used_qubits += set(comp.sslbls)
        
        idle_qubits = set(line_labels) - used_qubits
        
        padded_layer = layer
        for q in idle_qubits:
            idle_lbl = self.idle_lbl
            idle_lbl.sslbls = (q,)
            padded_layer.append(idle_lbl)
        
        return padded_layer

    def _to_nice_serialization(self):
        state = super()._to_nice_serialization()

        state['idle_lbl'] = self.idle_lbl

        rand_state = list(self.rand_state.get_state())
        rand_state[1] = rand_state[1].tolist()
        state['rand_state'] = rand_state

        return state

    @classmethod
    def _from_nice_serialization(cls, state):
        rand_state = _np.random.RandomState()
        rand_state.set_state(state['rand_state'])
        return cls(state['idle_lbl'], seed=rand_state)

