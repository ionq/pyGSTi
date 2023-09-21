"""
Random compilation routines for U3+CNOT layers by Pauli frames
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

from pygsti.baseobjs import Label as _Label
from pygsti.extras.randomcompiler.randomcompiler import RandomCompiler as _RandomCompiler


class U3LayerPauliFrameRC(_RandomCompiler):
    """A random compiler for U3 circuits being randomized by Pauli frames.

    The circuits should consist of layers that are either:
        - Only 'Gu3' gates (described below)
        - Only CNOT or CPHASE gates

    The Gu3 label should match Qiskit U3 convention with theta, phi, and lambda args such that:
        Gu3(\theta, \phi, \lambda) = \cos(\theta/2)          & -e^{i\lambda}\sin(\theta/2)\\
                                     e^{i\phi}\sin(\theta/2) & e^{i(\phi+\lambda)\cos(\theta/2)
    Alternatively, this is also (in matmul order) RZ(\phi + \pi) SX RZ(\theta + \pi) SX RZ(\lambda)
    using a ZXZXZ decomposition with Z rotations (RZ) and sqrt(X) (SX) gates.

    The pre/post frames are given as Nx2 array where each row corresponds to the symplectic phase
    vector of the Pauli operation, i.e. I = [0,0], X=[0,2], Y=[2,2], Z=[2,0].
    """
    def get_default_frame(self, layer, line_labels):
        """Definition of a default frame to use given a padded layer.

        For Pauli frames, this is a layer-independent Nx2 array of zeros.

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
        return _np.zeros((len(line_labels), 2))

    def get_next_frame(self, layer, line_labels, preframe):
        """Compute or sample the next frame given a layer and preframe.

        For a layer of U3 gates, this is a randomly sampled set of symplectic
        phase vectors for each qubit.
        For a layer of two-qubit gates, the preframe is commuted through the
        respective two-qubit gate.

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
        num_qubits = len(line_labels)
        label_names = [c.name for c in layer.components]
        if all([ln == 'Gu3' for ln in label_names]):
            # Randomly sample a frame
            postframe = 2*self.rand_state.randint(0, 2, (num_qubits, 2))
        elif all([ln in ['Gcnot'] for ln in label_names]): # TODO: Add Gcphase in here too
            # Pass the frame through the 2-qubit gate
            postframe = preframe
            for comp in layer:
                icontrol, itarget = [line_labels.index(q) for q in comp.qubits]
                postframe[icontrol, :] = (preframe[icontrol, :] + preframe[itarget, :]) % 4
        else:
            raise RuntimeError("Circuit must consist on only Gu3 or Gcnot layers") # TODO: Gcphase
        
        return postframe
    
    def layer_compilation(self, layer, line_labels, preframe, postframe):
        """Recompile a layer based on random frames sandwiching the layer.

        For a layer of U3 gates, the U3 parameters are adjusted based on the
        surrounding frames. For a layer of two-qubit gates, the gate is passed
        through because the postframe has already been updated.

        Paramaters:
        -----------
        layer: list
            A Circuit layer to recompile.

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
        used_qubits = []

        new_layer = []

        def mod_2pi(theta):
            while (theta > _np.pi or theta <= -1 * _np.pi):
                if theta > _np.pi:
                    theta = theta - 2 * _np.pi
                elif theta <= -1 * _np.pi:
                    theta = theta + 2 * _np.pi
            return theta

        for comp in layer:
            if comp.name == 'Gu3':
                # Update U3 parameters based on frames
                theta, phi, lamb = comp.args
                iqubit = line_labels.index(comp.qubits[0])
                if preframe[iqubit][0] == 2:    # Z gate preceeding the layer
                    lamb = lamb + _np.pi
                if postframe[iqubit][0] == 2:   # Z gate following the layer
                    phi = phi + _np.pi
                if preframe[iqubit][1] == 2:    # X gate preceeding the layer
                    theta = theta - _np.pi
                    phi = phi
                    lamb = -lamb - _np.pi
                if postframe[iqubit][1] == 2:   # X gate following the layer
                    theta = theta - _np.pi
                    phi = -phi - _np.pi
                    lamb = lamb

                new_args = (mod_2pi(theta), mod_2pi(phi), mod_2pi(lamb))
                new_label = _Label('Gu3', comp.qubits, args=new_args)
                new_layer.append(new_label)
            else:
                # This is a 2Q gate that only updates the frames
                new_layer.append(comp)
            used_qubits.extend(comp.qubits)

        assert set(used_qubits) == set(line_labels), \
            "Random compilation must be performed on a idle-padded layer"

        return new_layer

    