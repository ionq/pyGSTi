"""
Random compilation routines
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

from pygsti.baseobjs.nicelyserializable import NicelySerializable as _NicelySerializable


class LayerRandomCompilationFunction(_NicelySerializable):
    """
    A convenient base class for building serializable "functions" for the
    random compilation of layers in a circuit.
    """
    def __call__(self, layer, preframe, postframe):
        """Recompiles a layer based on random frames of reference sandwiching the layer.

        Paramaters:
        -----------
        layer: list
            A Circuit layer to recompile.

        preframe: list or dict
            Contains the information describing the random operation before the layer
            for each line_label in the Circuit.

        postframe: list or dict
            Contains the information describing the random operation after the layer
            for each line_label in the Circuit.
        
        Returns
        -------
        rc_layer: list
            The recompiled layer.
        """
        raise NotImplementedError("Derived classes should implement this!")

    def __init__(self):
        super().__init__()

class U3LayerPauliRCFunction(LayerRandomCompilationFunction):
    """A rule performing Pauli random compiling on layers with U3, CNOT, and/or CPHASE gates.

    The Gu3 label matches Qiskit U3 convention with theta, phi, and lambda args such that:
        Gu3(\theta, \phi, \lambda) = \cos(\theta/2)          & -e^{i\lambda}\sin(\theta/2)\\
                                     e^{i\phi}\sin(\theta/2) & e^{i(\phi+\lambda)\cos(\theta/2)
    Alternatively, this is also (in matmul order) RZ(\phi + \pi) SX RZ(\theta + \pi) SX RZ(\lambda)
    using a ZXZXZ decomposition with Z rotations (RZ) and sqrt(X) (SX) gates.

    The pre/post frames are given as a Nx2 array for N qubits, where the first and second column
    indicate Z or X components respectively, e.g. I is 00, X is 01, Y is 11, and Z is 10.
    """
    def __call__(self, layer, preframe, postframe):


# TODO: Port over/clean up random compilation code
class RandomCompiler(object):
    """A prescription for randomly compiling the gates of one circuit into other gates.

    Parameters
    ----------
    layer_comp_fn: LayerRandomCompilationFunction
        A function which takes a l
    """

