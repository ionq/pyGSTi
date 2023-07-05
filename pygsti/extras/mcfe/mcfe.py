"""
Mirror Circuit Fidelity Estimation objects
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

from pygsti.circuits import Circuit as _Circuit
from pygsti.protocols.protocol import CircuitListsDesign as _CListDesign
from pygsti.protocols import vb as _vb

class MCFEDesign(_vb.BenchmarkingDesign):
    """Experiment design that holds circuits for mirror circuit fidelity estimation.

    Parameters
    ----------
    """
    def __init__(self, circuits_or_edesign, num_mirror_samples=300, num_ref_samples=100,
                 forward_compilation_fn=None, mirror_compilation_fn=None, mirror_compilation_src='original'):
        if isinstance(circuits_or_edesign, _Circuit):
            base_edesign = _CListDesign([circuits_or_edesign])
        elif not isinstance(circuits_or_edesign, _CListDesign):
            base_edesign = _CListDesign(circuits_or_edesign)
        else:
            base_edesign = circuits_or_edesign

        # TODO: Build circuit mirroring edesign

# TODO: Specific subclasses

class MCFESummaryStatistics(_vb.SummaryStatistics):
    """Extension of summary statistics for mirror circuit fidelity estimation
    """
    # TODO: Add MCFE-specific things to *_statistics members
    def __init__(self, name):
        super().__init__(name)

