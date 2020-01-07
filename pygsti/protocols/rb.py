""" RB Protocol objects """
#***************************************************************************************************
# Copyright 2015, 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#***************************************************************************************************

import time as _time
import os as _os
import numpy as _np
import pickle as _pickle
import collections as _collections
import warnings as _warnings
import copy as _copy
import scipy.optimize as _spo
from scipy.stats import chi2 as _chi2

from . import protocol as _proto
from .modeltest import ModelTest as _ModelTest
from .. import objects as _objs
from .. import algorithms as _alg
from .. import construction as _construction
from .. import io as _io
from .. import tools as _tools

from ..objects import wildcardbudget as _wild
from ..objects.profiler import DummyProfiler as _DummyProfiler
from ..objects import objectivefns as _objfns

from ..extras import rb as _rb


#Useful to have a base class?
#class RBInput(_proto.ProtocolInput):
#    pass

#Structure:
# MultiInput -> specifies multiple circuit structures on (possibly subsets of) the same data (e.g. collecting into one large dataset the data for multiple protocols)
# MultiProtocol -> runs, on same input circuit structure & data, multiple protocols (e.g. different GST & model tests on same GST data)
#   if that one input is a MultiInput, then it must have the same number of inputs as there are protocols and each protocol is run on the corresponding input.
#   if that one input is a normal input, then the protocols can cache information in a Results object that is handed down.
# SimultaneousInput -- specifies a "qubit structure" for each sub-input
# SimultaneousProtocol -> runs multiple protocols on the same data, but "trims" circuits and data before running sub-protocols
#  (e.g. Volumetric or randomized benchmarks on different subsets of qubits) -- only accepts SimultaneousInputs.

#Inputs:
# Simultaneous: (spec1)
#    Q1: ByDepthData
#    Q2,Q3: ByDepthData
#    Q4: ByDepthData

#Protocols:
# Simultaneous:
#    Q1: Multi
#      VB (aggregate)
#      PredictedModelA
#      PredictedModelB
#      Datasetcomp (between passes)
#    Q2,Q3: VB
#    Q4: VB

#OR: so that auto-running performs the above protocols:
#Inputs:
# Simultaneous: (spec1)
#    Q1: MultiInput(VB, PredictedModelA, PredictedModelB) - or MultiBenchmark?
#       ByDepthData
#    Q2,Q3: VB-ByDepthData
#    Q4: VB-ByDepthData


class ByDepthInput(_proto.CircuitListsInput):
    def __init__(self, depths, circuit_lists, qubit_labels=None):
        assert(len(depths) == len(circuit_lists)), \
            "Number of depths must equal the number of circuit lists!"
        super().__init__(circuit_lists, qubit_labels=qubit_labels)
        self.depths = depths

        
class BenchmarkingInput(ByDepthInput):
    def __init__(self, depths, circuit_lists, ideal_outs, qubit_labels=None):
        assert(len(depths) == len(ideal_outs))
        super().__init__(depths, circuit_lists, qubit_labels)
        self.idealout_lists = ideal_outs
        self.auxfile_types['idealout_lists'] = 'json'


class CliffordRBInput(BenchmarkingInput):

    def __init__(self, pspec, depths, circuits_per_depth, qubit_labels=None, randomizeout=False,
                 citerations=20, compilerargs=[], descriptor='A Clifford RB experiment',
                 verbosity=1):
        #Translated from clifford_rb_experiment
        if qubit_labels is None: qubit_labels = tuple(pspec.qubit_labels)
        circuit_lists = []
        ideal_outs = []
    
        for lnum, l in enumerate(depths):
            if verbosity > 0:
                print('- Sampling {} circuits at CRB length {} ({} of {} depths)'.format(circuits_per_depth, l,
                                                                                         lnum + 1, len(depths)))
                print('  - Number of circuits sampled = ', end='')
            circuits_at_depth = []
            idealouts_at_depth = []
            for j in range(circuits_per_depth):
                c, iout = _rb.sample.clifford_rb_circuit(pspec, l, subsetQs=qubit_labels, randomizeout=randomizeout,
                                                  citerations=citerations, compilerargs=compilerargs)
                circuits_at_depth.append(c)
                idealouts_at_depth.append(iout)
                if verbosity > 0: print(j + 1, end=',')
            circuit_lists.append(circuits_at_depth)
            ideal_outs.append(idealouts_at_depth)
            if verbosity > 0: print('')

        super().__init__(depths, circuit_lists, ideal_outs, qubit_labels)
        self.circuits_per_depth = circuits_per_depth
        self.randomizeout = randomizeout
        self.citerations = citerations
        self.compilerargs = compilerargs
        self.descriptor = descriptor


class DirectRBInput(BenchmarkingInput):

    def __init__(self, pspec, depths, circuits_per_depth, qubit_labels=None, sampler='Qelimination', samplerargs=[],
                 addlocal=False, lsargs=[], randomizeout=False, cliffordtwirl=True, conditionaltwirl=True,
                 citerations=20, compilerargs=[], partitioned=False, descriptor='A DRB experiment',
                 verbosity=1):

        if qubit_labels is None: qubit_labels = tuple(pspec.qubit_labels)
        circuit_lists = []
        ideal_outs = []

        for lnum, l in enumerate(depths):
            if verbosity > 0:
                print('- Sampling {} circuits at DRB length {} ({} of {} depths)'.format(circuits_per_depth, l,
                                                                                         lnum + 1, len(depths)))
                print('  - Number of circuits sampled = ', end='')
            circuits_at_depth = []
            idealouts_at_depth = []
            for j in range(circuits_per_depth):
                c, iout = _rb.sample.direct_rb_circuit(
                    pspec, l, subsetQs=qubit_labels, sampler=sampler, samplerargs=samplerargs,
                    addlocal=addlocal, lsargs=lsargs, randomizeout=randomizeout,
                    cliffordtwirl=cliffordtwirl, conditionaltwirl=conditionaltwirl,
                    citerations=citerations, compilerargs=compilerargs,
                    partitioned=partitioned)
                circuits_at_depth.append(c)
                idealouts_at_depth.append((''.join(map(str, iout)),))
                if verbosity > 0: print(j + 1, end=',')
            circuit_lists.append(circuits_at_depth)
            ideal_outs.append(idealouts_at_depth)
            if verbosity > 0: print('')

        super().__init__(depths, circuit_lists, ideal_outs, qubit_labels)
        self.circuits_per_depth = circuits_per_depth
        self.sampler = sampler
        self.samplerargs = samplerargs
        self.addlocal = addlocal
        self.lsargs = lsargs
        self.randomizeout = randomizeout
        self.cliffordtwirl = cliffordtwirl
        self.conditionaltwirl = conditionaltwirl
        self.citerations = citerations
        self.compilerargs = compilerargs
        self.partitioned = partitioned
        self.descriptor = descriptor


#TODO: maybe need more input types for simultaneous RB and mirrorRB "experiments"


class Benchmark(_proto.Protocol):

    summary_datatypes = ('success_counts', 'total_counts', 'hamming_distance_counts', 'success_probabilities')
    dscmp_datatypes =  ('tvds', 'pvals', 'jsds', 'llrs', 'sstvds')

    def __init__(self, name):
        super().__init__(name)

    def compute_summary_data(self, data):

        def success_counts(dsrow, circ, idealout):
            if dsrow.total == 0: return 0  # shortcut?
            return dsrow[idealout]

        def hamming_distance_counts(dsrow, circ, idealout):
            nQ = len(circ.line_labels)  # number of qubits
            assert(nQ == len(idealout[-1]))
            hamming_distance_counts = _np.zeros(nQ + 1, float)
            if dsrow.total > 0:
                for outcome_lbl, counts in dsrow.counts.items():
                    outbitstring = outcome_lbl[-1]
                    hamming_distance_counts[_rb.analysis.hamming_distance(outbitstring, idealout[-1])] += counts
            return list(hamming_distance_counts)  # why a list?

        def get_summary_values(icirc, circ, dsrow, idealout):
            sc = success_counts(dsrow, circ, idealout)
            tc = dsrow.total
            ret = {'success_counts': sc,
                   'total_counts': tc,
                   'success_probabilities': _np.nan if tc == 0 else sc / tc,
                   'hamming_distance_counts': hamming_distance_counts(dsrow, circ, idealout)}
            return ret

        return self.compute_dict(data, self.summary_datatypes,
                                 get_summary_values, for_passes='all')

    def compute_circuit_data(self, data):
        names = ['success_counts', 'total_counts', 'hamming_distance_counts', 'success_probabilities']

        def get_circuit_values(icirc, circ, dsrow, idealout):
            ret = {'twoQgate_count': circ.twoQgate_count(),
                   'depth': circ.depth(),
                   'target': idealout,
                   'circuit_index': icirc,
                   'width': len(circ.line_labels)}
            ret.update(dsrow.aux)  # note: will only get aux data from *first* pass in multi-pass data
            return ret

        return self.compute_dict(data, names, get_circuit_values, for_passes="first")

    def compute_dscmp_data(self, data, dscomparator):

        def get_dscmp_values(icirc, circ, dsrow, idealout):
            ret = {'tvds':  dscomparator.tvds.get(circ, _np.nan),
                   'pvals': dscomparator.pVals.get(circ, _np.nan),
                   'jsds':  dscomparator.jsds.get(circ, _np.nan),
                   'llrs':  dscomparator.llrs.get(circ, _np.nan)}
            return ret

        return self.compute_dict(data, "dscmpdata", self.dsmp_datatypes, get_dscmp_values, for_passes="none")

    def compute_predicted_probs(self, data, model):

        def get_success_prob(icirc, circ, dsrow, idealout):
            #if set(circ.line_labels) != set(qubits):
            #    trimmedcirc = circ.copy(editable=True)
            #    for q in circ.line_labels:
            #        if q not in qubits:
            #            trimmedcirc.delete_lines(q)
            #else:
            #    trimmedcirc = circ
            return {'success_probabilities': model.probs(circ)[('success',)]}

        return self.compute_dict(data, ('success_probabilities',),
                                 get_success_prob, for_passes="none")

    #def compute_results_qty(self, results, qtyname, component_names, compute_fn, force=False, for_passes="all"):
    def compute_dict(self, data, component_names, compute_fn, for_passes="all"):

        inp = data.input
        ds = data.dataset

        depths = inp.depths
        qty_data = _proto.NamedDict('Datatype', 'category', None,
                                    {comp: _proto.NamedDict('Depth', 'int', 'float', {depth: [] for depth in depths})
                                     for comp in component_names})
        
        #loop over all circuits
        for depth, circuits_at_depth, idealouts_at_depth in zip(depths, inp.circuit_lists, inp.idealout_lists):
            for icirc, (circ, idealout) in enumerate(zip(circuits_at_depth, idealouts_at_depth)):
                dsrow = ds[circ] if (ds is not None) else None  # stripOccurrenceTags=True ??
                # -- this is where Tim thinks there's a bottleneck, as these loops will be called for each
                # member of a simultaneous experiment separately instead of having an inner-more iteration
                # that loops over the "structure", i.e. the simultaneous qubit sectors.
                #TODO: <print percentage>

                for component_name, val in compute_fn(icirc, circ, dsrow, idealout).items():
                    qty_data[component_name][depth].append(val)  # maybe use a pandas dataframe here?

        return qty_data

    def create_depthwidth_dict(self, depths, widths, fillfn, seriestype):
        return _proto.NamedDict(
            'Depth', 'int', seriestype, {depth: _proto.NamedDict(
                'Width', 'int', seriestype, {width: fillfn() for width in widths}) for depth in depths})


class PassStabilityTest(_proto.Protocol):
    pass


class VolumetricBenchmarkGrid(Benchmark):
    #object that can take, e.g. a multiinput or simultaneous input and create a results object with the desired grid of width vs depth numbers

    def __init__(self, depths='all', widths='all', datatype='success_probabilities',
                 paths='all', statistic='mean', aggregate=True, rescaler='auto',
                 dscomparator=None, name=None):

        super().__init__(name)
        self.depths = depths
        self.widths = widths
        self.datatype = datatype
        self.paths = paths if paths == 'all' else sorted(paths)  # need to ensure paths are grouped by common prefix
        self.statistic = statistic
        self.aggregate = aggregate
        self.dscomparator = dscomparator
        self.rescaler = rescaler

    def run(self, data):
        #Note: implement "run" here, since we want to deal with multi-pass and multi-inputs

        #Since we know that VolumetricBenchmark protocol objects Create a single results just fill
        # in data under the result object's 'volumetric_benchmarks' and 'failure_counts'
        # keys, and these are indexed by width and depth (even though each VolumetricBenchmark
        # only contains data for a single width), we can just "merge" the VB results of all
        # the underlying by-depth datas, so long as they're all for different widths.

        #data can be multi-pass and multi-input
        # when multi-input, can be a multi-input of simult-inputs of VBs
        protocol_info = {}  # TODO: ADD INFO
        if data.is_multipass():
            passes = list(data.items())
            results = _proto.ProtocolResults(data, self.name, protocol_info, 'Pass', 'category')
        else:
            passes = [(None, data)]
            results = _proto.ProtocolResults(data, self.name, protocol_info, 'Qty', 'category')

        # Trim down data based on specs/paths given to describe
        # which data should be used for each depth,width pair, e.g.
        # ('SpecA', ('Q0',Q1','Q2')) for width=3, etc.
        for passname, passdata in passes:
            if self.paths == 'all':
                paths = passdata.get_tree_paths()
                trimmed_passdata = passdata
            else:
                paths = self.paths
                trimmed_passdata = passdata.filter_paths(self.paths)

            #Then run resulting data normally, giving a results object
            # with "top level" dicts correpsonding to different paths
            VB = VolumetricBenchmark(self.depths, self.datatype, self.statistic,
                                     self.rescaler, self.dscomparator, name=self.name)
            separate_results = VB.run(trimmed_passdata)

            #Then we merge/flatten the data from different paths into one depth vs width grid
            vb = _proto.NamedDict('Depth', 'int', None)
            fails = _proto.NamedDict('Depth', 'int', None)
            path_for_gridloc = {}
            for path in paths:
                #TODO: need to be able to filter based on widths... - maybe replace .update calls
                # with something more complicated when width != 'all'
                #print("Aggregating path = ", path)  #TODO - show progress something like this later?
                
                #Traverse path to get to root of VB data
                root = separate_results
                for key in path:
                    root = root[key]
                assert(isinstance(root.data.input, ByDepthInput)), \
                    "All paths must lead to by-depth inputs, not %s!" % str(type(root.data.input))

                #Get the list of depths we'll extract from this (`root`) sub-results
                depths = root.data.input.depths if (self.depths == 'all') else \
                    filter(lambda d: d in self.depths, root.data.input.depths)
                width = len(root.data.input.qubit_labels)  # sub-results contains only a single width

                for depth in depths:
                    if depth not in vb:  # and depth not in fails
                        vb[depth] = _proto.NamedDict('Width', 'int', 'float')
                        fails[depth] = _proto.NamedDict('Width', 'int', None)
                        path_for_gridloc[depth] = {}  # just used for meaningful error message

                    if width in path_for_gridloc[depth]:
                        raise ValueError(("Paths %s and %s both give data for depth=%d, width=%d!  Set the `paths`"
                                          " argument of this VolumetricBenchmarkGrid to avoid this.") %
                                         (str(path_for_gridloc[depth][width]), str(path), depth, width))

                    vb[depth][width] = root['volumetric_benchmarks'][depth][width]
                    fails[depth][width] = root['failure_counts'][depth][width]
                    path_for_gridloc[depth][width] = path

            if self.statistic in ('minmin', 'maxmax') and not self.aggregate:
                self._update_vb_minmin_maxmax(vb)   # aggregate now since we won't aggregate over passes

            if passname is None:
                results.qtys['volumetric_benchmarks'] = vb
                results.qtys['failure_counts'] = fails
            else:
                results.qtys[passname] = _proto.NamedDict('Qty', 'category', None)
                results.qtys[passname]['volumetric_benchmarks'] = vb
                results.qtys[passname]['failure_counts'] = fails

        if self.aggregate and len(passes) > 1:  # aggregate pass data into a single set of qty dicts
            aggregated_results = _proto.ProtocolResults(data, self.name, protocol_info, 'Qty', 'category')
            agg_vb = _proto.NamedDict('Depth', 'int', None)
            agg_fails = _proto.NamedDict('Depth', 'int', None)
            passnames = [passname for passname, _ in passes]
            template = results.qtys[passnames[0]]['volumetric_benchmarks']  # to get widths and depths

            #Get function to aggregate the different per-circuit datatype values
            if self.statistic == 'max' or self.statistic == 'maxmax':
                np_fn = _np.nanmax
            elif self.statistic == 'mean':
                np_fn = _np.nanmean
            elif self.statistic == 'min' or self.statistic == 'minmin':
                np_fn = _np.nanmin
            elif self.statistic == 'dist':
                def np_fn(v): return v  # identity
            else: raise ValueError("Invalid statistic '%s'!" % self.statistic)

            def agg_fn(percircuitdata, width, rescale=True):
                """ Aggregates datatype-data for all circuits at same depth """
                rescaled = self.rescale_function(percircuitdata, width) if rescale else percircuitdata
                if _np.isnan(rescaled).all():
                    return _np.nan
                else:
                    return np_fn(rescaled)

            for depth, template_by_width_data in template.items():
                agg_vb[depth] = _proto.NamedDict('Width', 'int', 'float')
                agg_fails[depth] = _proto.NamedDict('Width', 'int', None)

                for width in template_by_width_data.keys():
                    # ppd = "per pass data"
                    vb_ppd = [results.qtys[passname]['volumetric_benchmarks'][depth][width]
                              for passname in passnames]
                    fail_ppd = [results.qtys[passname]['failure_counts'][depth][width]
                                for passname in passnames]

                    successcount = 0
                    failcount = 0
                    for (successcountpass, failcountpass) in fail_ppd:
                        successcount += successcountpass
                        failcount += failcountpass
                    agg_fails[depth][width] = (successcount, failcount)

                    if self.statistic == 'dist':
                        agg_vb[depth][width] = [item for sublist in vb_ppd for item in sublist]
                    else:
                        agg_vb[depth][width] = agg_fn(vb_ppd, width, rescale=False)

            if self.statistic in ('minmin', 'maxmax'):
                self._update_vb_minmin_maxmax(aggregated_results.qtys['volumetric_benchmarks'])
            results = aggregated_results  # replace per-pass results with aggregated results

        return results

    def _update_vb_minmin_maxmax(self, vb):
        for d in vb.keys():
            for w in vb[d].keys():
                for d2 in vb.keys():
                    for w2 in vb[d2].keys():
                        if self.statistic == 'minmin' and d2 <= d and w2 <= w and vb[d2][w2] < vb[d][w]:
                            vb[d][w] = vb[d2][w2]
                        if self.statistic == 'maxmax' and d2 >= d and w2 >= w and vb[d2][w2] > vb[d][w]:
                            vb[d][w] = vb[d2][w2]


class VolumetricBenchmark(Benchmark):

    def __init__(self, depths='all', datatype='success_probabilities',
                 statistic='mean', rescaler='auto', dscomparator=None, name=None):

        assert(statistic in ('max', 'mean', 'min', 'dist', 'maxmax', 'minmin'))

        super().__init__(name)
        self.depths = depths
        #self.widths = widths  # widths='all',
        self.datatype = datatype
        self.statistic = statistic
        self.dscomparator = dscomparator

        if isinstance(rescaler, str):
            if rescaler == 'auto':
                if datatype == 'success_probabilities':
                    def rescale_function(data, width):
                        return list((_np.array(data) - 1 / 2**width) / (1 - 1 / 2**width))
                else:
                    def rescale_function(data, width):
                        return data
            elif rescaler == 'none':
                def rescale_function(data, width):
                    return data
            else:
                raise ValueError("Unknown rescaling option!")
        else:
            rescale_function = rescaler
        self.rescale_function = rescale_function

    def _run(self, data):

        inp = data.input

        #Note: can only take/put things ("results") in data.cache that *only* depend on the input
        # and dataset (i.e. the DataProtocol object).  Protocols must be careful of this in their implementation!
        if self.datatype in self.summary_datatypes:
            if self.datatype not in data.cache:
                summary_data_dict = self.compute_summary_data(data)
                data.cache.update(summary_data_dict)
        elif self.datatype in self.dscmp_datatypes:
            if self.datatype not in data.cache:
                dscmp_data = self.compute_dscmp_data(data, self.dscomparator)
                data.cache.update(dscmp_data)
        else:
            raise ValueError("Invalid datatype: %s" % self.datatype)
        src_data = data.cache[self.datatype]

        #self.compute_circuit_data(results)
        #self.compute_predicted_probs(results, qtyname, model)

        #Get function to aggregate the different per-circuit datatype values
        if self.statistic == 'max' or self.statistic == 'maxmax':
            np_fn = _np.nanmax
        elif self.statistic == 'mean':
            np_fn = _np.nanmean
        elif self.statistic == 'min' or self.statistic == 'minmin':
            np_fn = _np.nanmin
        elif self.statistic == 'dist':
            def np_fn(v): return v  # identity
        else: raise ValueError("Invalid statistic '%s'!" % self.statistic)

        def agg_fn(percircuitdata, width, rescale=True):
            """ Aggregates datatype-data for all circuits at same depth """
            rescaled = self.rescale_function(percircuitdata, width) if rescale else percircuitdata
            if _np.isnan(rescaled).all():
                return _np.nan
            else:
                return np_fn(rescaled)

        def failcnt_fn(percircuitdata):
            """ Returns (nSucceeded, nFailed) for all circuits at same depth """
            nCircuits = len(percircuitdata)
            failcount = int(_np.sum(_np.isnan(percircuitdata)))
            return (nCircuits - failcount, failcount)

        #TODO REMOVE
        #BEFORE SimultaneousInputs: for qubits in inp.get_structure():        
        #width = len(qubits)

        data_per_depth = src_data
        if self.depths == 'all':
            depths = data_per_depth.keys()
        else:
            depths = filter(lambda d: d in data_per_depth, self.depths)
        width = len(inp.qubit_labels)

        vb = self.create_depthwidth_dict(depths, (width,), lambda: None, 'float')
        fails = self.create_depthwidth_dict(depths, (width,), lambda: None, None)

        for depth in depths:
            percircuitdata = data_per_depth[depth]
            fails[depth][width] = failcnt_fn(percircuitdata)
            vb[depth][width] = agg_fn(percircuitdata, width)

        protocol_info = {}  # TODO: ADD INFO
        results = _proto.ProtocolResults(data, self.name, protocol_info, 'Qty', 'category')
        results.qtys['volumetric_benchmarks'] = vb
        results.qtys['failure_counts'] = fails
        results.auxfile_types['volumetric_benchmarks'] = 'json'
        results.auxfile_types['failure_counts'] = 'json'
        return results


class PredictedData(_proto.Protocol):
    #maybe just a function??

    def __init__(self, name):
        super().__init__(name)

    def _run(self, data):

        for i, ((circ, dsrow), auxdict, (pcirc, pdsrow)) in enumerate(iterator):
            if pcirc is not None:
                if not circ == pcirc:
                    print('-{}-'.format(i))
                    pdsrow = predds[circ]
                    _warnings.warn("Predicted DataSet is ordered differently to the main DataSet!"
                                   + "Reverting to potentially slow dictionary hashing!")

        #Similar to above but only on first dataset... see create_summary_data
            

class RB(_proto.Protocol):
    pass

