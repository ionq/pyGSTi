""" RB-related functions of gates and gatesets """
from __future__ import division, print_function, absolute_import, unicode_literals
#*****************************************************************
#    pyGSTi 0.9:  Copyright 2015 Sandia Corporation
#    This Software is released under the GPL license detailed
#    in the file "license.txt" in the top-level pyGSTi directory
#*****************************************************************

from . import analysis as _analysis
from ... import tools as _tls
from ... import objects as _objs
from ... import construction as _cnst
from ... import algorithms as _algs

import numpy as _np
import warnings as _warnings

def average_gate_infidelity(A ,B, mxBasis="gm"):
    # Todo : docstring
    return 1 - _tls.average_gate_fidelity(A ,B, mxBasis)

def entanglement_infidelity(A, B, mxBasis=None):
    # Todo : docstring
    return 1 - float(_tls.process_fidelity(A, B, mxBasis))

def gateset_infidelity(gs, gs_target, itype = 'EI', 
                       weights=None, mxBasis=None):
    """
    Computes the average-over-gates of the average gate infidelity between 
    gates in `gs` and the gates in `gs_target`. This quantity is sometimes
    called the "average error rate" and Proctor et al Phys. Rev. Lett. 119, 
    130502 (2017) it is called the average gateset infidelity.
    
    Parameters
    ----------
    gs : GateSet
        The gateset to calculate the average infidelity, to `gs_target`, of.

    gs_target : GateSet
        The gateset to calculate the average infidelity, to `gs`, of.
        
    itype : str, optional
        The infidelity type. Either 'EI', corresponding to entanglement
        infidelity, or 'AGI', corresponding to average gate infidelity.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates 
        in `gs` and the values are, possibly unnormalized, probabilities. 
        These probabilities corresponding to the weighting in the average,
        so if the gateset contains gates A and B and weights[A] = 2 and
        weights[B] = 1 then the output is Inf(A)*2/3  + Inf(B)/3 where
        Inf(X) is the infidelity (to the corresponding element in the other
        gateset) of X. If None, a uniform-average is taken, equivalent to
        setting all the weights to 1.

    mxBasis : {"std","gm","pp"} or Basis object, optional
        The basis of the gatesets. If None, the basis is obtained from
        the gateset.
        
    Returns
    -------
    float
        The average infidelity between the two gatesets.
        
    """
    assert(itype == 'AGI' or itype == 'EI'), "The infidelity type must be `AGI` (average gate infidelity) or `EI` (entanglement infidelity)"
    
    if mxBasis is None: mxBasis = gs.basis
        
    sum_of_weights = 0
    I_list = []
    for gate in list(gs_target.gates.keys()):
        if itype == 'AGI':
            I = average_gate_infidelity(gs[gate],gs_target[gate], mxBasis=mxBasis)
        if itype == 'EI':
            I = entanglement_infidelity(gs[gate],gs_target[gate], mxBasis=mxBasis)
        if weights is None:
            w = 1
        else:
            w = weights[gate]
       
        I_list.append(w*I)
        sum_of_weights += w
        
    assert(sum_of_weights > 0), "The sum of the weights should be positive!"
    AI = _np.sum(I_list)/sum_of_weights
        
    return AI

def predicted_RB_number(gs, gs_target, weights=None, d=None, rtype='EI'):
    """
    Predicts the RB number (aka, the RB error rate) from a gateset, 
    using the "L-matrix" theory from Proctor et al Phys. Rev. Lett. 119, 
    130502 (2017). Note that this gives the same predictions as the 
    theory by Wallman Quantum 2, 47 (2018).
    
    This theory is valid for various types of RB, including standard 
    Clifford RB -- i.e., it will accurately predict the per-Clifford 
    error rate reported by standard Clifford RB. It is also valid for
    "direct RB" under broad circumstances.
    
    For this function to be valid the gateset should be trace preserving 
    and completely positive in some representation, but the particular 
    representation of the gateset used is irrelevant, as the predicted RB 
    error rate is a gauge-invariant quantity. The function is likely reliable 
    when complete positivity is slightly violated, although the theory on
    which it is based assumes complete positivity.
    
    Parameters
    ----------
    gs : GateSet
        The gateset to calculate the RB number of. This gateset is the 
        gateset randomly sampled over, so this is not necessarily the 
        set of physical primitives. In Clifford RB this is a set of 
        Clifford gates; in "direct RB" this normally would be the 
        physical primitives.

    gs_target: GateSet
        The target gateset, corresponding to `gs`. This function is not invariant 
        under swapping `gs` and `gs_target`: this GateSet must be the target gateset,
        and should consistent of perfect gates.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates 
        in `gs` and the values are the unnormalized probabilities to apply
        each gate at each stage of the RB protocol. If not None, the values 
        in weights must all be non-negative, and they must not all be zero. 
        Because, when divided by their sum, they must be a valid probability
        distribution. If None, the weighting defaults to an equal weighting 
        on all gates, as this is used in many RB protocols (e.g., Clifford RB).
        But, this weighting is flexible in the "direct RB" protocol.
    
    d : int, optional
        The Hilbert space dimension.  If None, then sqrt(gs.dim) is used.
        
    rtype : str, optional
        The type of RB error rate, either "EI" or "AGI", corresponding to 
        different dimension-dependent rescalings of the RB decay constant
        p obtained from fitting to Pm = A + Bp^m. "EI" corresponds to
        an RB error rate that is associated with entanglement infidelity, which
        is the probability of error for a gate with stochastic errors. This is 
        the RB error rate defined in the "direct RB" protocol. "AGI" corresponds
        to an RB error rate that is associated with average gate infidelity. This 
        is the more standard (but perhaps less well motivated) definition of the 
        RB error rate.

    Returns
    -------
    
    r : float.
        The predicted RB number.
        
    """
    if d is None: d = int(round(_np.sqrt(gs.dim)))
    p = predicted_RB_decay_parameter(gs, gs_target, weights=weights)
    r =  _analysis.p_to_r(p, d=d, rtype=rtype)
    return r

def predicted_RB_decay_parameter(gs, gs_target, weights=None):
    """
    Computes the second largest eigenvalue of the 'L matrix' (see the `L_matrix`
    function). For standard Clifford RB and direct RB, this corresponds to the 
    RB decay parameter p in Pm = A + Bp^m for "reasonably low error" trace 
    preserving and completely positive gates. See also the `predicted_RB_number` 
    function.
    
    Parameters
    ----------
    gs : Gateset
        The gateset to calculate the RB decay parameter of. This gateset is the 
        gateset randomly sampled over, so this is not necessarily the 
        set of physical primitives. In Clifford RB this is a set of 
        Clifford gates; in "direct RB" this normally would be the 
        physical primitives.
  
    gs_target : Gateset
        The target gateset corresponding to gs. This function is not invariant under
        swapping `gs` and `gs_target`: this GateSet must be the target gateset, and
        should consistent of perfect gates.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates 
        in `gs` and the values are the unnormalized probabilities to apply
        each gate at each stage of the RB protocol. If not None, the values 
        in weights must all be non-negative, and they must not all be zero. 
        Because, when divided by their sum, they must be a valid probability
        distribution. If None, the weighting defaults to an equal weighting 
        on all gates, as this is used in many RB protocols (e.g., Clifford RB).
        But, this weighting is flexible in the "direct RB" protocol.

    Returns
    -------
    
    p : float.
        The second largest eigenvalue of L. This is the RB decay parameter
        for various types of RB.
    """
    L = L_matrix(gs, gs_target, weights=weights)
    E = _np.absolute(_np.linalg.eigvals(L))
    E = _np.flipud(_np.sort(E))
    if abs(E[0] - 1) > 10**(-12):
        _warnings.warn("Output may be unreliable because the gateset is not approximately trace-preserving.")

    if E[1].imag > 10**(-10):
        _warnings.warn("Output may be unreliable because the RB decay constant has a significant imaginary component.")
    p = E[1]
    return p

def rb_gauge(gs, gs_target, weights=None, mxBasis=None, eigenvector_weighting=1.0):
    """
    Computes the gauge transformation required so that, when the gateset is transformed
    via this gauge-transformation, the RB number -- as predicted by the function 
    `predicted_RB_number` -- is the average gateset infidelity between the transformed
    `gs` gateset and the target gateset `gs_target`. This transformation is defined
    Proctor et al Phys. Rev. Lett. 119, 130502 (2017), and see also Wallman Quantum 2, 
    47 (2018).
    
    Parameters
    ----------
    gs : Gateset
        The RB gateset. This is not necessarily the set of physical primitives -- it 
        is the gateset randomly sampled over in the RB protocol (e.g., the Cliffords).
  
    gs_target : Gateset
        The target gateset corresponding to gs. This function is not invariant under
        swapping `gs` and `gs_target`: this GateSet must be the target gateset, and
        should consistent of perfect gates.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates 
        in `gs` and the values are the unnormalized probabilities to apply
        each gate at each stage of the RB protocol. If not None, the values 
        in weights must all be non-negative, and they must not all be zero. 
        Because, when divided by their sum, they must be a valid probability
        distribution. If None, the weighting defaults to an equal weighting 
        on all gates, as this is used in many RB protocols (e.g., Clifford RB).
        But, this weighting is flexible in the "direct RB" protocol.
                
    mxBasis : {"std","gm","pp"}, optional
        The basis of the gatesets. If None, the basis is obtained from the gateset.
        
    eigenvector_weighting : float, optional
        Must be non-zero. A weighting on the eigenvector with eigenvalue that
        is the RB decay parameter, in the sum of this eigenvector and the
        eigenvector with eigenvalue of 1 that defines the returned matrix `l_operator`. 
        The value of this factor does not change whether this `l_operator` transforms into
        a gauge in which r = AGsI, but it may impact on other properties of the
        gates in that gauge. It is irrelevant if the gates are unital.
        
    Returns
    -------    
    l_operator: array
        The matrix defining the gauge-transformation.
        
    """                    
    gam, vecs = _np.linalg.eig(L_matrix(gs,gs_target,weights=weights))
    absgam = abs(gam)
    index_max = _np.argmax(absgam)
    gam_max = gam[index_max]

    if abs(gam_max - 1) > 10**(-12):
        _warnings.warn("Output may be unreliable because the gateset is not approximately trace-preserving.")        
      
    absgam[index_max] = 0.0
    index_2ndmax = _np.argmax(absgam)
    decay_constant = gam[index_2ndmax]
    if decay_constant.imag > 10**(-12):
        _warnings.warn("Output may be unreliable because the RB decay constant has a significant imaginary component.")

    vec_l_operator = vecs[:,index_max] + eigenvector_weighting*vecs[:,index_2ndmax]
    
    if mxBasis is None:
        mxBasis = gs.basis.name
    assert(mxBasis=='pp' or mxBasis=='gm' or mxBasis=='std'), "mxBasis must be 'gm', 'pp' or 'std'."
    
    if mxBasis is 'pp' or 'gm':
        assert(_np.amax(vec_l_operator.imag) < 10**(-15)), "If 'gm' or 'pp' basis, RB gauge matrix should be real."
        vec_l_operator = vec_l_operator.real
        
    vec_l_operator[abs(vec_l_operator) < 10**(-15)] = 0.
    l_operator = _tls.unvec(vec_l_operator) 
    
    return l_operator

def transform_to_rb_gauge(gs, gs_target, weights=None, mxBasis=None, eigenvector_weighting=1.0):
    """
    Transforms a GateSet into the "RB gauge" (see the `RB_gauge` function), as 
    introduced in Proctor et al Phys. Rev. Lett. 119, 130502 (2017). This gauge 
    is a function of both the gateset and its target. These may be input in any 
    gauge, for the purposes of obtaining "r = average gateset infidelity" between 
    the output GateSet and gs_target.
    
    Parameters
    ----------
    gs : Gateset
        The RB gateset. This is not necessarily the set of physical primitives -- it 
        is the gateset randomly sampled over in the RB protocol (e.g., the Cliffords).
  
    gs_target : Gateset
        The target gateset corresponding to gs. This function is not invariant under
        swapping `gs` and `gs_target`: this GateSet must be the target gateset, and
        should consistent of perfect gates.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates 
        in `gs` and the values are the unnormalized probabilities to apply
        each gate at each stage of the RB protocol. If not None, the values 
        in weights must all be non-negative, and they must not all be zero. 
        Because, when divided by their sum, they must be a valid probability
        distribution. If None, the weighting defaults to an equal weighting 
        on all gates, as this is used in many RB protocols (e.g., Clifford RB).
        But, this weighting is flexible in the "direct RB" protocol.
                
    mxBasis : {"std","gm","pp"}, optional
        The basis of the gatesets. If None, the basis is obtained from the gateset.
        
    eigenvector_weighting : float, optional
        Must be non-zero. A weighting on the eigenvector with eigenvalue that
        is the RB decay parameter, in the sum of this eigenvector and the
        eigenvector with eigenvalue of 1 that defines the returned matrix `l_operator`. 
        The value of this factor does not change whether this `l_operator` transforms into
        a gauge in which r = AGsI, but it may impact on other properties of the
        gates in that gauge. It is irrelevant if the gates are unital.
        
    Returns
    -------   
    gs_in_RB_gauge: GateSet
        The gateset `gs` transformed into the "RB gauge".
        
    """            
    l = rb_gauge(gs,gs_target,weights=weights,mxBasis=mxBasis,
                 eigenvector_weighting=eigenvector_weighting)
    gs_in_RB_gauge = gs.copy()
    S = _objs.FullGaugeGroupElement( _np.linalg.inv(l) )
    gs_in_RB_gauge.transform( S )         
    return gs_in_RB_gauge

def L_matrix(gs, gs_target, weights=None):
    """
    Constructs a generalization of the 'L-matrix' linear operator on superoperators,
    from Proctor et al Phys. Rev. Lett. 119, 130502 (2017), represented as a 
    matrix via the "stack" operation. This eigenvalues of this matrix 
    describe the decay constant (or constants) in an RB decay curve for an 
    RB protocol whereby random elements of the provided gateset are sampled 
    according to the `weights` probability distribution over the
    gateset. So, this facilitates predictions of Clifford RB and direct RB 
    decay curves.
    
    Parameters
    ----------
    gs : Gateset
        The RB gateset. This is not necessarily the set of physical primitives -- it 
        is the gateset randomly sampled over in the RB protocol (e.g., the Cliffords).
  
    gs_target : Gateset
        The target gateset corresponding to gs. This function is not invariant under
        swapping `gs` and `gs_target`: this GateSet must be the target gateset, and
        should consistent of perfect gates.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates 
        in `gs` and the values are the unnormalized probabilities to apply
        each gate at each stage of the RB protocol. If not None, the values 
        in weights must all be non-negative, and they must not all be zero. 
        Because, when divided by their sum, they must be a valid probability
        distribution. If None, the weighting defaults to an equal weighting 
        on all gates, as this is used in many RB protocols (e.g., Clifford RB).
        But, this weighting is flexible in the "direct RB" protocol.
                
    Returns
    -------
    L : float
        A weighted version of the L operator from Proctor et al Phys. Rev. Lett. 
        119, 130502 (2017), represented as a matrix using the 'stacking' convention.  
        
    """
    if weights is None:
        weights = {}
        for key in list(gs_target.gates.keys()):
            weights[key] = 1.
            
    normalizer = _np.sum(_np.array([weights[key] for key in list(gs_target.gates.keys())]))           
    dim = len(list(gs_target.gates.keys()))
    L_matrix = (1 / normalizer) * _np.sum(weights[key]*_np.kron(gs.gates[key].T,
                 _np.linalg.inv(gs_target.gates[key])) for key in list(gs_target.gates.keys()))
    
    return L_matrix

def R_matrix_predicted_RB_decay_parameter(gs, group, group_to_gateset=None, weights=None):
    """
    Returns the second largest eigenvalue of a generalization of the 'R-matrix' [see the 
    `R_matrix` function] introduced in Proctor et al Phys. Rev. Lett. 119, 130502 (2017).
    This number is a prediction of the RB decay parameter for trace-preserving gates and 
    a variety of forms of RB, including Clifford and direct RB. This function creates a
    matrix which scales super-exponentially in the number of qubits.
    
    Parameters
    ----------
    gs : Gateset
        The gateset to predict the RB decay paramter for. If `subset_sampling` is None 
        or `group_to_gateset` is None, the labels of the gates in `gs` should be the 
        same as the labels of the group elements in `group`. For Clifford RB this 
        would be the clifford gateset, for direct RB it would be the primitive gates.
            
    group : MatrixGroup
        The group that the `gs` gateset contains gates from (`gs` does not
        need to be the full group, and could be a subset of `group`). For
        Clifford RB and direct RB, this would be the Clifford group.
        
    subset_sampling : list, optional
        If not None, a list of gate labels from `gs`, for which random sequences of 
        this subset of gates are implemented in the RB protocol. Even if this is 
        all of the gates of `gs`, this list needs to be specified if `gs` and 
        `group` are either (1) not labelled the same (and so `group_to_gateset` is 
        not None), or (2) `gs` is a subset of group.
        
    group_to_gateset : dict, optional
        If not None, a dictionary that maps labels of group elements to labels
        of `gs`. Only used if `subset_sampling` is not None. If `subset_sampling` is 
        not None and the `gs` and `group` elements have the same labels, this dictionary
        is not required. Otherwise it is necessary.
        
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates in `gs` 
        and the values are the unnormalized probabilities to apply each gate at 
        each stage of the RB protocol. If not None, the values in weights must all
        be positive or zero, and they must not all be zero (because, when divided by 
        their sum, they must be a valid probability distribution). If None, the
        weighting defaults to an equal weighting on all gates, as used in most RB
        protocols.
      
    d : int, optional
        The Hilbert space dimension. If None, then sqrt(gs.dim) is used.
    
    Returns
    -------
    p : float
        The predicted RB decay parameter. Valid for standard Clifford RB or direct RB
        with trace-preserving gates, and in a range of other circumstances.
        
    """ 
    R = R_matrix(gs, group, group_to_gateset=group_to_gateset, weights=weights)
    E = _np.absolute(_np.linalg.eigvals(R))
    E = _np.flipud(_np.sort(E))
    p = E[1]
    return p

def R_matrix(gs, group, group_to_gateset=None, weights=None):
    """
    Constructs a generalization of the 'R-matrix' of Proctor et al Phys. 
    Rev. Lett. 119, 130502 (2017). This matrix described the exact behaviour 
    of the average success probablities of RB sequences. This matrix is 
    super-exponentially large in the number of qubits, but can be 
    constructed for 1-qubit gatesets.
    
    Parameters
    ----------
    gs : Gateset
        The noisy gateset (e.g., the Cliffords) to calculate the R matrix of.
        The correpsonding `target` gateset (not required in this function)
        must be equal to or a subset of (a faithful rep of) the group `group`. 
        If group_to_gateset is None, the labels of the gates in gs should be 
        the same as the labels of the corresponding group elements in `group`. 
        For Clifford RB `gs` should be the clifford gateset; for direct RB 
        this should be the native gateset.
            
    group : MatrixGroup
        The group that the `gs` gateset contains gates from. For Clifford RB
        or direct RB, this would be the Clifford group.
        
    group_to_gateset : dict, optional
        If not None, a dictionary that maps labels of group elements to labels
        of gs. This is required if the labels of the gates in `gs` are different
        from the labels of the corresponding group elements in `group`.
      
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates in gs
        and the values are the unnormalized probabilities to apply each gate at 
        for each layer of the RB protocol. If None, the weighting defaults to an 
        equal weighting on all gates, as used in most RB protocols (e.g., Clifford
        RB).
    
    Returns
    -------
    R : float
        A weighted, a subset-sampling generalization of the 'R-matrix' from Proctor 
        et al Phys. Rev. Lett. 119, 130502 (2017).
        
    """    
    if group_to_gateset is None:
        for key in list(gs.gates.keys()):
            assert(key in group.labels), "Gates labels are not in `group`!"
    else: 
        for key in list(gs.gates.keys()):
            assert(key in group_to_gateset.values()), "Gates labels are not in `group_to_gateset`!"
    
    d = int(round(_np.sqrt(gs.dim)))
    group_dim = len(group)
    R_dim = group_dim * d**2
    R = _np.zeros([R_dim,R_dim],float)
    
    if weights is None:
        weights = {}
        for key in list(gs.gates.keys()):
            weights[key] = 1.
        
    normalizer = _np.sum(_np.array([weights[key] for key in list(gs.gates.keys())]))
    
    for i in range(0,group_dim):
        for j in range(0,group_dim):
            label_itoj = group.labels[group.product([group.get_inv(i),j])]
            if group_to_gateset is not None:
                if label_itoj in group_to_gateset:
                    gslabel = group_to_gateset[label_itoj]
                    R[j*d**2:(j+1)*d**2,i*d**2:(i+1)*d**2] = weights[gslabel]*gs.gates[gslabel]
            else:
                if label_itoj in list(gs.gates.keys()):
                    gslabel = label_itoj                
                    R[j*d**2:(j+1)*d**2,i*d**2:(i+1)*d**2] = weights[gslabel]*gs.gates[gslabel]
            
    R = R/normalizer
    
    return R

def exact_RB_ASPs(gs, group, m_max, m_min=0, m_step=1, success_outcomelabel=('0',),
                  group_to_gateset=None, weights=None, compilation=None, group_twirled=False):
    """
    Calculates the exact RB average success probablilites (ASP), using some
    generalizations of the formula given Proctor et al Phys. Rev. Lett. 119, 
    130502 (2017). This formula does not scale well with group size
    and qubit number, and for the Clifford group it is likely only practical for 
    a single qubit.
    
    Parameters
    ----------
    gs : Gateset
        The noisy gateset (e.g., the Cliffords) to calculate the R matrix of.
        The correpsonding `target` gateset (not required in this function)
        must be equal to or a subset of (a faithful rep of) the group `group`. 
        If group_to_gateset is None, the labels of the gates in gs should be 
        the same as the labels of the corresponding group elements in `group`. 
        For Clifford RB `gs` should be the clifford gateset; for direct RB 
        this should be the native gateset.
            
    group : MatrixGroup
        The group that the `gs` gateset contains gates from. For Clifford RB
        or direct RB, this would be the Clifford group.
      
    m_max : int
        The maximal sequence length of the random gates, not including the
        inversion gate.
        
    m_min : int, optional
        The minimal sequence length. Defaults to the smallest valid value of 0.
        
    m_step : int, optional
        The step size between sequence lengths. Defaults to the smallest valid 
        value of 1.
       
    success_outcomelabel : str or tuple, optional
        The outcome label associated with success.
        
   group_to_gateset : dict, optional
        If not None, a dictionary that maps labels of group elements to labels
        of gs. This is required if the labels of the gates in `gs` are different
        from the labels of the corresponding group elements in `group`.
      
    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates in gs
        and the values are the unnormalized probabilities to apply each gate at 
        for each layer of the RB protocol. If None, the weighting defaults to an 
        equal weighting on all gates, as used in most RB protocols (e.g., Clifford
        RB).

    compilation : dict, optional
        If `gs` is not the full group `group` (with the same labels), then a 
        compilation for the group elements, used to implement the inversion gate 
        (and the initial randomgroup element, if `group_twirled` is True). This 
        is a dictionary with the group labels as keys and a gate sequence of the 
        elements of `gs` as values.
    
    group_twirled : bool, optional
        If True, the random sequence starts with a single uniformly random group
        element before the m random elements of `gs`.
     
    Returns
    -------
    m : float
        Array of sequence length values that the ASPs have been calculated for.
        
    P_m : float
        Array containing ASP values for the specified sequence length values.
        
    """    
    if compilation is None:
        for key in list(gs.gates.keys()):
            assert(key in group.labels), "Gates labels are not in `group`, so `compilation must be specified."
        for label in group.labels:
            assert(label in list(gs.gates.keys())), "Some group elements not in `gs`, so `compilation must be specified."
    
    d = int(round(_np.sqrt(gs.dim)))
    i_max = _np.floor((m_max - m_min ) / m_step).astype('int')
    m = _np.zeros(1+i_max,int)
    P_m = _np.zeros(1+i_max,float)
    group_dim = len(group)
    R = R_matrix(gs, group, group_to_gateset=group_to_gateset, weights=weights)
    success_prepLabel = list(gs.preps.keys())[0] #just take first prep
    success_effectLabel = success_outcomelabel[-1] if isinstance(success_outcomelabel,tuple) else success_outcomelabel
    extended_E = _np.kron(_tls.column_basis_vector(0,group_dim).T,gs.povms['Mdefault'][success_effectLabel].T)
    extended_rho = _np.kron(_tls.column_basis_vector(0,group_dim),gs.preps[success_prepLabel])
    
    if compilation is None:
        extended_E = group_dim*_np.dot(extended_E, R)
        if group_twirled is True:  
            extended_rho = _np.dot(R,extended_rho)
    else:
        full_gateset = _cnst.build_alias_gateset(gs,compilation)
        R_fullgroup = R_matrix(full_gateset,group)
        extended_E = group_dim*_np.dot(extended_E, R_fullgroup)
        if group_twirled is True:        
            extended_rho = _np.dot(R_fullgroup, extended_rho)
               
    Rstep = _np.linalg.matrix_power(R,m_step)
    Riterate =  _np.linalg.matrix_power(R,m_min)
    for i in range (0,1+i_max):
        m[i] = m_min + i*m_step
        P_m[i] = _np.dot(extended_E,_np.dot(Riterate,extended_rho))
        Riterate = _np.dot(Rstep,Riterate)

    return m, P_m


def L_matrix_ASPs(gs, gs_target, m_max, m_min=0, m_step=1, success_outcomelabel=('0',),
                  compilation=None, group_twirled=False, weights=None, gauge_optimize=True, 
                  return_error_bounds=False, norm='diamond'):
    """
    Computes RB average survival probablities, as predicted by the 'L-matrix'
    theory of Proctor et al Phys. Rev. Lett. 119, 130502 (2017). Within the function, 
    the gs is gauge-optimized to gs_target. This is *not* optimized to the gauge specified 
    by Proctor et al, but instead performs the standard pyGSTi gauge-optimization (using the 
    frobenius distance). In most cases, this is likely to be a reasonable proxy for the gauge 
    optimization perscribed by Proctor et al.
    
    Parameters
    ----------
    gs : Gateset
        The noisy gateset.
           
    gs_target : Gateset
        The target gateset.
        
    m_max : int
        The maximal sequence length of the random gates, not including the inversion gate.
        
    m_min : int, optional
        The minimal sequence length. Defaults to the smallest valid value of 0.
        
    m_step : int, optional
        The step size between sequence lengths.
        
    success_outcomelabel : str or tuple, optional
        The outcome label associated with success.

    weights : dict, optional
        If not None, a dictionary of floats, whereby the keys are the gates in gs
        and the values are the unnormalized probabilities to apply each gate at 
        for each layer of the RB protocol. If None, the weighting defaults to an 
        equal weighting on all gates, as used in most RB protocols (e.g., Clifford
        RB).

    compilation : dict, optional
        If `gs` is not the full group, then a compilation for the group elements, 
        used to implement the inversion gate (and the initial random group element, 
        if `group_twirled` is True). This is a dictionary with the group labels as
        keys and a gate sequence of the elements of `gs` as values.
    
    group_twirled : bool, optional
        If True, the random sequence starts with a single uniformly random group
        element before the m random elements of `gs`.
        
    gauge_optimize : bool, optional
        If True a gauge-optimization to the target gateset is implemented before 
        calculating all quantities. If False, no gauge optimization is performed. 
        Whether or not a gauge optimization is performed does not affect the rate of 
        decay but it will generally affect the exact form of the decay. E.g., if a 
        perfect gateset is given to the function -- but in the "wrong" gauge -- no
        decay will be observed in the output P_m, but the P_m can be far from 1 (even
        for perfect SPAM) for all m. The gauge optimization is optional, as it is
        not guaranteed to always improve the accuracy of the reported P_m, although when
        gauge optimization is performed this limits the possible deviations of the 
        reported P_m from the true P_m.
        
    return_error_bounds : bool, optional
        Sets whether or not to return error bounds for how far the true ASPs can deviate
        from the values returned by this function.
        
    norm : str, optional
        The norm used in the error bound calculation. Either 'diamond' for the diamond
        norm (the default) or '1to1' for the Hermitian 1 to 1 norm.
    
    Returns
    -------
    m : float
        Array of sequence length values that the ASPs have been calculated for.
        
    P_m : float
        Array containing predicted ASP values for the specified sequence length values.
    
    if error_bounds is True:
    
        lower_bound: float
            Array containing lower bounds on the possible ASP values

        upper_bound: float
            Array containing upper bounds on the possible ASP values
            
    """    
    d = int(round(_np.sqrt(gs.dim)))
        
    if gauge_optimize:
        gs_go = _algs.gaugeopt_to_target(gs,gs_target)
    else:
        gs_go = gs.copy()
    L = L_matrix(gs_go,gs_target,weights=weights)
    success_prepLabel = list(gs.preps.keys())[0] #just take first prep
    success_effectLabel = success_outcomelabel[-1] if isinstance(success_outcomelabel,tuple) else success_outcomelabel
    identity_vec = _tls.vec(_np.identity(d**2,float))
    
    if compilation is not None:
        gs_group = _cnst.build_alias_gateset(gs_go,compilation)
        gs_target_group = _cnst.build_alias_gateset(gs_target,compilation)
        delta = gate_dependence_of_errormaps(gs_group,gs_target_group,norm=norm)
        emaps = errormaps(gs_group,gs_target_group)
        E_eff = _np.dot(gs_go.povms['Mdefault'][success_effectLabel].T,emaps.gates['Gavg'])
        
        if group_twirled is True:
            L_group = L_matrix(gs_group,gs_target_group)
        
    if compilation is None:
        delta = gate_dependence_of_errormaps(gs_go,gs_target,norm=norm)
        emaps = errormaps(gs_go,gs_target)
        E_eff = _np.dot(gs_go.povms['Mdefault'][success_effectLabel].T,emaps.gates['Gavg'])
    
    i_max = _np.floor((m_max - m_min ) / m_step).astype('int')
    m = _np.zeros(1+i_max,int)
    P_m = _np.zeros(1+i_max,float)
    upper_bound = _np.zeros(1+i_max,float)
    lower_bound = _np.zeros(1+i_max,float)
    
    Lstep = _np.linalg.matrix_power(L,m_step)
    Literate =  _np.linalg.matrix_power(L,m_min)
    for i in range (0,1+i_max):
        m[i] = m_min + i*m_step
        if group_twirled:
            L_m_rdd = _tls.unvec(_np.dot(L_group,_np.dot(Literate,identity_vec)))
        else:
            L_m_rdd = _tls.unvec(_np.dot(Literate,identity_vec))
        P_m[i] = _np.dot(E_eff,_np.dot(L_m_rdd,gs_go.preps[success_prepLabel]))
        Literate = _np.dot(Lstep,Literate)
        upper_bound[i] = P_m[i] + delta/2
        lower_bound[i] = P_m[i] - delta/2
        if upper_bound[i] > 1:
            upper_bound[i]=1.
        if lower_bound[i] < 0:
            lower_bound[i]=0.
    if return_error_bounds:    
        return m, P_m, lower_bound, upper_bound
    else:
        return m, P_m
    
    
def errormaps(gs, gs_target):
    """
    Computes the 'left-multiplied' error maps associated with a noisy gate 
    set, along with the average error map. This is the gate set [E_1,...] 
    such that 
    
    G_i = E_iT_i, 
    
    where T_i is the gate which G_i is a noisy 
    implementation of. There is an additional gate in the set, that has 
    the key 'Gavg' and which is the average of the error maps.
    
    Parameters
    ----------
    gs : GateSet
        The imperfect gateset.
    
    gs_target : GateSet
        The target gateset.
    
    Returns
    -------
    errormaps : GateSet
        The left multplied error gates, along with the average error map,
        with the key 'Gavg'.  
        
    """    
    errormaps_gate_list = []
    errormaps = gs.copy()
    for gate in list(gs_target.gates.keys()):
        errormaps.gates[gate] = _np.dot(gs.gates[gate], 
                               _np.transpose(gs_target.gates[gate]))     
        errormaps_gate_list.append(errormaps.gates[gate])
        
    errormaps.gates['Gavg'] = _np.mean( _np.array([ i for i in errormaps_gate_list]), 
                                        axis=0, dtype=_np.float64)           
    return errormaps

def gate_dependence_of_errormaps(gs, gs_target, norm='diamond', mxBasis=None):
    """
    Computes the "gate-dependence of errors maps" parameter defined by
    
    delta_avg = avg_i|| E_i - avg_i(E_i) ||, 
    
    where E_i are the error maps, and the norm is either the diamond norm 
    or the 1-to-1 norm. This quantity is defined in Magesan et al PRA 85 
    042311 2012.
    
    Parameters
    ----------
    gs : GateSet
        The actual gateset
    
    gs_target : GateSet
        The target gateset.
        
    norm : str, optional
        The norm used in the calculation. Can be either 'diamond' for
        the diamond norm, or '1to1' for the Hermitian 1 to 1 norm.
        
    mxBasis : {"std","gm","pp"}, optional
        The basis of the gatesets. If None, the basis is obtained from
        the gateset.
 
    Returns
    -------
    delta_avg : float
        The value of the parameter defined above.
        
    """
    error_gs = errormaps(gs, gs_target)
    delta = []
    
    if mxBasis is None:
        mxBasis = gs.basis.name
    assert(mxBasis=='pp' or mxBasis=='gm' or mxBasis=='std'), "mxBasis must be 'gm', 'pp' or 'std'."
    
    for gate in list(gs_target.gates.keys()):
        if norm=='diamond':
            delta.append(_tls.diamonddist(error_gs.gates[gate],error_gs.gates['Gavg'],
                                          mxBasis=mxBasis))            
        elif norm=='1to1': 
            gate_dif = error_gs.gates[gate]-error_gs.gates['Gavg']
            delta.append(_tls.norm1to1(gate_dif,n_samples=1000, mxBasis=mxBasis,return_list=False))            
        else:
            raise ValueError("Only diamond or 1to1 norm available.")  
            
    delta_avg = _np.mean(delta)    
    return delta_avg

#def Magesan_theory_predicted_decay(gs, gs_target, mlist, success_outcomelabel=('0',), 
#                                   norm='1to1', order='zeroth', return_all = False):
#    
#    assert(order == 'zeroth' or order == 'first')
# 
#    d = int(round(_np.sqrt(gs.dim)))
#    MTPs = {}
#    MTPs['r'] = gateset_infidelity(gs,gs_target,itype='AGI')    
#    MTPs['p'] = _analysis.r_to_p(MTPs['r'],d,rtype='AGI')
#    MTPs['delta'] = gate_dependence_of_errormaps(gs, gs_target, norm)
#    error_gs = errormaps(gs, gs_target)   
#       
#    R_list = []
#    Q_list = []
#    for gate in list(gs_target.gates.keys()):
#        R_list.append(_np.dot(_np.dot(error_gs.gates[gate],gs_target.gates[gate]),
#              _np.dot(error_gs.gates['Gavg'],_np.transpose(gs_target.gates[gate]))))
#        Q_list.append(_np.dot(gs_target.gates[gate],
#              _np.dot(error_gs.gates[gate],_np.transpose(gs_target.gates[gate]))))
#    
#    error_gs.gates['GR'] = _np.mean(_np.array([ i for i in R_list]),axis=0)
#    error_gs.gates['GQ'] = _np.mean(_np.array([ i for i in Q_list]),axis=0)    
#    error_gs.gates['GQ2'] = _np.dot(error_gs.gates['GQ'],error_gs.gates['Gavg'])   
#    error_gs.preps['rhoc_mixed'] = 1./d*_cnst.basis_build_identity_vec(error_gs.basis)#
#
#    #Assumes standard POVM labels
#    povm = _objs.UnconstrainedPOVM( [('0_cm', gs_target.povms['Mdefault']['0']),
#                                     ('1_cm', gs_target.povms['Mdefault']['1'])] )
#    ave_error_gsl = _cnst.gatestring_list([('rho0','Gavg'),('rho0','GR'),('rho0','Gavg','GQ')])
#    data = _cnst.generate_fake_data(error_gs, ave_error_gsl, nSamples=1, sampleError="none")#

#    pr_L_p = data[('rho0','Gavg')][success_outcomelabel]
#    pr_L_I = data[('rho0','Gavg')][success_outcomelabel_cm]
#    pr_R_p = data[('rho0','GR')][success_outcomelabel]
#    pr_R_I = data[('rho0','GR')][success_outcomelabel_cm]
#    pr_Q_p = data[('rho0','Gavg','GQ')][success_outcomelabel]
#    p = MTPs['p']    
#    B_1 = pr_R_I
#    A_1 = (pr_Q_p/p) - pr_L_p + ((p -1)*pr_L_I/p) + ((pr_R_p - pr_R_I)/p)
#    C_1 = pr_L_p - pr_L_I
#    q = average_gate_infidelity(error_gs.gates['GQ2'],_np.identity(d**2,float))
#    q = _analysis.r_to_p(q,d,rtype='AGI')
#    
#    if order  == 'zeroth':
#        MTPs['A'] = pr_L_I
#        MTPs['B'] = pr_L_p - pr_L_I  
#    if order == 'first':
#        MTPs['A'] = B_1
#        MTPs['B'] = A_1 - C_1*(q - 1)/p**2
#        MTPs['C'] = C_1*(q- p**2)/p**2
#    
#    if order == 'zeroth':
#        Pm = MTPs['A'] +  MTPs['B']*MTPs['p']**_np.array(mlist)
#    if order == 'first':
#        Pm = MTPs['A'] +  (MTPs['B'] + _np.array(mlist)*MTPs['C'])*MTPs['p']**_np.array(mlist)
#    
#   sys_eb = (MTPs['delta'] + 1)**(_np.array(mlist)+1) - 1   
#    if order == 'first':
#        sys_eb = sys_eb - (_np.array(mlist)+1)*MTPs['delta']
#    
#    upper = Pm + sys_eb
#    upper[upper > 1]=1.
#
#    lower = Pm - sys_eb
#    lower[lower < 0]=0.
#
#   return mlist, Pm, upper, lower, MTPs