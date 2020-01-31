""" Specification of  IBM Q Yorktown (aka ibmqx2) """
#***************************************************************************************************
# Copyright 2015, 2019 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License.  You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0 or in the LICENSE file in the root pyGSTi directory.
#***************************************************************************************************

qubits = ['Q' + str(x) for x in range(5)]

twoQgate = 'Gcnot'

edgelist = [('Q1', 'Q0'),
            ('Q0', 'Q1'),
            ('Q0', 'Q2'),
            ('Q2', 'Q0'),
            ('Q2', 'Q1'),
            ('Q1', 'Q2'),
            ('Q3', 'Q2'),
            ('Q2', 'Q3'),
            ('Q3', 'Q4'),
            ('Q4', 'Q3'),
            ('Q2', 'Q4'),
            ('Q4', 'Q2')]

spec_format = 'ibmq_v2019'
