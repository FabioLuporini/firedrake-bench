from firedrake_forms import FiredrakeForms
from firedrake import *

import os

opt_name = ['plain', 'quadrature-O', 'tensor', 'uflacs', 'coffee-base', 'coffee-auto']
speedup_opt_name = [i for i in opt_name if i not in ['plain']]
form = 'hyperelasticity'
form_max_nf = {
    'mass': 3,
    'helmholtz': 3,
    'poisson': 3,
    'elasticity': 3,
    'hyperelasticity': 0
}


class FiredrakeFormsCoffee(FiredrakeForms):
    name = 'FiredrakeFormsCoffee'
    series = {}
    params = [('q', [1, 2, 3, 4]),
              ('p', [1, 2, 3, 4]),
              ('form', [form]),
              ('opt', opt_name)
              ]

    def __init__(self):
        super(FiredrakeFormsCoffee, self).__init__()
        self.ffc_failures = {}

    def forms(self, q=1, p=1, dim=3, max_nf=3, form='mass', opt='plain'):
        if opt in ["plain"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = {
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["quadrature-O"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = True
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = {
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["tensor"]:
            parameters["form_compiler"]["representation"] = "tensor"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = False
            parameters["coffee"] = {
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["uflacs"]:
            parameters["form_compiler"]["representation"] = "uflacs"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = False
            parameters["coffee"] = {
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["coffee-base"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = {
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse'),
                "O1": True
            }
        if opt in ["coffee-auto"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = {
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse'),
                "autotune": True
            }
        super(FiredrakeFormsCoffee, self).forms(q, p, dim, form_max_nf[form], form, opt=opt)

if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    # Benchmark
    b = FiredrakeFormsCoffee()
    # Load the test cases that are known to fail due to FFC failure
    if os.path.exists("ffc_failures.dat"):
        with open('ffc_failures.dat', 'r') as f:
            b.ffc_failures = eval(f.read())
    b.main(load=None)
    # Store the test cases that are known to fail due to FFC failure
    with open('ffc_failures.dat', 'w+') as f:
        f.write(str(b.ffc_failures))

    # Plot
    regions = ['nf %d' % i for i in range(form_max_nf[form]+1)]
    b.plot(xaxis='opt', regions=regions, kinds='bar',
           xlabel='Assembly implementation',
           xticklabels=speedup_opt_name,
           format='pdf',
           ylabel='Speedup over unoptimised code',
           speedup=('plain',),
           title='%(form)s (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
