from firedrake_forms import FiredrakeForms
from firedrake import *

import os

opt_name = ['quadrature-O', 'tensor', 'coffee-base', 'coffee-auto']

class FiredrakeFormsCoffee(FiredrakeForms):
    name = 'FiredrakeFormsCoffee'
    series = {}
    params = [('q', [1, 2, 3, 4]),
              ('p', [1, 2, 3, 4]),
              ('form', ['mass', 'elasticity', 'mixed_poisson', 'poisson', 'helmholtz']),
              ('opt', opt_name)
              ]

    def forms(self, q=1, p=1, dim=3, max_nf=3, form='mass', opt='plain'):
        if opt in ["plain"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = { \
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["quadrature-O"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = True
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = { \
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["tensor"]:
            parameters["form_compiler"]["representation"] = "tensor"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = False
            parameters["coffee"] = { \
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse')
            }
        if opt in ["coffee-base"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = { \
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse'),
                "licm": 1,
                "ap": True
            }
        if opt in ["coffee-auto"]:
            parameters["form_compiler"]["representation"] = "quadrature"
            parameters["form_compiler"]['optimize'] = False
            parameters["form_compiler"]["pyop2-ir"] = True
            parameters["coffee"] = { \
                "compiler": os.environ.get('PYOP2_BACKEND_COMPILER', 'gnu'),
                "simd_isa": os.environ.get('PYOP2_SIMD_ISA', 'sse'),
                "autotune": True
            }
        super(FiredrakeFormsCoffee, self).forms(q, p, dim, max_nf, form)

if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    # Benchmark
    b = FiredrakeFormsCoffee()
    b.main(load=None)

    # Plot
    regions = ['nf %d' % i for i in range(4)]
    b.plot(xaxis='opt', regions=regions, kinds='bar,barlog',
           xlabel='Local Assembly Execution Mode',
           xticklabels=opt_name,
           title='%(form)s (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
    b.plot(xaxis='opt', regions=regions, kinds='bar',
           xlabel='Local Assembly Execution Mode',
           xticklabels=opt_name,
           ylabel='Speedup over unoptimised code', speedup=('plain',),
           title='%(form)s (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
    for i, r in enumerate(regions):
        b.plot(xaxis='opt', regions=[r], kinds='bar,barlog',
               xlabel='Local Assembly Execution Mode', groups=['form'],
               xticklabels=opt_name,
               figname='FiredrakeFormsCoffee_nf%d' % i,
               title=str(i) + ' premultiplying functions (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
