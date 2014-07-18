from firedrake_forms import FiredrakeForms
from firedrake import *
from pyop2.coffee.ast_plan import V_OP_UAJ


class FiredrakeFormsCoffee(FiredrakeForms):
    name = 'FiredrakeFormsCoffee'
    series = {}
    params = [('q', [1, 2, 3, 4]),
              ('p', [1, 2, 3, 4]),
              ('form', ['mass', 'elasticity', 'mixed_poisson', 'poisson', 'helmholtz']),
              ('opt', [(0, False, False, False), (0, False, False, True),
                       (1, True, False, False), (0, False, True, False)])
              ]

    def forms(self, q=1, p=1, dim=3, max_nf=3, form='mass', opt=(0, False, False, False)):
        parameters["coffee"]["licm"] = opt[0]
        parameters["coffee"]["ap"] = opt[1]
        parameters["coffee"]["autotune"] = opt[2]
        parameters["coffee"]["ffc-opt"] = opt[3]
        parameters["form_compiler"]['optimize'] = opt[3]
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
           xlabel='COFFEE Optimisations (BASE, FFC-OPT, LICM1, AUTO)',
           xvalues=['0/n/n/n', '0/n/n/y', '1/y/n/n', '0/n/y/n'],
           title='%(form)s (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
    b.plot(xaxis='opt', regions=regions, kinds='bar',
           xlabel='COFFEE Optimisations (BASE, FFC-OPT, LICM1, AUTO)',
           ylabel='Speedup over unoptimised baseline', speedup=((0, False, False, False),),
           xvalues=['0/n/n/n', '0/n/n/y', '1/y/n/n', '0/n/y/n'],
           title='%(form)s (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
    for i, r in enumerate(regions):
        b.plot(xaxis='opt', regions=[r], kinds='bar,barlog',
               xlabel='COFFEE Optimisations (BASE, FFC-OPT, LICM1, AUTO)', groups=['form'],
               figname='FiredrakeFormsCoffee_nf%d' % i,
               xvalues=['0/n/n/n', '0/n/n/y', '1/y/n/n', '0/n/y/n'],
               title=str(i) + ' premultiplying functions (single core, 3D, degree q = %(q)d, premultiplying degree p = %(p)d)')
