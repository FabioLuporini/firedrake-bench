from ordering import Ordering
from firedrake import *
from pyop2.profiling import get_timers

from firedrake_common import FiredrakeBenchmark

parameters["coffee"]["licm"] = True
parameters["coffee"]["ap"] = True

class FiredrakeOrdering(FiredrakeBenchmark, Ordering):

    def reorder(self, size=64, dim=2, degree=1, reorder=True, refine=0, scale=1.0):
        self.series['size'] = size
        self.series['scale'] = scale
        self.series['dim'] = dim
        self.series['reorder'] = reorder
        self.params = [('degree', degree)]
        self.meta['cells'] = (2 if dim == 2 else 6)*size**dim
        self.meta['dofs'] = (size+1)**dim

        mesh = Mesh("meshes/box_cylinder_%dd_%s.msh" % (dim, scale), reorder=reorder)

        with self.timed_region('setup'):
            V = FunctionSpace(mesh, "CG", degree)
            u = Function(V)

        with self.timed_region('cell_integral'):
            assemble(u*dx)

        with self.timed_region('facet_integral'):
            assemble(u('+')*dS)

        for task, timer in get_timers(reset=True).items():
            self.register_timing(task, timer.total)


if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    FiredrakeOrdering().main()
