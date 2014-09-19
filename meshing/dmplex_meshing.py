from pybench import Benchmark
from firedrake import *
from firedrake.petsc import PETSc

make_sizes = lambda nmin, nmax: [int(2**x) for x in range(nmin, nmax)]

make_mesh = {2: lambda x: UnitSquareMesh(x, x),
             3: lambda x: UnitCubeMesh(x, x, x)}

num_cells = {2: lambda s: [2*x**2 for x in s],
             3: lambda s: [6*x**3 for x in s]}

regions = ['Generate', 'Distribute']

class DMPlexMeshing(Benchmark):
    benchmark = 'DMPlex_UnitMesh'
    description = 'DMPlex mesh generation benchmark'
    method = 'meshing'
    warmups = 0
    repeats = 3
    profileregions = regions

    def __init__(self):
        super(DMPlexMeshing, self).__init__()
        args, _ = self.parser().parse_known_args()
        self.meta['dim'] = args.dim
        self.meta['sizes'] = make_sizes(args.minsize, args.maxsize)
        self.series = {'dim' : self.meta['dim'], 'np': op2.MPI.comm.size}
        self.params = [('dim', [self.meta['dim']]),
                       ('size', self.meta['sizes']),
                       ('fs', ['scalar'])]

    def parser(self, **kwargs):
        p = super(DMPlexMeshing, self).parser(**kwargs)
        p.add_argument('--dim', type=int, default=2,
                       help='Dimension of generated DMPlex (default=2)')
        p.add_argument('--minsize', '-nmin', type=int, default=1,
                       help='Minimum n, where cube size is (2^n, 2^n)')
        p.add_argument('--maxsize', '-nmax', type=int, default=6,
                       help='Maximum n, where cube size is (2^n, 2^n)')
        return p

    def meshing(self, size=32, degree=1, dim=2, fs='scalar'):
        with self.timed_region('Generate'):
            boundary = PETSc.DMPlex().create(op2.MPI.comm)
            boundary.setDimension(dim-1)
            if dim == 2:
                boundary.createSquareBoundary([0., 0.], [1., 1.], [size, size])
            elif dim == 3:
                boundary.createCubeBoundary([0., 0., 0.], [1., 1., 1.], [size, size, size])
            plex = PETSc.DMPlex().generate(boundary)

        with self.timed_region('Distribute'):
            parallel_sf = plex.distribute(overlap=1)


if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    # Benchmark
    b = DMPlexMeshing()
    b.main()

    b.plot(xaxis='size', regions=regions, xlabel='mesh size (cells)',
           xvalues=num_cells[b.meta['dim']](b.meta['sizes']),
           figname=b.name, wscale=0.7, format='pdf')
