from pybench import Benchmark, parser
from firedrake import op2

sizes = [int(2**x) for x in range(1, 6)]

num_cells = {2: lambda s: [2*x**2 for x in s],
             3: lambda s: [6*x**3 for x in s]}

regions = ['Generate', 'Distribute']

class Meshing(Benchmark):
    warmups = 0
    repeats = 3

    def __init__(self, **kwargs):
        super(Meshing, self).__init__(**kwargs)
        args, _ = self.parser().parse_known_args()
        self.meta['dim'] = args.dim
        self.meta['sizes'] = args.size
        self.series = {'dim' : self.meta['dim'],
                       'np': op2.MPI.comm.size,
                       'variant': args.branch}
        self.params = [('dim', [self.meta['dim']]),
                       ('size', self.meta['sizes']),
                       ('fs', ['scalar'])]

    def parser(self, **kwargs):
        p = super(Meshing, self).parser(**kwargs)
        p.add_argument('--dim', type=int, default=2,
                       help='Dimension of generated DMPlex (default=2)')
        p.add_argument('-m', '--size', type=int, nargs='+',
                       help='Mesh sizes to use')
        p.add_argument('--branch', default='master',
                       help='PETSc branch used')
        return p

if __name__ == '__main__':
    p = parser(description="Plot results for meshing benchmark")
    p.add_argument('--dim', type=int, default=2,
                   help='Dimension of generated DMPlex (default=2)')
    p.add_argument('-m', '--size', type=int, nargs='+',
                   help='mesh sizes to plot')
    p.add_argument('--branch', nargs='+', default=['master'],
                   help='PETSc branches to plot')
    args = p.parse_args()
    variants = args.branch or ['master']
    groups = ['variant'] if len(args.branch) > 1 else []

    b = Meshing(resultsdir=args.resultsdir, plotdir=args.plotdir)
    if args.weak:
        b.combine_series([('dim', [args.dim]), ('np', args.weak), ('variant', variants)],
                         filename='DMPlex_UnitMesh')
        b.plot(xaxis='size', regions=regions, xlabel='mesh size (cells)', groups=groups,
               xvalues=num_cells[b.meta['dim']](b.meta['sizes']),
               figname=b.name, wscale=0.7, format='pdf')

    if args.parallel:
        b.combine_series([('np', args.parallel), ('dim', [args.dim]), ('variant', variants)],
                         filename='DMPlex_UnitMesh')
        b.plot(xaxis='np', regions=regions, xlabel='Number of processors', groups=groups,
               figname=b.name, wscale=0.7, format='pdf')
