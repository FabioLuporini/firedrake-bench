from pybench import Benchmark, parser
from firedrake import op2

sizes = [int(2**x) for x in range(1, 6)]

num_cells = {2: lambda s: [2*x**2 for x in s],
             3: lambda s: [6*x**3 for x in s]}

regions = ['Distribute', 'DistributeOverlap']
petsc_events = { 'Distribute': ['Mesh Partition', 'Mesh Migration'],
                 'Overlap' : ['Mesh Partition', 'Mesh Migration']}

for stage, events in petsc_events.iteritems():
    for event in events:
        regions.append( "%s::%s" % (stage, event) )

class Meshing(Benchmark):
    warmups = 0
    repeats = 3

    def __init__(self, **kwargs):
        super(Meshing, self).__init__(**kwargs)
        args, _ = self.parser().parse_known_args()
        self.meta['dim'] = args.dim
        self.meta['sizes'] = args.size
        self.meta['refine'] = args.refine
        self.series = {'dim' : self.meta['dim'],
                       'np': op2.MPI.comm.size,
                       'refine': self.meta['refine'],
                       'variant': args.branch}
        self.params = [('dim', [self.meta['dim']]),
                       ('size', self.meta['sizes']),
                       ('refine', [self.meta['refine']]),
                       ('fs', ['scalar'])]

    def parser(self, **kwargs):
        p = super(Meshing, self).parser(**kwargs)
        p.add_argument('--dim', type=int, default=2,
                       help='Dimension of generated DMPlex (default=2)')
        p.add_argument('-m', '--size', type=int, nargs='+',
                       help='Mesh sizes to use')
        p.add_argument('--branch', default='master',
                       help='PETSc branch used')
        p.add_argument('--refine', type=int, default=0,
                       help='Refine level (regular)')
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
        b.combine_series([('np', args.weak), ('dim', [args.dim]), ('variant', variants)],
                         filename='DMPlex_UnitMesh')
        b.plot(xaxis='size', regions=regions, xlabel='Mesh size (cells)', groups=groups,
               xvalues=num_cells[b.meta['dim']](b.meta['sizes']),
               figname='FixedRanks', wscale=0.7, format='pdf',
               title='DMPlex_UnitMesh: dim=%(dim)d, nprocs=%(np)d' )

    if args.parallel:
        b.combine_series([('np', args.parallel), ('dim', [args.dim]), ('variant', variants)],
                         filename='DMPlex_UnitMesh')
        b.plot(xaxis='np', regions=regions, groups=groups, kinds='plot,loglog',
               xlabel='Number of processors', xticklabels=args.parallel,
               figname='DMPlexDistribute', wscale=0.7, format='pdf',
               title='DMPlex_UnitMesh: dim=%(dim)d, size=2^%(size)d')
