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
        self.meta['np'] = op2.MPI.comm.size
        self.meta['partitioner'] = [args.partitioner]
        self.meta['redistribute'] = [args.redistribute]
        self.series = {'dim' : self.meta['dim'],
                       'np': self.meta['np'],
                       'variant': args.branch}
        self.params = [('size', self.meta['sizes']),
                       ('refine', self.meta['refine']),
                       ('partitioner', self.meta['partitioner']),
                       ('redistribute', self.meta['redistribute'])]

    def parser(self, **kwargs):
        p = super(Meshing, self).parser(**kwargs)
        p.add_argument('--dim', type=int, default=2,
                       help='Dimension of generated DMPlex (default=2)')
        p.add_argument('-m', '--size', type=int, nargs='+',
                       help='Mesh sizes to use')
        p.add_argument('--branch', default='master',
                       help='PETSc branch used')
        p.add_argument('--refine', type=int, default=[0], nargs='+',
                       help='Refine level (regular)')
        p.add_argument('--partitioner', default='chaco',
                       help='Partitioner used for initial distribution')
        p.add_argument('--redistribute', default=0,
                       help='Re-distribute DMPlex again after initial distribution')
        return p

    def num_cells(self, size):
        return num_cells[self.meta['dim']](size)

if __name__ == '__main__':
    p = parser(description="Plot results for meshing benchmark")
    p.add_argument('--dim', type=int, default=2,
                   help='Dimension of generated DMPlex (default=2)')
    p.add_argument('-m', '--size', type=int, nargs='+',
                   help='mesh sizes to plot')
    p.add_argument('--branch', nargs='+', default=['master'],
                   help='PETSc branches to plot')
    p.add_argument('--refine', type=int, default=[0], nargs='+',
                   help='Refine level (regular)')
    p.add_argument('--partitioner', default='chaco',
                   help='Partitioner used for initial distribution')
    args = p.parse_args()
    variants = args.branch or ['master']
    groups = ['variant'] if len(args.branch) > 1 else []

    b = Meshing(resultsdir=args.resultsdir, plotdir=args.plotdir)

    if len(args.refine) > 1:
        # Parallel refinement plots take the first "size" param and
        # plot timings for the equivalent meshes with the specified
        # refinement levels, eg. -m 8 --refine 0 1 2 will plot:
        #     (m=8, r=0), (m=4, r=1) and (m=2, r=2)
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from os import path
        import numpy as np

        b.combine_series([('np', args.parallel), ('dim', [args.dim]),
                          ('variant', variants)], filename='DMPlex_UnitMesh')

        # Create a figure
        figsize = (9, 6)
        fig = plt.figure("DMPlexDistribute_refine.pdf", figsize=figsize, dpi=300)
        ax = fig.add_subplot(111)
        # Precompute colormap
        cmap = mpl.cm.get_cmap("Set1")
        colors = [cmap(i) for i in np.linspace(0, 0.9, len(regions)*len(args.refine))]

        # Manually overlay multiple subplots to insert group-subsets,
        # ie. combinations of size/refine that create the same mesh
        for i, r in enumerate(args.refine):
            regionstyles = {'Distribute': {'linestyle': 'solid',
                                           'color': colors[i]},
                            'DistributeOverlap': {'linestyle': 'dashed',
                                                  'color': colors[i]},
                            'Generate': {'linestyle': 'dashdot',
                                         'color': colors[i]}}
            if r > 0:
                regionstyles['Refine'] = {'linestyle': 'dotted',
                                          'color': colors[i]}

            regions = regionstyles.keys()

            m = args.size[0] / 2**r
            params = {'dim': args.dim, 'variant': args.branch[0],
                      'partitioner': args.partitioner}
            groups = {'refine': [r], 'size': [m]}
            labels = {(r, m): "size %d, refine %d" % (m, r)}
            b.subplot(ax, xaxis='np', kind='loglog', xvals=args.parallel,
                      xlabel='Number of processors', xticklabels=args.parallel,
                      title='DMPlexDistribute with Parallel Refinement',
                      regions=regions, groups=groups, params=params,
                      plotstyle=regionstyles, axis='tight',
                      labels='long', legend={'loc': 'best'})

        fname = "ParallelRefine_loglog_dim%d_variant%s.pdf" % (args.dim, 'master')
        fpath = path.join(args.plotdir, fname)
        fig.savefig(path.join(args.plotdir, fname),
                    orientation='landscape', format="pdf",
                    transparent=True, bbox_inches='tight')

    elif args.parallel:
        b.combine_series([('np', args.parallel), ('dim', [args.dim]), ('variant', variants)],
                         filename='DMPlex_UnitMesh')

        b.plot(xaxis='np', regions=regions, groups=groups, kinds='plot,loglog',
               xlabel='Number of processors', xticklabels=args.parallel,
               colormap = "Set1", axis='tight',
               figname='DMPlexDistribute', wscale=0.7, format='pdf',
               title='DMPlex_UnitMesh: dim=%(dim)d, size=2^%(size)d')
