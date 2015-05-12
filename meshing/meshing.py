from pybench import Benchmark, parser
from firedrake import op2

sizes = [int(2**x) for x in range(1, 6)]

num_cells = {2: lambda s: [2*x**2 for x in s],
             3: lambda s: [6*x**3 for x in s]}

def generate_regions(event_dict, regions=[]):
    for stage, events in event_dict.iteritems():
        for event in events:
            if isinstance(event, tuple):
                regions.append( "%s::%s::%s" % (stage, event[0], event[1]) )
            else:
                regions.append( "%s::%s" % (stage, event) )
    return regions


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
        self.series = {'np': self.meta['np'],
                       'variant': args.branch}
        self.params = [('dim', [self.meta['dim']]),
                       ('size', self.meta['sizes']),
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
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from os import path
    import numpy as np

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
    p.add_argument('--redistribute', default=0,
                   help='Redisitrbute mesh for load balance')
    p.add_argument('--notitle', action='store_true', default=False,
                   help='Remove titles from generated plots')
    p.add_argument('--messages', action='store_true', default=False,
                   help='Plot message sizes rather than runtime')
    args = p.parse_args()
    variants = args.branch or ['master']
    groups = ['variant'] if len(args.branch) > 1 else []
    figsize = (6, 5)

    b = Meshing(resultsdir=args.resultsdir, plotdir=args.plotdir)

    if len(args.refine) > 1:
        # Parallel refinement plots take the first "size" param and
        # plot timings for the equivalent meshes with the specified
        # refinement levels, eg. -m 8 --refine 0 1 2 will plot:
        #     (m=8, r=0), (m=4, r=1) and (m=2, r=2)
        b.combine_series([('np', args.parallel), ('dim', [args.dim]),
                          ('variant', variants)], filename='DMPlex_UnitMesh')

        # Create a figure
        figsize = (9, 5)
        fig = plt.figure("DMPlexDistribute_refine.pdf", figsize=figsize, dpi=300)
        ax = fig.add_subplot(111)
        # Precompute colormap
        cmap = mpl.cm.get_cmap("Set1")
        colors = [cmap(i) for i in np.linspace(0, 0.9, 4*len(args.refine))]

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
                      'partitioner': args.partitioner, 'redistribute': 0}
            groups = {'refine': [r], 'size': [m]}
            labels = {(r, m): "size %d, refine %d" % (m, r)}

            regionlabels = {'DistributeOverlap': "Overlap"}
            title = '' if args.notitle else 'Mesh Distribution with Parallel Refinement'
            b.subplot(ax, xaxis='np', kind='loglog', xvals=args.parallel,
                      xlabel='Number of processors', xticklabels=args.parallel,
                      regions=regions, groups=groups, params=params,
                      regionlabels=regionlabels,
                      plotstyle=regionstyles, axis='tight',
                      labels='long', title=title, legend={'loc': 'best', 'ncol': 2})

        fname = "Refine_time_loglog_dim%d_partitioner%s_variant%s.pdf" % (args.dim, args.partitioner, 'master')
        fpath = path.join(args.plotdir, fname)
        fig.savefig(path.join(args.plotdir, fname),
                    orientation='landscape', format="pdf",
                    transparent=True, bbox_inches='tight')

    elif args.redistribute > 0:
        # Precompute colormap
        cmap = mpl.cm.get_cmap("Set1")
        colors = [cmap(i) for i in np.linspace(0, 0.9, 6)]

        if args.messages:
            # Plot average message size rather than runtime
            figname = 'Redistribute_msgs'
            ylabel = "Avg. message size [bytes]"
            colors = colors[1:]
            regions = []
            petsc_events = { 'Redistribute': [('Mesh Partition', 'messageLength'),
                                            ('Mesh Migration', 'messageLength')]}
            regionlabels = {'Redistribute::Mesh Partition::messageLength': "Redistribute: Partition",
                            'Redistribute::Mesh Migration::messageLength': "Redistribute: Migration"}
        else:
            figname = 'Redistribute_time'
            ylabel = 'Time [sec]'
            regions = ['Redistribute']
            petsc_events = { 'Redistribute': [('Mesh Partition', 'time'),
                                              ('Mesh Migration', 'time')]}
            regionlabels = {'Redistribute::Mesh Partition::time': "Redistribute: Partition",
                            'Redistribute::Mesh Migration::time': "Redistribute: Migration"}

        regions = generate_regions(petsc_events, regions)
        b.combine_series([('np', args.parallel), ('dim', [args.dim]),
                          ('variant', variants)], filename='DMPlex_UnitMesh')

        title = '' if args.notitle else 'Mesh distribution (all-to-all)'
        b.plot(xaxis='np', regions=regions, groups=groups, kinds='loglog',
               xlabel='Number of processors', xticklabels=args.parallel,
               ylabel=ylabel, regionlabels=regionlabels,
               colors=colors, axis='tight', figsize=figsize,
               figname=figname, title=title, format='pdf')

    elif args.parallel:
        # Precompute colormap
        cmap = mpl.cm.get_cmap("Set1")
        colors = [cmap(i) for i in np.linspace(0, 0.9, 6)]

        if args.messages:
            # Plot average message size rather than runtime
            figname = 'Distribute_msgs'
            ylabel = "Avg. message size [bytes]"
            colors = colors[2:]
            regions = []
            petsc_events = { 'Distribute': [('Mesh Partition', 'messageLength'),
                                            ('Mesh Migration', 'messageLength')],
                             'Overlap' : [('Mesh Partition', 'messageLength'),
                                          ('Mesh Migration', 'messageLength')]}
            regionlabels = {'Distribute::Mesh Partition::messageLength': "Distribute: Partition",
                            'Distribute::Mesh Migration::messageLength': "Distribute: Migration",
                            'Overlap::Mesh Partition::messageLength': "Overlap: Partition",
                            'Overlap::Mesh Migration::messageLength': "Overlap: Migration"}
        else:
            figname = 'Distribute_time'
            ylabel = 'Time [sec]'
            regions = ['Distribute', 'DistributeOverlap']
            petsc_events = { 'Distribute': [('Mesh Partition', 'time'),
                                            ('Mesh Migration', 'time')],
                             'Overlap' : [('Mesh Partition', 'time'),
                                          ('Mesh Migration', 'time')]}
            regionlabels = {'DistributeOverlap': "Overlap",
                            'Distribute::Mesh Partition::time': "Distribute: Partition",
                            'Distribute::Mesh Migration::time': "Distribute: Migration",
                            'Overlap::Mesh Partition::time': "Overlap: Partition",
                            'Overlap::Mesh Migration::time': "Overlap: Migration"}

        b.combine_series([('np', args.parallel), ('dim', [args.dim]),
                         ('variant', variants)], filename='DMPlex_UnitMesh')
        regions = generate_regions(petsc_events, regions)
        title = '' if args.notitle else 'Mesh distribution (one-to-all)'
        b.plot(xaxis='np', regions=regions, groups=groups, kinds='loglog',
               xlabel='Number of processors', xticklabels=args.parallel,
               ylabel=ylabel, regionlabels=regionlabels,
               colors=colors, axis='tight', figsize=figsize,
               figname=figname, title=title, format='pdf')
