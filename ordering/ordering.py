from pybench import Benchmark, parser

# Create a series of meshes that roughly double in number of DOFs
sizes = lambda dim: [int((1e4*2**x)**(1./dim)) + 1 for x in range(4)]
cells = {2: [2*x**2 for x in sizes(2)],
         3: [6*x**3 for x in sizes(3)]}
regions = ['cell_integral', 'facet_integral']


class Ordering(Benchmark):

    params = [('degree', range(1, 4))]
    method = 'reorder'
    benchmark = 'Reorder'
    plotstyle = {'total': {'marker': '*'},
                 'mesh': {'marker': '+'},
                 'setup': {'marker': 'x'},
                 'cell_integral': {'marker': '>'},
                 'facet_integral': {'marker': '<'}}
    profilegraph = {'format': 'pdf',
                    'node_threshold': 2.0}
    profileregions = regions

if __name__ == '__main__':
    p = parser(description="Plot results for advection-diffusion benchmark")
    p.add_argument('--dim', type=int, default=3,
                   help='problem dimension to plot')
    p.add_argument('-d', '--degree', type=int, nargs='+',
                   help='polynomial degrees to plot')
    p.add_argument('-m', '--size', type=int, nargs='+',
                   help='mesh sizes to plot')
    p.add_argument('-v', '--variant', nargs='+', help='variants to plot')
    args = p.parse_args()
    variants = args.variant or ['Firedrake']
    groups = ['reorder']
    degrees = args.degree or [1, 2, 3]
    dim = args.dim
    if args.parallel:
        b = Ordering(benchmark='ReorderingParallel', 
                    resultsdir=args.resultsdir, plotdir=args.plotdir)
        b.combine_series([('np', args.parallel), ('variant', variants),
                          ('dim', [dim]), ('reorder', [True, False]), 
                          ('size', args.size or sizes(dim))], filename='Reorder')
        b.plot(xaxis='np', regions=regions, xlabel='Number of processors',
               kinds='plot,loglog', groups=groups, format='pdf',
               title='RCM Reordering (strong scaling, %(dim)dD, degree %(degree)d, mesh size %(size)s**2)')

