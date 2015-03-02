from meshing import Meshing
from firedrake import op2
from firedrake.petsc import PETSc

regions = ['Generate', 'Distribute', 'Refine', 'DistributeOverlap']
petsc_events = { 'Distribute': ['Mesh Partition', 'Mesh Migration'],
                 'Overlap' : ['Mesh Partition', 'Mesh Migration']}

for stage, events in petsc_events.iteritems():
    for event in events:
        regions += "%s::%s" % (stage, event)

class DMPlexMeshing(Meshing):
    benchmark = 'DMPlex_UnitMesh'
    description = 'DMPlex mesh generation benchmark'

    method = 'meshing'
    profileregions = regions

    def meshing(self, size=32, degree=1, dim=2, fs='scalar', refine=0):
        log = PETSc.Log()
        stage_dist = log.Stage("Distribute")
        stage_ol = log.Stage("Overlap")
        log.begin()

        with self.timed_region('Generate'):
            boundary = PETSc.DMPlex().create(op2.MPI.comm)
            boundary.setDimension(dim-1)
            if dim == 2:
                boundary.createSquareBoundary([0., 0.], [1., 1.], [size, size])
                boundary.setTriangleOptions("pqezQYSl")
            elif dim == 3:
                boundary.createCubeBoundary([0., 0., 0.], [1., 1., 1.], [size, size, size])
            plex = PETSc.DMPlex().generate(boundary)

        with self.timed_region('Distribute'):
            stage_dist.push()
            plex.distribute(overlap=0)
            stage_dist.pop()

        with self.timed_region('Refine'):
            plex.setRefinementUniform(True)
            for i in range(refine):
                plex = plex.refine()

        with self.timed_region('DistributeOverlap'):
            stage_ol.push()
            plex.distributeOverlap(overlap=1)
            stage_ol.push()

        for stage, events in petsc_events.iteritems():
            for event in events:
                info = log.Event(event).getPerfInfo(log.Stage(stage))
                self.register_timing("%s::%s" % (stage, event), info['time'])


if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    # Benchmark
    b = DMPlexMeshing()
    b.main()
