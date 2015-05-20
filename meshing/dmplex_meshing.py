from meshing import Meshing
from firedrake import op2
from firedrake.petsc import PETSc
import numpy as np
import uuid
import mpi4py

regions = ['Generate', 'Distribute', 'Refine', 'DistributeOverlap']
petsc_events = { 'Distribute': ['Mesh Partition', 'Mesh Migration', 'SFBcastBegin', 'SFReduceBegin'],
                 'Overlap': ['Mesh Partition', 'Mesh Migration', 'SFBcastBegin', 'SFReduceBegin'],
                 'Redistribute': ['Mesh Partition', 'Mesh Migration', 'SFBcastBegin', 'SFReduceBegin']}
global_sum_keys = ['messageLength', 'numMessages']

seed = 8957382

for stage, events in petsc_events.iteritems():
    for event in events:
        regions += "%s::%s" % (stage, event)

class DMPlexMeshing(Meshing):
    benchmark = 'DMPlex_UnitMesh'
    description = 'DMPlex mesh generation benchmark'

    method = 'meshing'
    profileregions = regions

    def meshing(self, dim=2, size=32, refine=0, partitioner="chaco", redistribute=0):
        # HACK-ALERT: Resetting counters via PETSc.Log().destroy()
        # causes errors, so we add a unique ID to the stages.
        unique = str(uuid.uuid1())
        log = PETSc.Log()
        stage_ol = log.Stage(unique+"Overlap")
        stage_dist = log.Stage(unique+"Distribute")
        stage_redist = log.Stage(unique+"Redistribute")
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

        # Set prescribed partitioner
        part = plex.getPartitioner()
        part.setType(partitioner)
        part.setUp()

        if partitioner == "shell":
            # Create a random (bad!) partitioning
            nprocs = self.meta['np']
            ncells = self.num_cells([size])[0]
            if op2.MPI.comm.rank == 0:
                np.random.seed(seed)
                points = np.random.choice(np.arange(ncells, dtype=np.int32),
                                          ncells, replace=False)
                sizes = [ncells / nprocs for _ in range(nprocs)]
            else:
                points = np.zeros(0, dtype=np.int32)
                sizes = np.zeros(nprocs, dtype=np.int32)
            part.setShellPartition(nprocs, sizes, points)

        stage_dist.push()
        with self.timed_region('Distribute'):
            plex.distribute(overlap=0)
        stage_dist.pop()

        with self.timed_region('Refine'):
            plex.setRefinementUniform(True)
            for i in range(refine):
                plex = plex.refine()

        stage_ol.push()
        with self.timed_region('DistributeOverlap'):
            overlap = 0 if redistribute > 0 else 1
            plex.distributeOverlap(overlap=overlap)
        stage_ol.push()

        if redistribute > 0:
            # Switch to parmetis for parallel re-partitioning
            part = plex.getPartitioner()
            part.setType("parmetis")
            part.setUp()

            # Re-distribute plex with new partitioning
            stage_redist.push()
            with self.timed_region('Redistribute'):
                plex.distribute(overlap)
            stage_redist.pop()

        # Extract petsc timings from log object
        for stage, events in petsc_events.iteritems():
            for event in events:
                info = log.Event(event).getPerfInfo(log.Stage(unique+stage))
                for key in info.keys():
                    value = info[key]
                    if key in global_sum_keys:
                        # Reduce global sum data, eg. total message volume
                        value = op2.MPI.comm.allreduce(value, op=mpi4py.MPI.SUM)
                    self.register_timing("%s::%s::%s" % (stage, event, key), value)


if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    # Benchmark
    b = DMPlexMeshing()
    b.main()
