from forms import Forms
from firedrake import *
from firedrake import __version__ as firedrake_version
from pyop2.profiling import get_timers
from pyop2 import __version__ as pyop2_version

parameters["assembly_cache"]["enabled"] = False

meshes = {2: UnitSquareMesh(31, 31), 3: UnitCubeMesh(9, 9, 9)}


def mass(q, p, dim, mesh, nf=0):
    V = FunctionSpace(mesh, 'CG', q)
    P = FunctionSpace(mesh, 'CG', p)
    u = TrialFunction(V)
    v = TestFunction(V)
    it = dot(v, u)
    f = [Function(P).assign(1.0) for _ in range(nf)]
    return reduce(inner, f + [it])*dx


def helmholtz(q, p, dim, mesh, nf=0):
    V = FunctionSpace(mesh, "CG", q)
    P = FunctionSpace(mesh, "CG", p)
    u = TrialFunction(V)
    v = TestFunction(V)
    f = [Function(P).assign(1.0) for _ in range(nf)]
    it = dot(grad(v), grad(u)) + 1.0*v*u
    return reduce(inner, f + [it])*dx


def elasticity(q, p, dim, mesh, nf=0):
    V = VectorFunctionSpace(mesh, 'CG', q)
    P = FunctionSpace(mesh, 'CG', p)
    u = TrialFunction(V)
    v = TestFunction(V)
    eps = lambda v: grad(v) + transpose(grad(v))
    it = 0.25*inner(eps(v), eps(u))
    f = [Function(P).assign(1.0) for _ in range(nf)]
    return reduce(inner, f + [it])*dx


def hyperelasticity(q, p, dim, mesh, nf=0):
    V = VectorFunctionSpace(mesh, 'CG', q)
    P = VectorFunctionSpace(mesh, 'CG', p)
    v = TestFunction(V)
    du = TrialFunction(V)  # Incremental displacement
    u = Function(V)        # Displacement from previous iteration
    B = Function(V)        # Body force per unit mass
    # Kinematics
    I = Identity(v.cell().topological_dimension())
    F = I + grad(u)        # Deformation gradient
    C = F.T*F              # Right Cauchy-Green tensor
    E = (C - I)/2          # Euler-Lagrange strain tensor
    E = variable(E)
    # Material constants
    mu = Constant(1.0)     # Lame's constants
    lmbda = Constant(0.001)
    # Strain energy function (material model)
    psi = lmbda/2*(tr(E)**2) + mu*tr(E*E)
    S = diff(psi, E)       # Second Piola-Kirchhoff stress tensor
    PK = F*S               # First Piola-Kirchoff stress tensor
    # Variational problem
    it = inner(PK, grad(v)) - inner(B, v)
    f = [Function(P).assign(1.0) for _ in range(nf)]
    return derivative(reduce(inner, map(div, f) + [it])*dx, u, du)


def poisson(q, p, dim, mesh, nf=0):
    V = VectorFunctionSpace(mesh, 'CG', q)
    P = VectorFunctionSpace(mesh, 'CG', p)
    u = TrialFunction(V)
    v = TestFunction(V)
    it = inner(grad(v), grad(u))
    f = [div(Function(P).assign(1.0)) for _ in range(nf)]
    return reduce(inner, f + [it])*dx


def mixed_poisson(q, p, dim, mesh, nf=0):
    BDM = FunctionSpace(mesh, "BDM", q)
    DG = FunctionSpace(mesh, "DG", q - 1)
    P = FunctionSpace(mesh, 'CG', p)
    W = BDM * DG
    sigma, u = TrialFunctions(W)
    tau, v = TestFunctions(W)
    it = dot(sigma, tau) + div(tau)*u + div(sigma)*v
    f = [Function(P).assign(1.0) for _ in range(nf)]
    return reduce(inner, f + [it])*dx


class FiredrakeForms(Forms):
    series = {'variant': 'Firedrake'}
    meta = {'coffee': parameters["coffee"],
            'firedrake': firedrake_version,
            'pyop2': pyop2_version}

    def forms(self, q=1, p=1, dim=3, max_nf=3, form='mass', dump_kernel=False):
        mesh = meshes[dim]
        A = assemble(eval(form)(q, p, dim, mesh))

        for nf in range(max_nf + 1):
            f = eval(form)(q, p, dim, mesh, nf)
            with self.timed_region('nf %d' % nf):
                assemble(f, tensor=A)
                A.M
            if dump_kernel:
                for i, k in enumerate(f._kernels):
                    with open('kernels/f_%s%d_q%d_p%d_dim%d_nf%d.c' % (form, i, q, p, dim, nf), 'w') as fil:
                        fil.write(k[-1].code)
        t = get_timers(reset=True)
        task = 'Assemble cells'
        self.register_timing(task, t[task].total)

if __name__ == '__main__':
    op2.init(log_level='WARNING')
    from ffc.log import set_level
    set_level('ERROR')

    FiredrakeForms().main()
