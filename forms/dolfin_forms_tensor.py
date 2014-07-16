from forms import Forms
from dolfin import *

# Form compiler options
parameters["form_compiler"]["optimize"] = False
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "tensor"
parameters["form_compiler"]["cpp_optimize_flags"] = '-xAVX'


def mass(q, p, dim, mesh, nf=0):
    V = FunctionSpace(mesh, 'CG', q)
    P = FunctionSpace(mesh, 'CG', p)
    u = TrialFunction(V)
    v = TestFunction(V)
    it = dot(v, u)
    f = [Function(P) for _ in range(nf)]
    for f_ in f:
        f_.interpolate(Expression('1.0'))
    return reduce(inner, f + [it])*dx


def helmholtz(q, p, dim, mesh, nf=0):
    V = FunctionSpace(mesh, "CG", q)
    P = FunctionSpace(mesh, "CG", p)
    u = TrialFunction(V)
    v = TestFunction(V)
    it = dot(grad(v), grad(u)) + 1.0*v*u
    f = [Function(P) for _ in range(nf)]
    for f_ in f:
        f_.interpolate(Expression('1.0'))
    return reduce(inner, f + [it])*dx


def elasticity(q, p, dim, mesh, nf=0):
    V = VectorFunctionSpace(mesh, 'CG', q)
    P = FunctionSpace(mesh, 'CG', p)
    u = TrialFunction(V)
    v = TestFunction(V)
    eps = lambda v: grad(v) + transpose(grad(v))
    it = 0.25*inner(eps(v), eps(u))
    f = [Function(P) for _ in range(nf)]
    for f_ in f:
        f_.interpolate(Expression('1.0'))
    return reduce(inner, f + [it])*dx


def hyperelasticity(q, p, dim, mesh, nf=0):
    V = VectorFunctionSpace(mesh, 'CG', q)
    P = VectorFunctionSpace(mesh, 'CG', p)
    v = TestFunction(V)
    du = TrialFunction(V)  # Incremental displacement
    u = Function(V)        # Displacement from previous iteration
    B = Function(V)        # Body force per unit mass
    # Kinematics
    d = len(u)
    I = Identity(d)
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
    f = [Function(P) for _ in range(nf)]
    for f_ in f:
        f_.interpolate(Expression(('1.0',)*dim))
    return derivative(reduce(inner, map(div, f) + [it])*dx, u, du)


def poisson(q, p, dim, mesh, nf=0):
    V = VectorFunctionSpace(mesh, 'CG', q)
    P = VectorFunctionSpace(mesh, 'CG', p)
    u = TrialFunction(V)
    v = TestFunction(V)
    it = inner(grad(v), grad(u))
    f = [Function(P) for _ in range(nf)]
    for f_ in f:
        f_.interpolate(Expression(('1.0',)*dim))
    return reduce(inner, map(div, f) + [it])*dx


def mixed_poisson(q, p, dim, mesh, nf=0):
    BDM = FunctionSpace(mesh, "BDM", q)
    DG = FunctionSpace(mesh, "DG", q - 1)
    P = FunctionSpace(mesh, 'CG', p)
    W = BDM * DG
    sigma, u = TrialFunctions(W)
    tau, v = TestFunctions(W)
    it = dot(sigma, tau) + div(tau)*u + div(sigma)*v
    f = [Function(P) for _ in range(nf)]
    for f_ in f:
        f_.interpolate(Expression('1.0'))
    return reduce(inner, f + [it])*dx


class DolfinForms(Forms):
    series = {'variant': 'DOLFIN_TENSOR'}
    params = [('q', [1, 2, 3, 4]),
              ('p', [1, 2]),
              ('form', ['helmholtz'])]

    def forms(self, q=1, p=1, dim=3, form='mass'):
        mesh = UnitSquareMesh(31, 31) if dim == 2 else UnitCubeMesh(9, 9, 9)
        A = assemble(eval(form)(q, p, dim, mesh))

        # Set, for each problem, the maximum number of pre-multiplying
        # functions, to prevent FFC from crashing
        if form == 'mass':
            _nf = 3
        elif form == 'helmholtz':
            _nf = 4
        elif form == 'poisson':
            _nf = 3
        elif form == 'elasticity':
            _nf = 3
        elif form == 'mixed_poisson':
            _nf = 3

        for nf in range(_nf):
            with self.timed_region('nf %d' % nf):
                print "    Benchmark nf = ", nf
                assemble(eval(form)(q, p, dim, mesh, nf), tensor=A)
        t = timings(True)
        task = 'Assemble cells'
        self.register_timing(task, float(t.get(task, 'Total time')))

if __name__ == '__main__':
    set_log_active(False)
    from ffc.log import set_level
    set_level('ERROR')

    DolfinForms().main()
