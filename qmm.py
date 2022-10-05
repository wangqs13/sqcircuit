import numpy as np
import time
import queue

complex64_type = "complex128"

class Orthogonal:

    def __init__(self, equal):
        self.equal = equal
        self.basis = []

    def add(self, a):
        for b in self.basis:
            a -= np.vdot(b, a) * b
        d = np.sqrt(np.real(np.vdot(a, a)))
        a /= d
        assert(equal(np.sqrt(np.real(np.vdot(a, a))), 1))
        self.basis.append(a)

    def isIndependent(self, a):
        res = 0
        for b in self.basis:
            tmp = np.vdot(b, a)
            res += np.real(tmp) * np.real(tmp) + np.imag(tmp) * np.imag(tmp)
        return not equal(res, np.vdot(a, a).real)

class QuantumMealyMachineWithQuantumInput:

    def __init__(self, dim_in, dim_s, unitary, measure):
        self.dim_in = dim_in
        self.dim_s = dim_s
        self.unitary = np.matrix(unitary, dtype = complex64_type, copy = True)
        self.gamma = len(measure)
        self.measure = [np.matrix(x, dtype = complex64_type, copy = True) for x in measure]
        
        dim_x, dim_y = unitary.shape
        assert(dim_x == dim_y and dim_x == dim_in * dim_s)
        for i in range(self.gamma):
            dim_x, dim_y = measure[i].shape
            assert(dim_x == dim_y and dim_x == dim_in * dim_s)

    def trace(self, rho):
        dim_x, dim_y = rho.shape
        assert(dim_x == dim_y and dim_x % self.dim_s == 0)

        res = np.matrix(np.zeros((self.dim_s, self.dim_s), dtype = complex64_type))
        for x in range(self.dim_s):
            for y in range(self.dim_s):
                res[x, y] = sum(rho[x + d, y + d] for d in range(0, dim_x, self.dim_s))
        return res

    def apply(self, rho, sigma, a):
        dim_x, dim_y = rho.shape
        assert(dim_x == dim_y and dim_x == self.dim_s)
        
        dim_x, dim_y = sigma.shape
        assert(dim_x == dim_y and dim_x == self.dim_in)

        return self.trace(self.measure[a] * self.unitary * np.kron(sigma, rho) * self.unitary.H * self.measure[a].H)

def check(m1, rho_in1, m2, rho_in2, basis, equal):
    assert(m1.dim_in == m2.dim_in and m1.gamma == m2.gamma)
    
    dim_x, dim_y = rho_in1.shape
    assert(dim_x == dim_y and dim_x == m1.dim_s)

    dim_x, dim_y = rho_in2.shape
    assert(dim_x == dim_y and dim_x == m2.dim_s)

    for x in basis:
        dim_x, dim_y = x.shape
        assert(dim_x == dim_y and dim_x == m1.dim_in)

    dim_in, gamma = m1.dim_in, m1.gamma
    dim = len(basis)
    q = queue.Queue()
    q.put([[], [], rho_in1, rho_in2])
    ortho = Orthogonal(equal)
    while q.qsize() > 0:
        pi, pa, rho1, rho2 = q.get()
        rho = np.concatenate((rho1.A1, -rho2.A1), axis = None)
        if ortho.isIndependent(rho):
            if not equal(np.trace(rho1), np.trace(rho2)):
                return [False, pi, pa, rho1, rho2]
            ortho.add(rho)
            for sigma in range(dim):
                for a in range(gamma):
                    q.put([pi + [sigma], pa + [a], m1.apply(rho1, basis[sigma], a), m2.apply(rho2, basis[sigma], a)])
    return [True, [], [], rho_in1, rho_in2]

equal = lambda x, y: abs(x - y) < 1e-8

I = np.matrix([[1, 0], [0, 1]], dtype = complex64_type)
H = np.matrix([[1, 1], [1, -1]], dtype = complex64_type) / 2**0.5
T = np.matrix([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype = complex64_type)
state = [np.matrix([[1, 0], [0, 0]], dtype = complex64_type), np.matrix([[0, 0], [0, 1]], dtype = complex64_type)]

def direct_sum(a, b):
    na, ma = a.shape
    nb, mb = b.shape
    c = np.matrix(np.zeros((na + nb, ma + mb), dtype = complex64_type))
    for x in range(na):
        for y in range(ma):
            c[x, y] = a[x, y]
    for x in range(nb):
        for y in range(mb):
            c[x + na, y + ma] = b[x, y]
    return c

def CNOT(dim, ctrl, targ):
    res = np.matrix(np.zeros((2 ** dim, 2 ** dim), dtype = complex64_type))
    for i in range(2 ** dim):
        for j in range(2 ** dim):
            if i == (j ^ (((j>>(dim-ctrl))&1) << (dim-targ))):
                res[i, j] = 1
    return res

def Toffoli(dim, ctrl1, ctrl2, targ):
    res = np.matrix(np.zeros((2 ** dim, 2 ** dim), dtype = complex64_type))
    for i in range(2 ** dim):
        for j in range(2 ** dim):
            if i == (j ^ ( (((j>>(dim-ctrl1))&1) & ((j>>(dim-ctrl2))&1)) << (dim-targ))):
                res[i, j] = 1
    return res

def testcase1():
    dim_in = 2
    dim_s = 8
    C = H
    shift0 = np.matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]], dtype = complex64_type)
    shift1 = np.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]], dtype = complex64_type)
    S = direct_sum(shift0, shift1)
    detection = Toffoli(4, 3, 4, 1)
    unitary = detection * np.kron(I, S) * np.kron(np.kron(np.kron(I, C), I), I)
    measure = [np.kron(np.kron(np.kron(state[0], I), I), I), np.kron(np.kron(np.kron(state[1], I), I), I)]
    machine = QuantumMealyMachineWithQuantumInput(dim_in, dim_s, unitary, measure)
    rho1 = np.kron(np.kron(state[0], state[0]), state[0])
    rho2 = np.kron(np.kron(state[0], state[1]), state[0])
    basis = [
            np.matrix([[1, 0], [0, 0]], dtype = complex64_type),
            np.matrix([[0, 0], [0, 1]], dtype = complex64_type),
            np.matrix([[0.5, 0.5], [0.5, 0.5]], dtype = complex64_type),
            np.matrix([[0.5, -0.5j], [0.5j, 0.5]], dtype = complex64_type)
        ]
    flag, pi, pa, rho1, rho2 = check(machine, rho1, machine, rho2, basis, equal)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

def tensor_product(a):
    res = a[0]
    for i in range(1, len(a)):
        res = np.kron(res, a[i])
    return res

def testcase_huge_walk(n):
    # n <= 6
    # testcase1: n = 2
    dim_in = 2
    dim_s = 2 ** (n + 1)
    C = H
    shift0 = np.matrix(np.zeros((2 ** n, 2 ** n), dtype = complex64_type))
    for x in range(2 ** n):
        shift0[(x + 1) % (2 ** n), x] = 1
    shift1 = np.matrix(np.zeros((2 ** n, 2 ** n), dtype = complex64_type))
    for x in range(2 ** n):
        shift1[(x - 1) % (2 ** n), x] = 1
    S = direct_sum(shift0, shift1)
    detection = np.matrix(np.zeros((2 ** (n + 2), 2 ** (n + 2)), dtype = complex64_type))
    for x in range(2 ** (n + 2)):
        y = x
        if (x & (2 ** n - 1)) == (2 ** n - 1): y ^= 2 ** (n + 1)
        detection[y, x] = 1
    unitary = detection * np.kron(I, S) * tensor_product([I, C] + [I for i in range(n)])
    measure = [tensor_product([state[0]] + [I for i in range(n + 1)]), tensor_product([state[1]] + [I for i in range(n + 1)])]
    machine = QuantumMealyMachineWithQuantumInput(dim_in, dim_s, unitary, measure)
    rho1 = tensor_product([state[0] for i in range(n + 1)]) # 000..0
    rho2 = tensor_product([state[0] if i == n else state[1] for i in range(n + 1)]) # 010...0
    basis = [
            np.matrix([[1, 0], [0, 0]], dtype = complex64_type),
            np.matrix([[0, 0], [0, 1]], dtype = complex64_type),
            np.matrix([[0.5, 0.5], [0.5, 0.5]], dtype = complex64_type),
            np.matrix([[0.5, -0.5j], [0.5j, 0.5]], dtype = complex64_type)
        ]
    time_start = time.time()
    flag, pi, pa, rho1, rho2 = check(machine, rho1, machine, rho2, basis, equal)
    time_end = time.time()
    print("time: ", time_end - time_start)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

def testcase_huge_qft(n):
    FT = np.matrix(np.zeros((2 ** n, 2 ** n), dtype = complex64_type))
    for j in range(2 ** n):
        for k in range(2 ** n):
            FT[j, k] = np.exp(2j * np.pi * j * k / (2 ** n)) / (2 ** (n/2))
    unitary1 = CNOT(n+2, 5, 1) * np.kron(I, direct_sum(FT, tensor_product([T] + [I for i in range(n - 1)])))
    #unitary2 = CNOT(n+2, n+2, 1) * np.kron(I, direct_sum(FT, tensor_product([I for i in range(n)])))
    measure = []
    for a in range(2):
        for b in range(2):
            measure.append(tensor_product([state[a], state[b]] + [I for i in range(n)]))
    machine1 = QuantumMealyMachineWithQuantumInput(4, 2 ** n, unitary1, measure)
    #machine2 = QuantumMealyMachineWithQuantumInput(4, 2 ** n, unitary2, measure)
    basis = [tensor_product([state[0], state[0]]), tensor_product([state[0], state[1]])]
    rho1 = tensor_product([state[0] for i in range(n)]) # 000..0
    rho2 = tensor_product([state[0] if i != 0 else state[1] for i in range(n)]) # 000..0
    time_start = time.time()
    flag, pi, pa, rho1, rho2 = check(machine1, rho1, machine1, rho2, basis, equal)
    time_end = time.time()
    print("time: ", time_end - time_start)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

def testcase7():
    H1 = np.kron(np.kron(H, I), I)
    H2 = np.kron(np.kron(I, H), I)
    H3 = np.kron(np.kron(I, I), H)
    U = CNOT(6, 4, 1) * np.kron(I, direct_sum(direct_sum(direct_sum(H1, H2), H3), Toffoli(3, 1, 2, 3)))
    basis = []
    for a in range(2):
        for b in range(2):
            basis.append(np.kron(np.kron(state[0], state[a]), state[b]))
    measure = []
    for a in range(2):
        for b in range(2):
            for c in range(2):
                measure.append(np.kron(np.kron(np.kron(state[a], state[b]), state[c]), np.eye(8, dtype = complex64_type)))
    machine = QuantumMealyMachineWithQuantumInput(8, 8, U, measure)
    rho1 = np.kron(np.kron(state[0], state[0]), state[0])
    rho2 = np.kron(np.kron(state[0], state[0]), state[1])
    time_start = time.time()
    flag, pi, pa, rho1, rho2 = check(machine, rho1, machine, rho2, basis, equal)
    time_end = time.time()
    print("time: ", time_end - time_start)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

def testcase8():
    H1 = np.kron(np.kron(H, I), I)
    H2 = np.kron(np.kron(I, H), I)
    H3 = np.kron(np.kron(I, I), H)
    U = CNOT(6, 4, 1) * np.kron(I, direct_sum(direct_sum(direct_sum(H1, H2), H3), Toffoli(3, 1, 2, 3)))
    basis = []
    for a in range(2):
        for b in range(2):
            basis.append(np.kron(np.kron(state[0], state[a]), state[b]))
    measure = []
    for a in range(2):
        for b in range(2):
            for c in range(2):
                measure.append(np.kron(np.kron(np.kron(state[a], state[b]), state[c]), np.eye(8, dtype = complex64_type)))
    machine = QuantumMealyMachineWithQuantumInput(8, 8, U, measure)
    rho1 = np.kron(np.kron(state[1], state[1]), state[0])
    rho2 = np.kron(np.kron(state[1], state[0]), state[1])
    time_start = time.time()
    flag, pi, pa, rho1, rho2 = check(machine, rho1, machine, rho2, basis, equal)
    time_end = time.time()
    print("time: ", time_end - time_start)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

def testcase9():
    basis = [np.kron(state[0], state[0]), np.kron(state[0], state[1])]
    FT = np.matrix(np.zeros((16, 16), dtype = complex64_type))
    for i in range(16):
        for j in range(16):
            FT[i, j] = np.exp(1j * np.pi * i * j / 8) / 4
    measure = []
    for a in range(2):
        for b in range(2):
            measure.append(np.kron(np.kron(state[a], state[b]), np.eye(16, dtype = complex64_type)))
    T1 = np.kron(T, np.eye(8, dtype = complex64_type))
    U = CNOT(6, 5, 1) * np.kron(I, direct_sum(FT, T1))
    machine = QuantumMealyMachineWithQuantumInput(4, 16, U, measure)
    rho1 = np.kron(np.kron(np.kron(state[0], state[0]), state[0]), state[0])
    rho2 = np.kron(np.kron(np.kron(state[1], state[0]), state[0]), state[0])
    time_start = time.time()
    flag, pi, pa, rho1, rho2 = check(machine, rho1, machine, rho2, basis, equal)
    time_end = time.time()
    print("time: ", time_end - time_start)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

def testcase10():
    basis = [np.kron(state[0], state[0]), np.kron(state[0], state[1])]
    FT = np.matrix(np.zeros((16, 16), dtype = complex64_type))
    for i in range(16):
        for j in range(16):
            FT[i, j] = np.exp(1j * np.pi * i * j / 8) / 4
    measure = []
    for a in range(2):
        for b in range(2):
            measure.append(np.kron(np.kron(state[a], state[b]), np.eye(16, dtype = complex64_type)))
    T1 = np.kron(T, np.eye(8, dtype = complex64_type))
    U = CNOT(6, 5, 1) * np.kron(I, direct_sum(FT, T1))
    machine = QuantumMealyMachineWithQuantumInput(4, 16, U, measure)
    rho1 = np.kron(np.kron(np.kron(state[0], state[0]), state[0]), state[0])
    rho2 = np.kron(np.kron(np.kron(state[0], state[1]), state[0]), state[0])
    time_start = time.time()
    flag, pi, pa, rho1, rho2 = check(machine, rho1, machine, rho2, basis, equal)
    time_end = time.time()
    print("time: ", time_end - time_start)
    print(flag)
    if not flag:
        print(pi)
        print(pa)

#testcase1()
#testcase9()
#print(8)
#testcase8()

#testcase_huge_walk(3) # n = 3, 4, 5, 6

#testcase_huge_qft(5) # n = 5, 6, 7

eval("testcase10()")