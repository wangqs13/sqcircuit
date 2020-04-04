#include <bits/stdc++.h>

#include "complex.h"
#include "matrix.h"
#include "almost.h"
#include "orthogonal.h"
#include "qmmwqi.h"

// using namespace std;

using std::cout;
using std::endl;

typedef double Real;
typedef Complex<Real> complex;
typedef Matrix<Real, complex> matrix;

const complex imag(0, 1);
const Real pi = 3.1415926535897932384626433;

const matrix I = matrix::eye(2);
const matrix O = matrix::zeros(2);
const matrix X = { {0, 1}, {1, 0} };
const matrix Y = { {0, -imag}, {imag, 0} };
const matrix Z = { {1, 0}, {0, -1} };
const matrix H = { {sqrt(2.0) / 2, sqrt(2.0) / 2}, {sqrt(2.0) / 2, -sqrt(2.0) / 2} };
const matrix T = { {1, 0}, {0, exp(imag * pi / 4)} };
const matrix CNOT = direct_sum(I, X);
const matrix SWAP = { {1, 0, 0, 0},
                     {0, 0, 1, 0},
                     {0, 1, 0, 0},
                     {0, 0, 0, 1} };
const matrix Toffoli = direct_sum(I, I, I, X);

std::vector<matrix> state = { matrix{{1, 0}, {0, 0}}, matrix{{0, 0}, {0, 1}} };

AlmostEqual<Real> equal(1e-8);

void testCase1()
{
    int dim_in = 2;
    int dim_s = 8;
    matrix C = H;
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
    matrix rho1 = tensor_product(state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[0], state[1], state[0]);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    printf("Case #%d\n", 1);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase2()
{
    int dim_in = 2;
    int dim_s = 8;
    matrix C = H;
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    matrix rho1 = tensor_product(basis[3], state[0], state[0]);
    matrix rho2 = tensor_product(basis[3], state[1], state[0]);
    printf("Case #%d\n", 2);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase3()
{
    int dim_in = 2;
    int dim_s = 8;
    matrix C = H;
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    matrix rho1 = tensor_product(state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[1], state[1], state[0]);
    printf("Case #%d\n", 3);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase4()
{
    int dim_in = 2;
    int dim_s = 8;
    matrix C = H;
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    matrix rho1 = tensor_product(state[0], state[0], state[1]);
    matrix rho2 = tensor_product(state[1], state[0], state[1]);
    printf("Case #%d\n", 4);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase5()
{
    int dim_in = 2;
    int dim_s = 8;
    matrix C = H;
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    matrix rho1 = tensor_product(state[0], state[1], state[0]);
    matrix rho2 = tensor_product(state[1], state[0], state[0]);
    printf("Case #%d\n", 5);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase6()
{
    int dim_in = 2;
    int dim_s = 8;
    matrix C = H;
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine1(dim_in, dim_s, unitary, measure);
    C = { {sqrt(2.0) / 2, sqrt(2.0) / 2 * imag}, {sqrt(2.0) / 2 * imag, sqrt(2.0) / 2} };
    unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    QuantumMealyMachineWithQuantumInput<Real> machine2(dim_in, dim_s, unitary, measure);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    matrix rho1 = tensor_product(state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[0], state[0], state[0]);
    printf("Case #%d\n", 6);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine1, rho1, machine2, rho2, basis, equal);
    auto ret = check(machine1, rho1, machine2, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase7()
{
    matrix H1 = tensor_product(H, I, I);
    matrix H2 = tensor_product(I, H, I);
    matrix H3 = tensor_product(I, I, H);
    matrix U = matrix::CNOT(6, 4, 1) * tensor_product(I, direct_sum(H1, H2, H3, Toffoli));
    std::vector<matrix> basis;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
        {
            basis.push_back(tensor_product(state[0], state[a], state[b]));
        }
    std::vector<matrix> measure;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
            for (int c = 0; c <= 1; ++c)
            {
                measure.push_back(tensor_product(state[a], state[b], state[c], I, I, I));
            }
    QuantumMealyMachineWithQuantumInput<Real> machine = { 8, 8, U, measure };
    matrix rho1 = tensor_product(state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[0], state[0], state[1]);
    printf("Case #%d\n", 7);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase8()
{
    matrix H1 = tensor_product(H, I, I);
    matrix H2 = tensor_product(I, H, I);
    matrix H3 = tensor_product(I, I, H);
    matrix U = matrix::CNOT(6, 4, 1) * tensor_product(I, direct_sum(H1, H2, H3, Toffoli));
    std::vector<matrix> basis;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
        {
            basis.push_back(tensor_product(state[0], state[a], state[b]));
        }
    std::vector<matrix> measure;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
            for (int c = 0; c <= 1; ++c)
            {
                measure.push_back(tensor_product(state[a], state[b], state[c], I, I, I));
            }
    QuantumMealyMachineWithQuantumInput<Real> machine = { 8, 8, U, measure };
    matrix rho1 = tensor_product(state[1], state[1], state[0]);
    matrix rho2 = tensor_product(state[1], state[0], state[1]);
    printf("Case #%d\n", 8);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase9()
{
    std::vector<matrix> basis;
    basis.push_back(tensor_product(state[0], state[0]));
    basis.push_back(tensor_product(state[0], state[1]));
    matrix FT = matrix(16, 16, [](int i, int j)
        {
            return exp(imag * pi * i * j / 8) / (4);
        });
    std::vector<matrix> measure;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
        {
            measure.push_back(tensor_product(state[a], state[b], I, I, I, I));
        }
    matrix T1 = tensor_product(T, I, I, I);
    matrix U = matrix::CNOT(6, 5, 1) * tensor_product(I, direct_sum(FT, T1));
    QuantumMealyMachineWithQuantumInput<Real> machine = { 4, 16, U, measure };
    matrix rho1 = tensor_product(state[0], state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[1], state[0], state[0], state[0]);
    printf("Case #%d\n", 9);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase10()
{
    std::vector<matrix> basis;
    basis.push_back(tensor_product(state[0], state[0]));
    basis.push_back(tensor_product(state[0], state[1]));
    matrix FT = matrix(16, 16, [](int i, int j)
        {
            return exp(imag * pi * i * j / 8) / (4);
        });
    std::vector<matrix> measure;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
        {
            measure.push_back(tensor_product(state[a], state[b], I, I, I, I));
        }
    matrix T1 = tensor_product(T, I, I, I);
    matrix U = matrix::CNOT(6, 5, 1) * tensor_product(I, direct_sum(FT, T1));
    QuantumMealyMachineWithQuantumInput<Real> machine = { 4, 16, U, measure };
    matrix rho1 = tensor_product(state[0], state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[0], state[1], state[0], state[0]);
    printf("Case #%d\n", 10);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine, rho1, machine, rho2, basis, equal);
    auto ret = check(machine, rho1, machine, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase11()
{
    auto Rx = [&](Real theta) -> matrix
    {
        return { {cos(theta / 2), -imag * sin(theta / 2)},{-imag * sin(theta / 2), cos(theta / 2)} };
    };
    int dim_in = 2;
    int dim_s = 8;
    matrix C = Rx(sqrt(2));
    matrix shift0
    {
        {0, 0, 0, 1},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0}
    };
    matrix shift1
    {
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 0, 0}
    };
    matrix S = direct_sum(shift0, shift1);
    matrix detection = matrix::Toffoli(4, 3, 4, 1);
    matrix unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    std::vector<matrix> measure = { tensor_product(state[0], I, I, I), tensor_product(state[1], I, I, I) };
    QuantumMealyMachineWithQuantumInput<Real> machine1(dim_in, dim_s, unitary, measure);
    C = Rx(-sqrt(2));
    unitary = detection * tensor_product(I, S) * tensor_product(I, C, I, I);
    QuantumMealyMachineWithQuantumInput<Real> machine2(dim_in, dim_s, unitary, measure);
    std::vector<matrix> basis = {
            {{1, 0}, {0, 0}},
            {{0, 0}, {0, 1}},
            {{0.5, 0.5}, {0.5, 0.5}},
            {{0.5, -0.5 * imag}, {0.5 * imag, 0.5}}
    };
    matrix rho1 = tensor_product(state[0], state[0], state[0]);
    matrix rho2 = tensor_product(state[0], state[0], state[0]);
    printf("Case #%d\n", 11);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine1, rho1, machine2, rho2, basis, equal);
    auto ret = check(machine1, rho1, machine2, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

void testCase12()
{
    matrix V = { {(1 + imag) / 2, (1 - imag) / 2}, {(1 - imag) / 2, (1 + imag) / 2} };
    matrix VH = V.hermitian();
    auto controlledGate = [](int dim, int ctrl, int targ, const matrix& u) -> matrix
    {
        assert(u.col == 2 && u.row == 2);
        assert(ctrl != targ);
        assert(1 <= ctrl && ctrl <= dim);
        assert(1 <= targ && targ <= dim);
        return matrix(1 << dim, 1 << dim, [&](int i, int j) -> complex
            {
                if (j >> (dim - ctrl) & 1)
                {
                    if (i >> (dim - ctrl) & 1)
                    {
                        for (int x = 1; x <= dim; ++x) if (x != ctrl && x != targ)
                            if ((j >> (dim - x) & 1) != (i >> (dim - x) & 1)) return 0;
                        return u[i >> (dim - targ) & 1][j >> (dim - targ) & 1];
                    }
                    else
                        return 0;
                }
                else
                {
                    if (i >> (dim - ctrl) & 1)
                        return 0;
                    else
                    {
                        return i == j;
                    }
                }
            });
    };
    matrix U1 = controlledGate(3, 2, 1, X) * controlledGate(3, 2, 1, V) * controlledGate(3, 3, 1, V) * controlledGate(3, 2, 3, X) * controlledGate(3, 3, 1, VH);
    matrix U2 = controlledGate(3, 2, 3, X) * matrix::Toffoli(3, 2, 3, 1);
    std::vector<matrix> basis;
    for (int a = 0; a < 2; ++a)
    {
        basis.push_back(tensor_product(state[0], state[a]));
    }
    std::vector<matrix> measure;
    for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
        {
            measure.push_back(tensor_product(state[a], state[b], I));
        }
    QuantumMealyMachineWithQuantumInput<Real> machine1 = { 4, 2, U1, measure };
    QuantumMealyMachineWithQuantumInput<Real> machine2 = { 4, 2, U2, measure };
    matrix rho1 = state[0];
    matrix rho2 = state[0];
    printf("Case #%d\n", 12);
    int64_t tmp = clock();
    //auto [flag, witness] = check(machine1, rho1, machine2, rho2, basis, equal);
    auto ret = check(machine1, rho1, machine2, rho2, basis, equal);
    auto flag = ret.first;
    auto witness = ret.second;
    printf("Time cost: %.3fs\n", double(clock() - tmp) / CLOCKS_PER_SEC);
    if (flag)
    {
        printf("Verdict: Yes.\n");
    }
    else
    {
        printf("Verdict: No.\n");
        std::cout << witness;
    }
    puts("____________________");
}

int main()
{
    testCase1();
    testCase2();
    testCase3();
    testCase4();
    testCase5();
    testCase6();
    testCase7();
    testCase8();
    testCase9();
    testCase10();
    testCase11();
    testCase12();

    return 0;
}
