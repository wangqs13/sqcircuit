# sqcircuit

This is the source code for checking equivalence of sequential quantum circuits [WY18].

## Testcases

There are 12 testcases provided in "main.cpp" with respect to [WY18].

## How to use

Please according to "main.cpp".

We define a quantum Mealy machine with input Hilbert space of dimension ```dim_in``` and state Hilbert space of dimension ```dim_s``` by
```
QuantumMealyMachineWithQuantumInput<Real> machine(dim_in, dim_s, unitary, measure);
```
Here, ```unitary``` is a unitary operator that describes the combinational part of a sequential quantum circuit, and ```measure``` is a description of the measurement.

After that, a call to ```check(machine1, rho1, machine2, rho2, basis, equal)``` will return a result indicating whether the two machines with their initial states, respectively, are equivalent or not, and, if not, a witness (that is, an input/output pair) will be given.

To run the code, please make sure all headers and "main.cpp" are in the same directory. Enter
```
g++ main.cpp -o main -std=c++11
```
on the command line to compile the code, and then run it by entering
```
./main
```
on the command line. 

## References

[WY18] Q. S. Wang and M. S. Ying. Equivalence checking of sequential quantum circuits. In: [arXiv:1811.07722](https://arxiv.org/abs/1811.07722).
