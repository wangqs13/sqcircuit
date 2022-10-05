# sqcircuit

This is the source code for checking equivalence of sequential quantum circuits [WY18], [WLY22].

## Testcases

There are 12 testcases provided in "main.cpp" with respect to [WY18].

There are 20 testcases in OpenQASM 3 provided in [Veri-Q Benchmark](https://github.com/Veri-Q/Benchmark) with respect to [WLY22].

## How to use CPP code

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

Some improvements were later made, see the source code of [qmm-benchmark](https://github.com/wangqs13/qmm-benchmark).

## How to use Python code

There is only one file "qmm.py", which is the whole Python code that contains the algorithm and the data. 

The implementation is efficient, and was used in the experiments of [WLY22]. 

## References

[WY18] Qisheng Wang and Mingsheng Ying. Equivalence checking of sequential quantum circuits. In: [arXiv:1811.07722v1](https://arxiv.org/abs/1811.07722v1).

[WLY22] Qisheng Wang, Riling Li and Mingsheng Ying. Equivalence checking of sequential quantum circuits. *IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems*, 41(9): 3143-3156, 2022. [arXiv:1811.07722](https://arxiv.org/abs/1811.07722). [doi:10.1109/TCAD.2021.3117506](https://doi.org/10.1109/TCAD.2021.3117506).
