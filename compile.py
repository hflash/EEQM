# This code is part of LINKEQ.
#
# (C) Copyright LINKE 2023.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
#
# -*- coding: utf-8 -*-
# @Time    : 2024/3/31 13:36
# @Author  : HFALSH @ LINKE
# @File    : compile.py
# @IDE     : PyCharm
import os
import re
import sys
import subprocess
import numpy as np
from quantumcircuit import QuantumCircuit
#from quantumcircuit.gate import *
from qiskit.transpiler.passes.layout import VF2Layout
from qiskit.transpiler import CouplingMap
from qiskit import QuantumCircuit as qiskitqc
from qiskit.converters import circuit_to_dag

from circuit_slices import get_gates


def initial_map(QASMfile:str, chipfile):
    coupling = []
    with open(chipfile) as f:
        lines = f.readlines()
        for i in range(1,len(lines)):
            q1,q2 = [int(x) for x in lines[i].strip().split()]
            coupling.append([q1,q2])
            coupling.append([q2,q1])
    coupling_map = CouplingMap(coupling)
    layout_instance = VF2Layout(coupling_map=coupling_map)
    qc = QuantumCircuit.from_QASM(QASMfile)
    dagtable = np.array(qc.to_dagtable())
    circuit = qiskitqc(qc.get_qubit_number())
    nowlayer = 0
    initial_map = [i for i in range(qc.get_qubit_number())]
    greedy_layout = None
    while nowlayer<dagtable.shape[1]:
        _,cnot_gates = get_gates(dagtable[:,nowlayer])
        if len(cnot_gates) != 0:
            for gate in cnot_gates:
                circuit.cz(gate[0], gate[1])
            layout_instance.run(circuit_to_dag(circuit))
            result = layout_instance.property_set['VF2Layout_stop_reason']
            if  result.value == "solution found":
                greedy_layout = layout_instance.property_set['layout']
            else:
                break
        nowlayer+=1
    if greedy_layout:
        for device_q,circuit_q in greedy_layout._v2p.items():
            initial_map[circuit_q] = device_q._index
    initial_mapping_str = '['
    for qubit in initial_map:
        initial_mapping_str += str(qubit)
        if qubit == initial_map[-1]:
            continue
        else:
            initial_mapping_str += ','
    initial_mapping_str += ']'
    return initial_mapping_str

# compile for SWin+
def compile(circuit_path, chip_path, result_path):
    circuit_path = os.path.abspath(circuit_path)
    chip_path = os.path.abspath(chip_path)
    result_path = os.path.abspath(result_path)
    cpp_program = "/home/edge/hflash/EEQM_t/cmake-build-debug/EEQM_scale"
    initial_mapping_str = initial_map(circuit_path, chip_path)

    result = subprocess.run([cpp_program, circuit_path, chip_path, result_path, "nogreedy", initial_mapping_str], capture_output=True, text=True)
    output = result.stdout

    if output:
        compiler_time = None
        compiler_depth = None
        lines = output.split('\n')
        for i in range(len(lines)):
            if "time = " in lines[i]:
                match = re.search(r'\d+\.\d+', lines[i])
                compiler_time = float(match.group())
            elif "cnot number = " in lines[i]:
                match = re.search(r'\d+', lines[i])
                extra_cnot_cnt = int(match.group())
            elif "depth = " in lines[i]:
                match = re.search(r'\d+', lines[i])
                compiler_depth = int(match.group())
        compile_result = {
            "compiler_time": compiler_time,
            "initial_mapping": initial_mapping_str,
            "compiler_depth": compiler_depth,
            "swap_count": int(extra_cnot_cnt/3),
        }
    else:
        print([circuit_path, chip_path, result_path, initial_mapping_str])
    return compile_result

# compile for SWin
def compile_without_slice(circuit_path, chip_path, result_path, strategy, initial_mapping_str):
    circuit_path = os.path.abspath(circuit_path)
    chip_path = os.path.abspath(chip_path)
    result_path = os.path.abspath(result_path)
    cpp_program = "/home/edge/hflash/EEQM_t/cmake-build-debug/EEQM_scale"
    initial_mapping_str = initial_mapping_str

    result = subprocess.run([cpp_program, circuit_path, chip_path, result_path, strategy, initial_mapping_str],
                            capture_output=True, text=True)
    output = result.stdout

    if output:
        compiler_time = None
        compiler_depth = None
        lines = output.split('\n')
        for i in range(len(lines)):
            if "time = " in lines[i]:
                match = re.search(r'\d+\.\d+', lines[i])
                compiler_time = float(match.group())
            elif "cnot number = " in lines[i]:
                match = re.search(r'\d+', lines[i])
                extra_cnot_cnt = int(match.group())
            elif "depth = " in lines[i]:
                match = re.search(r'\d+', lines[i])
                compiler_depth = int(match.group())
        compile_result = {
            "compiler_time": compiler_time,
            "initial_mapping": initial_mapping_str,
            "compiler_depth": compiler_depth,
            "swap_count": int(extra_cnot_cnt / 3),
        }
    else:
        print([circuit_path, chip_path, result_path, initial_mapping_str])
    return compile_result


# batch test SWin for small-scale circuits
# the window size are set as 8 in main func of C++ project
def batch_compile_circuit_small():
    circuit_path_all_small = '/home/edge/fzchen/swin_src/qasm-benchmark/cr_iccad_circuits/small'
    circuit_path_all_large = '/home/edge/fzchen/swin_src/qasm-benchmark/cr_iccad_circuits/large'
    # circuit_path_all_large = ''
    chip_paths = ["/home/edge/hflash/EEQM_t/couplings/t.txt", "/home/edge/hflash/EEQM_t/couplings/guadalupe.txt"]
    result_path = 'result.txt'
    strategies_small = ['search', 'greedy', 'nogreedy']
    strategies_large = ['search', 'greedy', 'nogreedy', 'fidelity']
    initial_mapping_strs = [str([i for i in range(5)]), str([i for i in range(16)])]
    # for root, dirs, files in os.walk(circuit_path_all_small):
    #     for file in files:
    #         for strategy in strategies_small:
    #             circuit_path = os.path.join(root, file)
    #             compile_result = compile_without_slice(circuit_path, chip_paths[0], result_path, strategy, initial_mapping_strs[0])
    #             print("filename->", file)
    #             print("strategy->", strategy)
    #             print("result->", compile_result)
    for root, dirs, files in os.walk(circuit_path_all_large):
        for file in files:
            for strategy in strategies_small:
                circuit_path = os.path.join(root, file)
                compile_result = compile_without_slice(circuit_path, chip_paths[1], result_path, strategy, initial_mapping_strs[1])
                print("filename->", file)
                print("strategy->", strategy)
                print("result->", compile_result)

def batch_compile_circuit_large():
    circuit_path_all_large = '/home/edge/fzchen/swin_src/qasm-benchmark/cr_iccad_circuits/large'
    # circuit_path_all_large = ''
    chip_paths = ["/home/edge/hflash/EEQM_t/couplings/t.txt", "/home/edge/hflash/EEQM_t/couplings/guadalupe.txt"]
    result_path = 'result.txt'
    strategies_large = ['search', 'greedy', 'nogreedy', 'fidelity']
    initial_mapping_strs = [str([i for i in range(5)]), str([i for i in range(16)])]
    for root, dirs, files in os.walk(circuit_path_all_large):
        for file in files:
            for strategy in strategies_large:
                circuit_path = os.path.join(root, file)
                compile_result = compile_without_slice(circuit_path, chip_paths[1], result_path, strategy, initial_mapping_strs[1])
                print("filename->", file)
                print("strategy->", strategy)
                print("result->", compile_result)

if __name__ == "__main__":
    batch_compile_circuit_small()
    batch_compile_circuit_large()
    # if len(sys.argv) != 4:
    #     exit(0)
    # circuit = sys.argv[1]
    # chip = sys.argv[2]
    # result_path = sys.argv[3]
    # compile_result = compile(circuit, chip, result_path)
    # print(compile_result)