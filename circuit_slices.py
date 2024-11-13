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
# @Time    : 2023/8/1 09:36
# @Author  : HFALSH @ LINKE
# @File    : circuit_slices.py
# @IDE     : PyCharm
import os
import numpy as np
from quantumcircuit import QuantumCircuit

def get_gates(layer):
    cnot_gates = []
    signle_gates = []
    have_checked = [0 for i in range(len(layer))]
    for i in range(len(layer)):
        if layer[i] == -1 or have_checked[i]:
            continue
        is_signle = True
        for j in range(i+1, len(layer)):
            if layer[i] == layer[j]:
                have_checked[j] = 1
                is_signle = False
                cnot_gates.append([i, j, layer[i]])
                break
        if is_signle:
            signle_gates.append([i, layer[i]])
    return signle_gates, cnot_gates


def cut_circuit_depth(file_path, write_path, slice_depth):
    """
    Args:
        file_path: read from qasm file
        write_path: write to a directory
        slice_depth: the depth of each circuit slice

    Returns:
        None

    The circuit slices is written to files in write_path
    """
    circuit = QuantumCircuit.from_QASM(file_path)
    filename = file_path.split('/')[-1].split('.qasm')[0]
    circuits_slices = []
    dag_table = circuit.to_dagtable()
    remaining_depth = circuit.get_circuit_depth() - 1
    """
        Todo: circuit depth update
    """
    slice_num = int(circuit.get_circuit_depth() / slice_depth) + 1
    for slice in range(slice_num):
        circuit_slice = QuantumCircuit(qubit_number=circuit.qubit_number)
        for i in range(slice_depth):
            if remaining_depth <= 0:
                break
            # print(remaining_depth)
            # print(i + slice * slice_depth)
            gate_time = set([row[i + slice * slice_depth] for row in dag_table])
            "set(), for the cnot gate id may occur twice."
            "gate_id '-1' indicate none gate in gate list while gate_list[-1] means the last gate."
            # print(gate_time)
            for gate_id in gate_time:
                # print(gate_id, end=' ')
                if gate_id == -1:
                    continue
                gate = circuit.gate_list[gate_id]
                circuit_slice.add_gate(gate)
                # print(circuit_slice.get_circuit_depth())
            remaining_depth -= 1
        # print('------------------')
        circuits_slices.append(circuit_slice)
    write_file_path = os.path.join(write_path, filename)
    if os.path.exists(write_file_path):
        pass
    else:
        os.mkdir(write_file_path)
    files_list = []
    for i in range(slice_num):
        filename_slice = os.path.join(write_file_path, filename+'_slice_'+str(i)+'.qasm')
        files_list.append(filename_slice)
        # print(filename_slice)
        # print(circuits_slices[i].get_circuit_depth())
        circuits_slices[i].to_QASM(filename_slice)
        # print(circuits_slices[i].to_dagtable())
    return files_list

def cut_circuit_gate(file_path, write_path, slices):
    circuit = QuantumCircuit.from_QASM(file_path)
    dagtbale = np.array(circuit.to_dagtable())
    total_gate_cnt = np.sum(np.where(dagtbale>0, 1, 0))
    average_slice_gates = total_gate_cnt/slices


    filename = file_path.split('/')[-1].split('.qasm')[0]
    circuits_slices = []
    nowlayer = 0
    single_gates = []
    cnot_gates = []
    bias = 0
    for slice in range(slices):
        circuit_slice = QuantumCircuit(qubit_number=circuit.qubit_number)
        tmp_gate_numbers = 0
        while tmp_gate_numbers<average_slice_gates-bias:
            if len(cnot_gates) != 0:
                gate = circuit.gate_list[cnot_gates.pop()[-1]]
                circuit_slice.add_gate(gate)
                tmp_gate_numbers += 2
            elif len(single_gates) != 0:
                gate = circuit.gate_list[single_gates.pop()[-1]]
                circuit_slice.add_gate(gate)
                tmp_gate_numbers += 1
            else:
                if nowlayer<dagtbale.shape[1]:
                    single_gates,cnot_gates = get_gates(dagtbale[:,nowlayer])
                    nowlayer += 1
                else:
                    break

        # print('------------------')
        if tmp_gate_numbers>0:
            circuits_slices.append(circuit_slice)
            bias += tmp_gate_numbers - average_slice_gates
    write_file_path = os.path.join(write_path, filename)
    if os.path.exists(write_file_path):
        pass
    else:
        os.mkdir(write_file_path)
    files_list = []
    for i in range(slices):
        filename_slice = os.path.join(write_file_path, filename+'_slice_'+str(i)+'.qasm')
        files_list.append(filename_slice)
        # print(filename_slice)
        # print(circuits_slices[i].get_circuit_depth())
        circuits_slices[i].to_QASM(filename_slice)
        # print(circuits_slices[i].to_dagtable())
    return files_list



def cut_circuit(file_path, write_path, method, args):
    methods = {
        'depth':cut_circuit_depth,
        'gatecnt':cut_circuit_gate
    }
    if method in methods:
        return methods[method](file_path, write_path, int(args))
    else:
        print('unknow cut method:', method)
        exit(0)





