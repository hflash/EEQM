import os
import sys
import random
import numpy as np
from compile import compile
from circuit_slices import cut_circuit, get_gates
from quantumcircuit import QuantumCircuit
from matplotlib import pyplot as plt
from multiprocessing import Pool, pool


def depth_sample(circuit, depth, startlayer):
    dagtable = np.array(circuit.to_dagtable())
    circuit_slice = QuantumCircuit(qubit_number=dagtable.shape[0])
    layer_numbers = dagtable.shape[1]
    now_depth = 0
    while startlayer < layer_numbers and now_depth < depth:
        gate_time = set(dagtable[:, startlayer])
        for gate_id in gate_time:
            if gate_id == -1:
                continue
            gate = circuit.gate_list[gate_id]
            circuit_slice.add_gate(gate)
        startlayer += 1
        now_depth += 1
    if now_depth == depth:
        return circuit_slice
    else:
        return None


def gate_sample(circuit, average_gate_cnt, startlayer):
    dagtable = np.array(circuit.to_dagtable())
    circuit_slice = QuantumCircuit(qubit_number=dagtable.shape[0])
    layer_numbers = dagtable.shape[1]
    tmp_gate_numbers = 0
    single_gates = []
    cnot_gates = []
    while tmp_gate_numbers < average_gate_cnt:
        if len(cnot_gates) != 0:
            gate = circuit.gate_list[cnot_gates.pop()[-1]]
            circuit_slice.add_gate(gate)
            tmp_gate_numbers += 2
        elif len(single_gates) != 0:
            gate = circuit.gate_list[single_gates.pop()[-1]]
            circuit_slice.add_gate(gate)
            tmp_gate_numbers += 1
        else:
            if startlayer < layer_numbers:
                single_gates, cnot_gates = get_gates(dagtable[:, startlayer])
                startlayer += 1
            else:
                break
    if tmp_gate_numbers >= average_gate_cnt:
        return circuit_slice
    else:
        return None


def sample(circuit_path, slice_number, method, arg, write_path):
    circuit = QuantumCircuit.from_QASM(circuit_path)
    if method == 'depth':
        sample_func = depth_sample
    elif method == 'gatecnt':
        sample_func = gate_sample
    else:
        print('unknow sample method:', method)
        exit(0)
    write_file_path = os.path.join(write_path, method + 'sample')
    if not os.path.exists(write_file_path):
        os.makedirs(write_file_path)
    has_get_circuit_number = 0
    files_list = []
    while has_get_circuit_number < slice_number:
        startlayer = random.randint(0, circuit.get_circuit_depth() - 1)
        circuit_slice = sample_func(circuit, arg, startlayer)
        if circuit_slice is not None:
            slice_filename = os.path.join(write_file_path, 'slice_' + str(has_get_circuit_number) + '.qasm')
            circuit_slice.to_QASM(slice_filename)
            files_list.append(slice_filename)
            has_get_circuit_number += 1
    return files_list


class Slice_info:
    def __init__(self, circuit_slices, method, arg):
        self.circuit_slices = circuit_slices
        self.result_dir = os.path.join(os.path.dirname(circuit_slices[0]), 'result')
        # if not os.path.exists(self.result_dir):
        #    os.mkdir(self.result_dir)
        self.method = method
        self.arg = arg
        self.processes: list[pool.ApplyResult] = []
        self.data = []

    def get_circuit_slices(self):
        return self.circuit_slices

    def set_processes(self, processes: list[pool.ApplyResult]):
        self.processes = processes

    def get_data(self):
        if len(self.data) == 0:
            for p in self.processes:
                result = p.get()
                self.data.append([result["compiler_time"], result["compiler_depth"]])
        return self.data

    def get_slice_number(self):
        return len(self.circuit_slices)

    def draw(self):
        data = np.array(self.get_data())
        x = [str(i) for i in range(len(self.circuit_slices))]
        plt.clf()
        plt.bar(x, data[:, 0])
        plt.savefig(os.path.join(os.path.dirname(self.circuit_slices[0]), 'time.png'))
        plt.clf()
        plt.bar(x, data[:, 1])
        plt.savefig(os.path.join(os.path.dirname(self.circuit_slices[0]), 'depth.png'))


def run_parallel(slice_infoes: list[Slice_info], coupling_file):
    p = Pool(processes=64)
    for slice_info in slice_infoes:
        slice_info.set_processes([p.apply_async(compile, (
        item, coupling_file, os.path.join(slice_info.result_dir, os.path.basename(item)[:-5])))
                                  for item in slice_info.get_circuit_slices()])
    p.close()
    p.join()


if __name__ == '__main__':
    result_dir = '/home/edge/fzchen/swin_src/result/'
    qasm_file = sys.argv[1]
    coupling_file = sys.argv[2]
    device_depth = int(sys.argv[3])
    coupling_name = os.path.basename(coupling_file)

    result_path = result_dir + coupling_name

    qasm_name = os.path.basename(qasm_file)
    circuit = QuantumCircuit.from_QASM(qasm_file)
    dagtable = np.array(circuit.to_dagtable())
    circuit_depth = dagtable.shape[1]
    total_gate_cnt = np.sum(np.where(dagtable != -1, 1, 0))

    initial_slice_number = int(2 * circuit_depth / device_depth)
    ideal_slice_number = int(circuit_depth / device_depth) + 1
    initial_gate_cnt = int(total_gate_cnt / initial_slice_number)
    initial_depth = int(device_depth / 2)

    write_path = os.path.join(result_path + str(device_depth), qasm_name[:-5])
    if not os.path.exists(write_path):
        os.makedirs(write_path)
    qasm_files = sample(qasm_file, 20, 'depth', initial_depth, write_path)
    depth_sample_slice_info = Slice_info(qasm_files, 'depth', 20)
    qasm_files = sample(qasm_file, 20, 'gatecnt', initial_gate_cnt, write_path)
    gatecnt_sample_slice_info = Slice_info(qasm_files, 'gatecnt', 20)
    run_parallel([depth_sample_slice_info, gatecnt_sample_slice_info], coupling_file)
    # slice_info.draw()
    # depth_sample_slice_info.draw()
    # gatecnt_sample_slice_info.draw()

    depth_sample_data = np.array(depth_sample_slice_info.get_data())
    gatecnt_sample_data = np.array(gatecnt_sample_slice_info.get_data())
    depth_sample_time = np.max(depth_sample_data[:, 0])
    gatecnt_sample_time = np.max(gatecnt_sample_data[:, 0])
    sample_time = np.max([depth_sample_time, gatecnt_sample_time])

    depth_sample_cv = np.std(depth_sample_data[:, 1]) / np.mean(depth_sample_data[:, 1])
    gatecnt_sampe_cv = np.std(gatecnt_sample_data[:, 1]) / np.mean(gatecnt_sample_data[:, 1])
    if depth_sample_cv < gatecnt_sampe_cv:
        better_method = 'depth'
        step = 1
        max_depth = np.max(depth_sample_data[:, 1])
        bound_depth = device_depth * initial_depth / max_depth
        bias = 2 * abs(bound_depth - initial_depth)
        if bias < 50:
            bias = 50
        if max_depth > device_depth:
            low_bound_depth = max(int(bound_depth - bias), 1)
            if max_depth - device_depth < 0.1 * device_depth:
                up_bound_depth = min(int(initial_depth + 0.5 * bias), 100)
            else:
                up_bound_depth = initial_depth
        else:
            up_bound_depth = min(int(bound_depth + bias), 100)
            if device_depth - max_depth < 0.1 * device_depth:
                low_bound_depth = max(int(initial_depth - 0.5 * bias), 1)
            else:
                low_bound_depth = initial_depth
        up_bound_slice_number = int(circuit_depth / low_bound_depth) + 1
        low_bound_slice_number = int(circuit_depth / up_bound_depth) + 1
        arg_list = []
        bias = up_bound_slice_number - low_bound_slice_number
        if bias < 10:
            up_bound_slice_number += 5
            low_bound_slice_number = max(1, low_bound_slice_number - 5)
        low_bound_slice_number = max(ideal_slice_number, low_bound_slice_number)
        for i in range(low_bound_slice_number, up_bound_slice_number):
            arg_list.append(int(circuit_depth / i) + 1)
        arg_list = list(set(arg_list))
    else:
        better_method = 'gatecnt'
        step = 1
        max_depth = np.max(gatecnt_sample_data[:, 1])
        # device_depth/(max_depth/gate_cnt)
        bound_gate_cnt = device_depth * initial_gate_cnt / max_depth
        bound_slice_number = int(total_gate_cnt / bound_gate_cnt)
        bias = abs(bound_slice_number - initial_slice_number)
        if bias < 5:
            bias = 5
        if max_depth > device_depth:
            up_bound_slice_number = int(bound_slice_number + bias)
            if max_depth - device_depth < 0.1 * device_depth:
                low_bound_slice_number = max(int(initial_slice_number - 0.5 * bias), 1)
            else:
                low_bound_slice_number = initial_slice_number
        else:
            low_bound_slice_number = max(int(bound_slice_number - bias), 1)
            if device_depth - max_depth < 0.1 * device_depth:
                up_bound_slice_number = int(initial_slice_number + 0.5 * bias)
            else:
                up_bound_slice_number = initial_slice_number
        low_bound_slice_number = max(ideal_slice_number, low_bound_slice_number)
        arg_list = [i for i in range(low_bound_slice_number, up_bound_slice_number, step)]
    print("better method:", better_method)
    print("arg list:", arg_list)
    slice_info_list = []
    for arg in arg_list:
        write_path = os.path.join(result_path + str(device_depth), qasm_name[:-5], str(arg))
        if not os.path.exists(write_path):
            os.mkdir(write_path)
        qasm_files = cut_circuit(qasm_file, write_path, better_method, arg)
        slice_info_list.append(Slice_info(qasm_files, better_method, arg))
    run_parallel(slice_info_list, coupling_file)
    log_file = write_path = os.path.join(result_path + str(device_depth), qasm_name[:-5], 'log')

    min_slice_number = 10000
    min_total_depth = 100000
    best_slice_max_depth = 10000
    best_slice_max_time = 10000
    each_slice_number_max_time = []
    with open(log_file, 'w') as f:
        f.write("better method: {} max depth: {} sample time: {}\n".format(better_method, max_depth, sample_time))
        for item in slice_info_list:
            # item.draw()
            data = item.get_data()
            npdata = np.array(data)
            slice_number = item.get_slice_number()
            total_depth = np.sum(npdata[:, 1])
            max_depth = np.max(npdata[:, 1])
            max_time = np.max(npdata[:, 0])
            each_slice_number_max_time.append(max_time)
            f.write("slice_number: {}  max depth: {} max time: {} total depth: {}\n".format(slice_number, max_depth,
                                                                                            max_time, total_depth))
            # f.write(str(data)+'\n')
            if max_depth <= device_depth:
                if slice_number < min_slice_number or (
                        slice_number == min_slice_number and total_depth < min_total_depth):
                    min_slice_number = slice_number
                    min_total_depth = total_depth
                    best_slice_max_depth = max_depth
                    best_slice_max_time = max_time
        max_max_time = np.max(each_slice_number_max_time[1:])
        f.write("best slice number: {} max depth: {} total time: {} total depth: {}".format(min_slice_number,
                                                                                            best_slice_max_depth,
                                                                                            sample_time + max_max_time,
                                                                                            min_total_depth))
