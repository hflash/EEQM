# -*- coding: UTF-8 -*-
import os
import copy
import random
import math
import importlib
import networkx as nx
import matplotlib.pyplot as plt
from quantumcircuit.gate import *
from quantumcircuit.register import *


class QuantumCircuit:
    def __init__(self, qubit_number=None, cbit_number=0):
        if qubit_number ==None:
            self.qubit_number = 0
        else:
            self.qubit_number = qubit_number
        self.qubit_registers = [None for _ in range(qubit_number)]
        self.cbit_number = cbit_number
        self.cbit_registers = [None for _ in range(cbit_number)]
        self.gate_number = 0
        self.gate_set = set()
        self.last_gate_on_qubit = [None for _ in range(qubit_number)]
        self.begin_gate_on_qubit = [None for _ in range(qubit_number)]
        self.graph_circuit = nx.DiGraph()
        self.gate_list = []
        self.initialize(Init(qubit_number))

    def add_gate(self, gate):
        qubits = gate.get_qubits()
        node_id = self.gate_number
        self.gate_number = self.gate_number + 1
        self.graph_circuit.add_node(node_id, gate=gate)
        for qubit in qubits:
            assert qubit < self.qubit_number
            if self.begin_gate_on_qubit[qubit] == None:
                self.last_gate_on_qubit[qubit] = node_id
                self.begin_gate_on_qubit[qubit] = node_id
            else:
                self.graph_circuit.add_edge(self.last_gate_on_qubit[qubit], node_id)
                self.last_gate_on_qubit[qubit] = node_id
        self.gate_set.add(gate.get_name())
        self.gate_list.append(gate)

    def get_gate_number(self):
        return self.gate_number

    def get_qubit_number(self):
        return self.qubit_number

    def get_gate_set(self):
        return self.gate_set

    def get_graph_list(self):
        return self.gate_list

    def get_circuit_depth(self):
        return len(nx.dag_longest_path(self.graph_circuit))

    def to_dag_graph(self):
        return self.graph_circuit

    def print_dag_circuit(self):
        dag = self.graph_circuit
        pos = nx.spring_layout(dag)
        # pos = nx.bipartite_layout(dag, nodes=dag.nodes)
        # node_colors = ['lightgray' for _ in node_names]
        nx.draw_networkx(dag, pos, with_labels=True)
        nx.draw_networkx(dag, pos, with_labels=True, node_color='lightgray', font_size=14)
        labels = nx.get_node_attributes(dag, 'name')
        nx.draw_networkx_labels(dag, pos, labels)
        plt.show()

    def get_qubit_duration_path(self, qubit):
        '''
        获取这个逻辑比特需要持续的时间：dag中，这个比特上执行的第一个门，到最后一个门，中间的所有路径
        only for transpiled circuit
        :param qubit:
        :return: paths of gates depend on this qubit
        '''
        begin_node = self.begin_gate_on_qubit[qubit]
        end_node = self.last_gate_on_qubit[qubit]
        if begin_node == end_node:
            path_nodes = []
            path_node = []
            path_node.append(self.gate_list[begin_node])
            path_nodes.append(path_node)
            return path_nodes
        else:
            paths = list(nx.all_simple_paths(self.graph_circuit, source=begin_node, target=end_node))
            path_nodes = []
            for path in paths:
                path_node = []
                for i in path:
                    path_node.append(self.gate_list[i])
                path_nodes.append(path_node)
        return path_nodes

    def initialize(self, gate):
        '''
        量子线路的初始化，后续会接入quantum state中量子态的生成等步骤
        :param gate:
        :return:
        '''
        node_id = self.gate_number
        self.gate_number = self.gate_number + 1
        self.graph_circuit.add_node(node_id, gate=gate)
        qubits = gate.get_qubits()
        for qubit in qubits:
            self.last_gate_on_qubit[qubit] = node_id
            self.begin_gate_on_qubit[qubit] = node_id
        self.gate_set.add(gate.get_name())
        self.gate_list.append(gate)

    def get_gate(self, gate_id):
        return self.gate_list[gate_id]

    def to_QASM(self, filename):
        try:
            with open(filename, "w") as f:
                f.write("OPENQASM 2.0;\n")
                f.write("include \"qelib1.inc\";\n")
                f.write("qreg q[" + str(self.qubit_number) + "];\n")
                f.write("creg c[" + str(self.cbit_number) + "];\n")
                for i in range(len(self.gate_list) - 1):
                    gate = self.gate_list[i + 1]
                    qubits = gate.get_qubits()
                    gate_name = gate.get_name()
                    gate_name_l = gate_name.lower()
                    if gate_name == "X" or gate_name == "Y" or gate_name == "Z" or gate_name == "H" or gate_name == "S" or gate_name == "T" or gate_name == "TDG":
                        qasm_s = gate_name_l + " q[" + str(qubits[0]) + "];\n"
                        f.write(qasm_s)
                    elif gate_name == "RX" or gate_name == "RZ":
                        theta = gate.get_para()
                        qasm_s = gate_name_l + "(" + str(theta) + ") q[" + str(qubits[0]) + "];\n"
                        f.write(qasm_s)
                    elif gate_name == "CX":
                        qasm_s = "cx q[" + str(qubits[0]) + "], q[" + str(qubits[1]) + "];\n"
                        f.write(qasm_s)
                    else:
                        print(gate_name, " is not supported")
        except Exception as e:
            print("Error in creating qasm file:", e)

    def get_circuit_gate_info(self):
        '''
        获取量子线路的信息，知道每个门有多少个
        Returns:{gate_name:gate_num} gate_name门的名字，gate_num这个门在线路总出现的次数
        '''
        gate_name_num = {}
        for i in range(len(self.gate_list)):
            gate = self.gate_list[i]
            gate_name = gate.get_name()
            if gate_name in gate_name_num.keys():
                gate_name_num[gate_name] = gate_name_num[gate_name]+1
            else:
                gate_name_num[gate_name] = 1
        return gate_name_num



    def to_QASM_list(self):
        '''
        转换成qasm list 一种输出到模拟器的格式
        :return:
        '''
        qasm_list = []
        qasm_list.append("OPENQASM 2.0;")
        qasm_list.append("include \"qelib1.inc\";")
        qasm_list.append("qreg q[" + str(self.qubit_number) + "];")
        qasm_list.append("creg c[" + str(self.cbit_number) + "];")
        for i in range(len(self.gate_list) - 1):
            gate = self.gate_list[i + 1]
            qubits = gate.get_qubits()
            gate_name = gate.get_name()
            gate_name_l = gate_name.lower()
            if gate_name == "X" or gate_name == "Y" or gate_name == "Z" or gate_name == "H" or gate_name == "S" or gate_name == "T" or gate_name == "TDG":
                qasm_s = gate_name_l + " q[" + str(qubits[0]) + "];"
                qasm_list.append(qasm_s)
            elif gate_name == "RX" or gate_name == "RZ":
                theta = gate.get_para()
                qasm_s = gate_name_l + "(" + str(theta) + ") q[" + str(qubits[0]) + "];"
                qasm_list.append(qasm_s)
            elif gate_name == "CX":
                qasm_s = "cx q[" + str(qubits[0]) + "], q[" + str(qubits[1]) + "];"
                qasm_list.append(qasm_s)
            else:
                print(gate_name, " is not supported")
        return qasm_list

    def add_QuantumCircuit(self, newQC):
        '''
        1.新的量子线路的量子比特数目小于等于老的量子线路的比特数目
        2.默认新的量子线路的比特i就是老量子线路的比特i
        :param newQC:
        :return:
        '''
        if self.qubit_number >= newQC.qubit_number:
            new_gate_list = newQC.get_graph_list()
            for index in range(len(new_gate_list)):
                new_gate = copy(new_gate_list[index])
                if new_gate.name != "Init":
                    self.add_gate(new_gate)
        else:
            raise ValueError(
                "The number of qubits in the new quantum circuit needs to be less than or equal to the number of "
                "qubits in the original circuit.")


    def add_QuantumCircuit(self, newQC,mapping):
        '''
        1.新的量子线路的量子比特数目小于等于老的量子线路的比特数目
        2.mapping[i]是指新线路上的逻辑比特i应该是老线路的哪个逻辑比特
        :param newQC:
        :return:
        '''
        if self.qubit_number >= newQC.qubit_number:
            new_gate_list = newQC.get_graph_list()
            for index in range(len(new_gate_list)):
                new_gate = copy(new_gate_list[index])
                if new_gate.name != "Init":
                    if new_gate.name == "CX":
                        control_qubit, target_qubit = new_gate.get_qubits()
                        new_control = mapping[control_qubit]
                        new_target = mapping[target_qubit]
                        new_gate.control_qubit = new_control
                        new_gate.target_qubit = new_target
                        self.add_gate(new_gate)
                    else:
                        old_qubit = new_gate.get_qubits()
                        new_qubit = mapping[old_qubit[0]]
                        new_gate.qubit = new_qubit
        else:
            raise ValueError(
                "The number of qubits in the new quantum circuit needs to be less than or equal to the number of "
                "qubits in the original circuit.")

    @classmethod
    def from_QASM(cls, filename):
        try:
            with open(filename, "r") as file:
                f = file.readlines()
                for i in range(len(f)):
                    if "OPENQASM" in f[i] or "include \"qelib1.inc\"" in f[i] or "creg" in f[i]:
                        continue
                    elif "qreg" in f[i]:
                        numbers = re.findall(r'\d+', f[i])
                        qubit_num = int(numbers[0])
                        qc = QuantumCircuit(qubit_num)
                    elif f[i] == '\n' or '//' in f[i]:
                        continue
                    else:
                        pattern = r'^[a-zA-Z]+'
                        match = re.match(pattern, f[i])
                        gate_name = match[0].upper()
                        if gate_name == "X" or gate_name == "Y" or gate_name == "Z" or gate_name == "H" or gate_name == "S" or gate_name == "T" or gate_name == "TDG":
                            numbers = re.findall(r'\d+', f[i])
                            qubit_num = int(numbers[0])
                            module = importlib.import_module("quantumcircuit.gate")
                            function = getattr(module, gate_name)
                            newgate = function(qubit_num)
                            qc.add_gate(newgate)
                        elif gate_name == "RX" or gate_name == "RZ":
                            pattern = r"(-?\d+\.?\d*)"
                            matches = re.findall(pattern, f[i])
                            qubit = int(matches[1])
                            theta = float(matches[0])
                            module = importlib.import_module("quantumcircuit.gate")
                            function = getattr(module, gate_name)
                            newgate = function(qubit, theta)
                            qc.add_gate(newgate)
                        elif gate_name == "CX":
                            numbers = re.findall(r'\d+', f[i])
                            q1 = int(numbers[0])
                            q2 = int(numbers[1])
                            newgate = CX(q1, q2)
                            qc.add_gate(newgate)
                        elif gate_name == "U":
                            pattern = r"(-?\d+\.?\d*)"
                            matches = re.findall(pattern, f[i])
                            qubit = int(matches[-1])
                            # theta = float(matches[0])
                            module = importlib.import_module("quantumcircuit.gate")
                            function = getattr(module, 'X')
                            newgate = function(qubit)
                            qc.add_gate(newgate)
                        else:
                            print("not supported gate")

                return qc
        except Exception as e:
            print("Error in reading qasm file:", e)

    @classmethod
    def load_alg(cls, filename):
        '''
        从已有的算法库里加载已有的量子线路
        :param filename:
        :return:
        '''
        files = []
        list_directory_contents("../qasm-benchmark", files)
        ready_file = []
        for i in files:
            if compare_strings(i, filename):
                ready_file.append(i)
        print(ready_file)
        data = input("Input filename: ")
        # data = "3_17_13"
        for i in ready_file:
            if compare_strings(i, data):
                filename = i
        with open(filename, "r") as file:
            f = file.readlines()
            qubit_number = 0
            for i in range(len(f)):
                if "qreg" in f[i]:
                    numbers = re.findall(r'\d+', f[i])
                    qubit_number = int(numbers[0])
            qc = QuantumCircuit(qubit_number)
            for i in range(len(f)):
                if "OPENQASM" in f[i] or "include \"qelib1.inc\"" in f[i] or "creg" in f[i] or "\n" == f[i]:
                    continue
                elif "qreg" in f[i]:
                    pass
                    # numbers = re.findall(r'\d+', f[i])
                    # qubit_num=int(numbers[0])
                    # qc=QuantumCircuit(qubit_num)
                else:
                    pattern = r'^[a-zA-Z]+'
                    match = re.match(pattern, f[i])
                    gate_name = match[0].upper()
                    if gate_name == "X" or gate_name == "Y" or gate_name == "Z" or gate_name == "H" or gate_name == "S" or gate_name == "T" or gate_name == "TDG":
                        numbers = re.findall(r'\d+', f[i])
                        qubit_num = int(numbers[0])
                        module = importlib.import_module("quantumcircuit.gate")
                        function = getattr(module, gate_name)
                        newgate = function(qubit_num)
                        qc.add_gate(newgate)
                    elif gate_name == "RX" or gate_name == "RZ":
                        pattern = r"(-?\d+\.?\d*)"
                        matches = re.findall(pattern, f[i])
                        qubit = int(matches[1])
                        theta = float(matches[0])
                        module = importlib.import_module("quantumcircuit.gate")
                        function = getattr(module, gate_name)
                        newgate = function(qubit, theta)
                        qc.add_gate(newgate)
                    elif gate_name == "CX":
                        numbers = re.findall(r'\d+', f[i])
                        q1 = int(numbers[0])
                        q2 = int(numbers[1])
                        newgate = CX(q1, q2)
                        qc.add_gate(newgate)
                    else:
                        print("not supported gate")
            return qc

    @classmethod
    def random_circuit(cls, num_qubits, depth, gate1p, gate2p, gate_set=None):
        '''
        generate a random circuit
        :param num_qubits: qubits number of circuits
        :param depth: depth of circuits
        :param gate1p: probability of single qubit gate
        :param gate2p: probability of two qubit gate
        :param gate_set: used gate set
        :return: quantum circuit
        '''

        if num_qubits < 0 or depth < 0:
            print("wrong qubit number or depth")
            return
        if gate1p < 0 or gate2p < 0:
            print("Probability less than 0")
            return
        if gate1p + gate2p * 2 > 1:
            print("Probability greater than 1")
            return
        qc = QuantumCircuit(num_qubits)
        if num_qubits == 0:
            return qc
        gate1list = ["X", "Y", "Z", "H", "S", "T", "RX", "RZ"]
        gate2list = ["CX"]
        paragateset = {"RX", "RZ"}

        # generate gate on critical path
        total_p = gate1p + gate2p * 2
        critical_gate_list = []
        used_qubit = []
        qubits = list(range(num_qubits))
        for _ in range(depth):
            random_number = random_float(total_p)
            # two qubit gate
            if random_number > gate1p:
                # first layer
                if len(used_qubit) == 0:
                    selected_elements = random.sample(qubits, 2)
                    new_gate = CX(selected_elements[0], selected_elements[1])
                    used_qubit.append(selected_elements[0])
                    used_qubit.append(selected_elements[1])
                    critical_gate_list.append(new_gate)
                else:
                    q1 = random.choice(used_qubit)
                    new_qubit_list = list(range(num_qubits))
                    new_qubit_list.remove(q1)
                    q2 = random.choice(new_qubit_list)
                    new_gate = CX(q1, q2)
                    used_qubit.append(q1)
                    used_qubit.append(q2)
                    critical_gate_list.append(new_gate)
            else:
                # single qubit gate
                # which qubit
                if len(used_qubit) == 0:
                    q1 = random.choice(qubits)
                else:
                    q1 = random.choice(used_qubit)
                # which gate
                gate_name = random.choice(gate1list)
                if gate_name in paragateset:
                    parameter = random.uniform(0, 2 * math.pi)
                    module = importlib.import_module("quantumcircuit.gate")
                    function = getattr(module, gate_name)
                    newgate = function(q1, parameter)
                    critical_gate_list.append(newgate)
                    used_qubit = [q1]
                else:
                    module = importlib.import_module("quantumcircuit.gate")
                    function = getattr(module, gate_name)
                    newgate = function(q1)
                    critical_gate_list.append(newgate)
                    used_qubit = [q1]

        # generate quantum circuit
        for index in range(depth):
            single_qubit_gate_num = 0
            two_qubit_gate_num = 0
            qubits = list(range(num_qubits))
            qc.add_gate(critical_gate_list[index])
            if critical_gate_list[0].name == "CX":
                control_target = critical_gate_list[0].get_qubits()
                qubits.remove(control_target[0])
                qubits.remove(control_target[1])
                two_qubit_gate_num = two_qubit_gate_num + 1
            else:
                used_qubit = critical_gate_list[0].get_qubits()
                qubits.remove(used_qubit[0])
                single_qubit_gate_num = single_qubit_gate_num + 1
            # add single qubit gate
            while ((single_qubit_gate_num + 1) / num_qubits < gate1p):
                q1 = random.choice(qubits)
                module = importlib.import_module("quantumcircuit.gate")
                function = getattr(module, gate_name)
                if gate_name in paragateset:
                    parameter = random.uniform(0, 2 * math.pi)
                    newgate = function(q1, parameter)
                else:
                    newgate = function(q1)
                qc.add_gate(newgate)
                qubits.remove(q1)
                single_qubit_gate_num = single_qubit_gate_num + 1
            # add two qubit gate
            while ((two_qubit_gate_num + 1) / num_qubits < gate2p):
                [q1, q2] = random.sample(qubits, 2)
                new_gate = CX(q1, q2)
                qc.add_gate(newgate)
                qubits.remove(q1)
                qubits.remove(q2)
                two_qubit_gate_num = two_qubit_gate_num + 1
        return qc

    def to_dagtable(self):
        '''
        dagtable: 2D list
        Each column of dagtable is a logical qubit,
        and each row is a time slice.
        If it is -1, it means that this time period is free,
        otherwise, it means the id of the gate to be executed.
        :return:
        '''
        # generate dagtable
        dagtable = [[] for _ in range(self.qubit_number)]
        # add gate from gate list
        for i in range(len(self.gate_list)):
            gate = self.gate_list[i]
            if gate.name == "Init":
                continue
            qubits = gate.get_qubits()
            if len(qubits) == 1:
                dagtable[qubits[0]].append(i)
            else:
                while (len(dagtable[qubits[0]]) < len(dagtable[qubits[1]])):
                    dagtable[qubits[0]].append(-1)
                while (len(dagtable[qubits[0]]) > len(dagtable[qubits[1]])):
                    dagtable[qubits[1]].append(-1)
                dagtable[qubits[0]].append(i)
                dagtable[qubits[1]].append(i)
        # complete
        max_length = max(len(row) for row in dagtable)
        for i in range(len(dagtable)):
            while len(dagtable[i]) < max_length:
                dagtable[i].append(-1)
        return dagtable

    def front_layer(self):
        dagtable = self.to_dagtable()
        front_set = set()
        for i in range(self.qubit_number):
            if dagtable[i][0] != -1:
                front_set.add(dagtable[i][0])
        gate_ids = list(front_set)
        return gate_ids

    def remove_gate(self, gate_id):
        pass


def random_float(p):
    return random.random() * p


def get_first_word(string):
    words = string.split()
    if len(words) > 0:
        return words[0]
    else:
        return ""


def compare_strings(a, b):
    a_lower = a.lower()
    b_lower = b.lower()
    if b_lower in a_lower:
        return True
    else:
        return False


def list_directory_contents(path, files):
    for child in os.listdir(path):
        child_path = os.path.join(path, child)
        if os.path.isfile(child_path):
            files.append(child_path)
        elif os.path.isdir(child_path):
            list_directory_contents(child_path, files)
