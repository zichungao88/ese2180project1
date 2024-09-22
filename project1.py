import numpy as np
import json

# 1 TODO: Write a system of equations Ax = b to solve for the node voltages
# DONE (More Details of the generic compositions of A & b in written report)
length = 25 # rows = columns
A = np.zeros((length, length))
# print(np.matrix(A))
b = np.zeros(length)
b.shape = (length, 1) # column vector
# print(np.matrix(b))
# print(A.dot(b))


# 2 TODO: Write a function that reads in a file to obtain the resistances on each link
# DONE
def read_resistances_json(file_name):
    resistances = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for value in data:
            node1 = value['node1']
            node2 = value['node2']
            resistance = value['resistance']
            resistances[(node1, node2)] = resistance
    return resistances

resistances = read_resistances_json('node_resistances.json')
# print(resistance)


# 3 TODO: Write a function that reads in a file to obtain the set of voltages at fixed nodes
# DONE
def read_voltages_json(file_name):
    voltages = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for value in data:
            node = value['node']
            voltage = value['voltage']
            voltages[node] = voltage
    return voltages

voltages = read_voltages_json('node_voltages.json')
# print(voltage)


# 4 TODO: Compute the matrix A using the data read in previously
# DONE
def get_neighbors(node_self):
    neighbors = []
    for key in resistances:
        if key[0] == node_self:
            neighbors.append(key[1])
        elif key[1] == node_self:
            neighbors.append(key[0])
    return neighbors

def calculate_A(resistances, voltages, length):
    for node_self in range(1, length + 1):
        if node_self in voltages:
            A[node_self - 1, node_self - 1] = 1
            b[node_self - 1] = voltages[node_self]
        else:
            neighbors = get_neighbors(node_self)
            # print('Node ' + str(node) + ' is connected to node(s) ' + str(neighbor))
            sum_resistances = 0
            for node_neighbor in neighbors:
                if (node_self, node_neighbor) in resistances:
                    resistance = resistances[(node_self, node_neighbor)]
                else:
                    resistance = resistances[(node_neighbor, node_self)]
                # print('The resistance between ' + str(node) + ' and ' + str(neighbor) + ' is ' + str(resistance))
                if resistance > 0:
                    A[node_self - 1, node_neighbor - 1] = -1 / resistance
                    sum_resistances += 1 / resistance
                A[node_self - 1, node_self - 1] = sum_resistances
                b[node_self - 1] = 0
    return A, b

A, b = calculate_A(resistances, voltages, length)
# print(np.matrix(np.round(A, 3)))
# print(np.matrix(np.round(b, 3)))


# 5 TODO: Compute the LU factorization of A
# IN PROGRESS
L = np.zeros((length, length))
U = np.zeros((length, length))

def LU_decomposition(A, length): # needs scrutiny
    for i in range(length - 1):
        for j in range(i + 1, length - 1):
            A[j, i] /= A[i, i]
            for k in range(i + 1, length - 1):
                A[j, k] -= A[j, i] * A[i, k]
    return A

U = LU_decomposition(A, length)
# print(np.matrix(np.round(U, 3)))


# 6 TODO: Compute & output the node voltages & currents through each link
# DONE
x = np.linalg.solve(A, b) # voltage at each node
node_voltages = [0] * length # 25 nodes in 5x5 grid
for i in range(len(node_voltages)):
    node_voltages[i] = float(x[i][0])
    # print('The voltage at Node ' + str(i + 1) + ' is ' + str(node_voltages[i]))

link_currents = [] # 40 resistors connecting all nodes (will append one by one)
for i in resistances:
    v1 = node_voltages[i[0] - 1]
    v2 = node_voltages[i[1] - 1]
    current = (v1 - v2) / resistances.get(i) # define current from 1st node to 2nd node in resistances dictionary key
    link_currents.append(current)
    # print('The current from Node ' + str(i[0]) + ' to ' + str(i[1]) + ' is ' + str(current))


# 7 TODO: Write the output of the previous three steps to a file


# 8 TODO: Repeat but with a tree/graph network