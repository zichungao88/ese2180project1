import numpy as np
import json

# TODO: Write a system of equations Ax = b to solve for the node voltages

rows = 25
columns = 25
A = np.zeros((rows, columns))
# print(np.matrix(A))
b = np.zeros(columns)
b.shape = (columns, 1)
# print(np.matrix(b))
# print(A.dot(b))


# TODO: Write a function that reads in a file to obtain the resistances on each link
# DONE
def read_resistances_json(file_name):
    resistances = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for entry in data:
            node1 = entry['node1']
            node2 = entry['node2']
            resistance = entry['resistance']
            resistances[(node1, node2)] = resistance
    return resistances

resistance = read_resistances_json('node_resistances.json')
# print(resistance)


# TODO: Write a function that reads in a file to obtain the set of voltages at fixed nodes
# DONE
def read_voltages_json(file_name):
    voltages = {}
    with open(file_name, 'r') as file:
        data = json.load(file)
        for entry in data:
            node = entry['node']
            voltage = entry['voltage']
            voltages[node] = voltage
    return voltages

voltage = read_voltages_json('node_voltages.json')
# print(voltage)


# TODO: Compute the matrix A using the data read in previously


# TODO: Compute the LU factorization of A


# TODO: Compute & output the node voltages & currents through each link


# TODO: Write the output of the previous three steps to a file


# TODO: Repeat but with a tree/graph network