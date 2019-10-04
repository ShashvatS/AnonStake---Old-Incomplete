from collections import namedtuple
import pprint

ASSERT = True

#read the constants file
file = open("constants.txt", "r")

constants = {}

for line in file:
    lst = line.split(' ')
    key, *data = lst

    if data[-1] == "\n":
        data = [int(i) for i in data[:-1]]
    else:
        data = [int(i) for i in data]
    
    constants[key] = data

constants["hashconstant"] = constants["hashconstant"][0]

if ASSERT:
    assert "hashconstant" in constants.keys()
    
    PRFs = ["hash", "addr", "sn", "tsn", "seed", "priority"]
    PRFs += ["binom" + str(i) for i in range(60)]
    
    for key in PRFs:
        assert(len(constants[key]) == 162)

    for i in range(60):
        assert("probvalue" + str(i) in constants.keys())

Circuit = namedtuple('Circuit', ['pub', 'aux', 'constraints'])
circuit = Circuit(pub=["ONE"], aux=set([]), constraints=[])

def LS(a, b):
    return [(a[i], b[i]) for i in range(len(a))]

def CUBE3(circuit, namespace, output, x, y, z):
    tmp = namespace + ".tmp"
    
    s1 = LS([1, 1, z], [x, y, "ONE"])
    s2 = LS([1], [tmp])
    
    c1 = (s1, s1, s2)
    c2 = (s1, s2, LS([1], [output]))

    circuit.aux.update([tmp])
    circuit.constraints.extend([c1, c2])

    return

def MiMC(circuit, namespace, output, x, k, mimc_constants):
    tmp1 = namespace + ".tmp1"
    tmp2 = namespace + ".tmp2"
    x1 = namespace + ".x1"
    x161 = namespace + ".x161"

    s1 = LS([1, 1], [x, k])
    s2 = LS([1], [tmp1])

    c1 = (s1, s1, s2)
    c2 = (s1, s2, LS([1], [x1]))

    s3 = LS([1, 1, mimc_constants[161]], [x161, k, "ONE"])
    s4 = LS([1], [tmp2])

    s5 = LS([1, -1], [output, k])

    c3 = (s3, s3, s4)
    c4 = (s3, s4, s5)

    for i in range(1, 161):
        newnamespace = namespace + ".cubespace" + str(i)
        xold = namespace + ".x" + str(i)
        xnew = namespace + ".x" + str(i + 1)

        circuit.aux.update([xnew, xold])

        CUBE3(circuit, newnamespace, xnew, xold, k, mimc_constants[i])

    circuit.aux.update([tmp1, tmp2, x1, x161])
    circuit.constraints.extend([c1, c2, c3, c4])

    return

def CRH(circuit, namespace, output, m1, m2):
    h1 = namespace + ".h1"

    circuit.aux.update([h1])
    
    hashconstant = constants["hashconstant"]
    mimc_constants = constants["hash"]
    
    MiMC(circuit, "mimcspace1", h1, hashconstant, m1, mimc_constants)
    a, b, c = circuit.constraints[-1]
    c.append([-1, hashconstant])
    circuit.constraints[-1] = (a, b, c)

    MiMC(circuit, "mimcspace2", output, h1, m2, mimc_constants)

    a, b, c = circuit.constraints[-1]
    c.append([-1, h1])
    circuit.constraints[-1] = (a, b, c)

def MerkleProof(circuit, namespace, root, node, path):
    pass


def eliminateConstantMultiplications(circuit):
    #replace constraints that have multiplication of two constants


    pass

CRH(circuit, "", "X", "m1", "m2")

pp = pprint.PrettyPrinter()
pp.pprint(circuit.constraints)