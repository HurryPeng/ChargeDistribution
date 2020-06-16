# test.py
# HurryPeng
# 2020.3.7

import numpy as np

def calcPointChargeField(coordCharge, q, coord):
    # np.array calcPointChargeField(np.array coordCharge, float q, np.array coord)
    K = 9.0E9
    deltaR = coord - coordCharge
    dist = np.linalg.norm(deltaR)
    if dist == 0.0:
        return np.array([0.0, 0.0])
    return K * q / dist**3 * deltaR

def calcPointChargeForce(coordSource, qSource, coord, q):
    # np.array calcPointChageForce(np.array coordSource, float qSource, np.array coord, float q)
    return q * calcPointChargeField(coordSource, qSource, coord)

class ElectricField(object):
    def __init__(self, _strength):
        # void ElectricField::__init__(np.array _strenth(np.array))
        self.strength = _strength
            # np.array ElectricField::strenth(np.array)

    def __add__(self, rhs):
        return ElectricField(lambda coord: self.strength(coord) + rhs.strength(coord))

    @classmethod
    def getUniformField(self, uniformField):
        # void ElectricField::getUniformField(np.array uniformField)
        return ElectricField(lambda coord: uniformField)

    @classmethod
    def getPointChargeField(self, coordCharge, q):
        # void ElectricField::getPointChargeField(np.array coordCharge, float q)
        return ElectricField(lambda coord: calcPointChargeField(coordCharge, q, coord))

print("########################\n")
field = ElectricField.getPointChargeField(np.array([1.50, 0.0]), 2e-10)
field += ElectricField.getPointChargeField(np.array([0.0, 1.50]), -1e-10)
field += ElectricField.getPointChargeField(np.array([-2.50, 0.0]), -1.5e-10)
field += ElectricField.getUniformField(np.array([0.02, 0.0]))

file = open("out.txt", "w")

for x in range(-50, 51, 1):
    for y in range(-50, 51, 1):
        vec = field.strength(np.array([float(x/10), float(y/10)]))
        print("{{%f, %f}, {%f, %f}}," %(x/10, y/10, vec[0], vec[1]), file = file)
            

file.close()

print(field.strength(np.array([0.0, 0.0])))
