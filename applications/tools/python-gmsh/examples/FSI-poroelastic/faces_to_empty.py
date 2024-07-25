
import os
import sys

directory = os.path.dirname(os.path.abspath(__file__))
filename = sys.argv[1]
filepath = os.path.join(directory, filename, "constant/polyMesh/boundary")
boundary = open(filepath)
boundary_list = boundary.readlines()
ind1 = boundary_list.index("    frontAndBack\n")
boundary_list[ind1+2] = "        type            empty;\n"
boundary_list[ind1+3] = "        physicalType    empty;\n"
boundary.close()
with open(filepath, 'w') as f:
    f.writelines(boundary_list)
