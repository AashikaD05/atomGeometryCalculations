from math import cos, acos, atan, sin, sqrt, radians, degrees, atan2
import numpy as np
import string
import sys
import re

#Tolerance
TOL = 0.05**(-10)

class blc:
    def __init__(self, coordFilename, atomRadFilename):
        self.coordFilename = coordFilename
        self.atomRadFilename = atomRadFilename
        
    #Return fileobject
    def openFile(self, filename):
        return open(filename, 'r')
    
    #Returning index of search string
    def indexSearch(self, s, file):
        for r in file:
            if r.startswith(s):
                return file.index(r)

    #Returning row of given index
    def rowReturn(self, a, file, atoms):
        return file[atoms[a]].split()
    
    #Checking if a string is a float or not
    def is_float(self, s):
        if '.' in s:
            return True
        else:
            return False
        
    #Returning all atoms present in the file
    def atomList(self, file):
        atoms = {}
        for idx, r in enumerate(file):
            rList = r.split()
            if len(rList) >=5:
                if (rList[0].isalnum() 
                        and rList[1].isdigit() 
                        and self.is_float(rList[2]) 
                        and self.is_float(rList[3]) 
                        and self.is_float(rList[4])) or (rList[-1] == "="):
                    atoms[rList[0]] = idx

        atomval = list(atoms.values())
        keyToDelete = []
        i = 0

        #Adding values to a new list which are to be deleted
        while i<len(atomval)-1:
            if (abs(atomval[i] - atomval[i+1]) not in [1, 2]):
                keyToDelete.append(atomval[i])
            i+=1

        #Deleting non atom lines
        for k in list(atoms.keys()):
            if atoms[k] in keyToDelete:
                del atoms[k]
        return atoms    

    #Computing G or M matrix
    def matrixCalculation(self, coordFile, Type):
        dimensions = coordFile[self.indexSearch("CELL", coordFile)].split()

        dimensions.pop(0)
        dimensions.pop(0)
        dimensions = list(map(float,dimensions))

        #Obtaining the matrix
        #               0  1  2    3     4      5
        # dimensions = [a, b, c, alpha, beta, gamma]

        a, b, c = dimensions[:3]
        alpha, beta, gamma = map(radians, dimensions[3:6])

        if Type == "G":
            G = [
                    [a**2, a*b*cos(gamma), a*c*cos(beta)],
                    [0, b**2, b*c*cos(alpha)],
                    [0, 0, c**2]
            ]

            return G
        else:
            V = c * sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*cos(alpha)*cos(beta)*cos(gamma)) / sin(gamma)
            M = [
                    [a, b * cos(gamma), c * cos(beta)],
                    [0, b * sin(gamma), c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)],
                    [0, 0, V]
            ]
            return M

    #Fractional to Cartesian Coordinates conversion
    def f_to_c(self, coordinates, coordfile, n):
        #Computing Transformation Matrix
        M = self.matrixCalculation(coordFile, "M")

        cartesian = [[] for i in range(n)]
        for i in range(n):
            for j in range(3):
                sum, l = 0, 2
                for k in range(3):
                    sum += float(coordinates[i][l]) * M[j][k]
                    l += 1
                cartesian[i].append(sum) 
        return cartesian

    #Bond Length Calculation
    def BLCalc(self, a1, a2, atoms):

        #Opening file name
        coordFile = self.openFile(self.coordFilename).readlines()
        atomRadFile = self.openFile(self.atomRadFilename).readlines()

        #Obtaining the matrix
        G = self.matrixCalculation(coordFile, "G")

        #Obtaining the coordinates from the file for the given atoms
        a1_row = self.rowReturn(a1, coordFile, atoms)
        a2_row = self.rowReturn(a2, coordFile, atoms)

        #Calculating Positional matrix
        v = [float(a1_row[i]) - float(a2_row[i]) for i in range(2,5)]

        bl_squared = 0
        for i in range(3):
            for j in range(3):
                bl_squared += v[j] * G[j][i] * v[i]
        bond_length = sqrt(bl_squared)

        #Searching for Atom radius 
        rad_a1, rad_a2 = -1,-1

        #Extracting the string part from the input
        a1_s = "".join(re.findall(r'[A-Za-z]',a1))
        a2_s = "".join(re.findall(r'[A-Za-z]',a2))

        #Determining Atom's Radius from the file
        for r in atomRadFile:
            r = r.strip().split()
            if r[0].upper() == a1_s:
                rad_a1 = float(r[2])
            if r[0].upper() == a2_s:
                rad_a2 = float(r[2])
            else:
                continue

        if rad_a1 == -1 or rad_a2 == -1 or rad_a1 + rad_a2 + TOL < bond_length:
            return False
        else:
            return bond_length            

    #Bond Angle Calculation
    def BACalc(self, a1, a2, a3, atoms):
        #Opening File
        coordFile = self.openFile(self.coordFilename).readlines()

        #Obtaining the coordinates from the file for the given atoms
        a1_row = self.rowReturn(a1, coordFile, atoms)
        a2_row = self.rowReturn(a2, coordFile, atoms)
        a3_row = self.rowReturn(a3, coordFile, atoms)

        #Obtaining the matrix
        G = self.matrixCalculation(coordFile, "G")

        #Computing vectors
        u = [float(a1_row[i]) - float(a2_row[i]) for i in range(2,5)]
        v = [float(a3_row[i]) - float(a2_row[i]) for i in range(2,5)]

        #Magnitude of vectors
        u_mag = sqrt(
                G[0][0]*u[0]**2 
            +   G[1][1]*u[1]**2 
            +   G[2][2]*u[2]**2 
            +   2*G[0][1]*u[0]*u[1]
            +   2*G[0][2]*u[0]*u[2]
            +   2*G[1][2]*u[1]*u[2]
        )
        v_mag = sqrt(
                G[0][0]*v[0]**2 
            +   G[1][1]*v[1]**2 
            +   G[2][2]*v[2]**2 
            +   2*G[0][1]*v[0]*v[1]
            +   2*G[0][2]*v[0]*v[2]
            +   2*G[1][2]*v[1]*v[2]
        )

        dot_product = 0
        for i in range(3):
            for j in range(3):
                dot_product += u[i] * G[i][j] * v[j]
        
        theta = degrees(acos(dot_product/(u_mag*v_mag)))
        return theta
    
    #Torsional Angle Calculation
    def TACalc(self, a1, a2, a3, a4, atoms):
        coordFile = self.openFile(self.coordFilename).readlines()

        #Obtaining coordinates for each atom
        coordinates = [self.rowReturn(a1, coordFile, atoms),    
                       self.rowReturn(a2, coordFile, atoms),    
                       self.rowReturn(a3, coordFile, atoms),    
                       self.rowReturn(a4, coordFile, atoms)]    

        #Converting Fractional to Cartesian Coordinates
        cartesian = self.f_to_c(coordinates, coordFile, 4)

        #Computing Vectors
        v1, v2, v3 = [], [], []

        for i in range(3):
            v1.append(float(cartesian[1][i]) - float(cartesian[0][i]))
            v2.append(float(cartesian[2][i]) - float(cartesian[1][i]))
            v3.append(float(cartesian[3][i]) - float(cartesian[2][i]))

        #Cross Product of Vectors  
        n1 = [ v1[1]*v2[2] - v1[2]*v2[1],
               v1[2]*v2[0] - v1[0]*v2[2],
               v1[0]*v2[1] - v1[1]*v2[0]]        

        n2 = [ v2[1]*v3[2] - v2[2]*v3[1],
               v2[2]*v3[0] - v2[0]*v3[2],
               v2[0]*v3[1] - v2[1]*v3[0]]  

        #Computing Torsional Angle 
        v2_mag = sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
        n1_v3 = n1[0]*v3[0] + n1[1]*v3[1] + n1[2]*v3[2]
        n1_n2 = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]

        theta = degrees(atan2(v2_mag*n1_v3, n1_n2))
        return theta
    

    #Least Square Lines Calculation
    def leastSquareCalc(self, inputted_atoms, n, atoms):
        coordFile = self.openFile(self.coordFilename).readlines()

        #Obtaining coordinates of inputted atoms only
        coordinates = []
        for i in inputted_atoms:
            coordinates.append(self.rowReturn(i, coordFile, atoms))

        #Converting Fractional to cartesian coordinates
        cartesian = self.f_to_c(coordinates, coordFile, n)

        #Calculating centroid 
        p_x, p_y, p_z = 0, 0, 0
        for r in cartesian:
            p_x += float(r[0])
            p_y += float(r[1])
            p_z += float(r[2])
        p_x, p_y, p_z = round(p_x/n, 2), round(p_y/n, 2), round(p_z/n, 2)

        #Building centered data matrix
        X = [[] for i in range(n)]
        for i in range(n):
            X[i].append(cartesian[i][0] - p_x)
            X[i].append(cartesian[i][1] - p_y)
            X[i].append(cartesian[i][2] - p_z)
        
        #Computing Covariance matrix 
        C = [[] for i in range(3)]
        for i in range(3):
            for j in range(3):
                sum = 0
                for k in range(n):
                    sum += X[k][i] * X[k][j]
                C[i].append(sum)  
        C = np.array(C)

        #Determining max eigenvalue and it's corresponding eigenvector
        eigenvalues, eigenvectors = np.linalg.eig(C)
        max_id = np.argmax(eigenvalues)
        dir_vector = eigenvectors[:, max_id]

        r_t = f"({p_x}, {p_y}, {p_z}) + t.({np.round(dir_vector[0], 2)}, {np.round(dir_vector[1], 2)}, {np.round(dir_vector[2], 2)})"
        return r_t

#main function 
if __name__ == "__main__":

    coordFilename = r"E:\Desktop\pdfs\IISc\Dinesh sir\Bond Length Calculator\sigi1.txt"
    atomRadFile = r"E:\Desktop\pdfs\IISc\Dinesh sir\Bond Length Calculator\atomcols.def"

    #Creating an object
    blcObj = blc(coordFilename, atomRadFile)

    ans = 'Y'
    while ans == 'Y' or ans == 'y':

        coordFile = (blcObj.openFile(coordFilename)).readlines()
        atoms = blcObj.atomList(coordFile)

        print("\n1. Bond Length\n" \
        "2. Bond Angle\n" \
        "3. Torsional Angle (dihedral)\n" \
        "4. Least Square Line Calculation\n")
        
        ch = int(input("Enter choice number: "))
        print(f"\nAtoms List present in this file :\n{list(atoms.keys())}\n")
        
        if ch == 1:
            a1 = input("Enter first atom: ").upper()
            a2 = input("Enter second atom: ").upper()

            bond_length = round(blcObj.BLCalc(a1, a2, atoms), 4)

            if bond_length == False:
                print("Either atom radius is not given or Bond length is too small..")
            else:
                print(f"Bond Length between {a1} and {a2} is {bond_length} Å\n")

        elif ch == 2:
            a1, a2, a3 = input("Enter any 3 atoms [eg. AS AG C]: ").upper().split()
            bond_angle = round(blcObj.BACalc(a1, a2, a3, atoms), 2)
            print(f"Bond Angle between {a1} {a2} {a3} is {bond_angle} ∘")

        elif ch == 3:
            a1, a2, a3, a4 = input("Enter any 4 atoms [eg. AS AG C N]: ").upper().split()
            torsional_angle = round(blcObj.TACalc(a1, a2, a3, a4, atoms), 2)
            print(f"Torsional angle between planes {a1}-{a2}-{a3} and {a2}-{a3}-{a4} is {torsional_angle} ∘")

        elif ch == 4:
            n = int(input("Enter number of atoms: "))
            inputted_atoms = []
            for i in range(n):
                a = input(f"Enter {i+1} atom: ").upper()
                inputted_atoms.append(a)
            r_t = blcObj.leastSquareCalc(inputted_atoms, n, atoms)
            print(f"Least Squares line through {inputted_atoms} is {r_t}")

        else:
            print("Incorrect choice..aborting")
            sys.exit(1)