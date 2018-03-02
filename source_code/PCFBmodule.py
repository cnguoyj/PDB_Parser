import math as m

#This function parses the given PDB file and returns a list with the coordinates
#of all atoms, a list with the name of all the residues and a list with all the
#lines that refer to HETATM.
def parsePDB(file_path):
    file = open(file_path)
    coordinates = []
    residues = []
    heterogens = []
    for line in file:
        if line.startswith("HETATM"):
            heterogens.append(line)
        line = line.split()
        temp_coordinates = []
        if line[0] == "ATOM":
            if line[2] == "CA":
                if not line[3].startswith("B"):
                    temp_coordinates.append(float(line[6]))
                    temp_coordinates.append(float(line[7]))
                    temp_coordinates.append(float(line[8]))
                    coordinates.append(temp_coordinates)
                    residues.append(line[3])
    return coordinates, residues, heterogens
        

#This function iterates over a given list with coordinates from atoms and calculates
#the centroid and the distance of all atoms to the centroid returning a list
#with all those distances.
def distance2Center(coordinates_list):
    x = 0
    y = 0
    z = 0
    counter = 0
    for coordinates in coordinates_list:
        x += coordinates[0]
        y += coordinates[1]
        z += coordinates[2]
        counter += 1
    x = x / counter
    y = y / counter
    z = z / counter
    distance2Center = []
    for coordinates in coordinates_list:
        x_temp = coordinates[0]
        y_temp = coordinates[1]
        z_temp = coordinates[2]
        distance = m.sqrt((x - x_temp)**2 + (y - y_temp)**2 + (z - z_temp)**2)
        distance2Center.append(distance)
    return distance2Center
    
#This function iterates over a given list with all HETATM lines from a PDB file
#and returns the coordinates for the given atom in the form of a list.    
def searchHeterogen(heterogens_list, atom):
    heterogens_coordinates = []
    for i in range(len(heterogens_list)):
        line = heterogens_list[i].split()
        temp_coordinates = []
        if line[2] == atom:
            temp_coordinates.append(float(line[6]))
            temp_coordinates.append(float(line[7]))
            temp_coordinates.append(float(line[8]))
            heterogens_coordinates.append(temp_coordinates)
    return heterogens_coordinates

