import PCFBmodule as PCFB
import math as m

#Here we save into variables the lists returned by the parsePDB function.
coordinates, residues, heterogens = PCFB.parsePDB('/Users/nuno_chicoria/Spyder_test/5kkk.pdb')

#Here we save into a variable the list returned by the distance2Center function.
distance2Center = PCFB.distance2Center(coordinates)

#Here we create a list with 3 letter code for all the hydrophobic aminoacids.
hydrophobic = ['ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'PRO', 'VAL', 'TRP']

#Here we create a new file where we save 1. the aminoacid name 2. its distance
#to the center and 3. if it is a hydrophobic (1) or non-hydrophobic (0) aminoacid.
file = open('distance2Center.csv', 'w')
header = "Residue name;Distance2Center;Hidrophobic\n"
file.write(header)
for i in range(len(residues)):
    name_temp = residues[i]
    dist_temp = distance2Center[i]
    if residues[i] in hydrophobic:
        str_temp = "%(name_temp)s;%(dist_temp)f;1\n" % {'name_temp': name_temp, "dist_temp": dist_temp}
        file.write(str_temp)
    else:
        str_temp = "%(name_temp)s;%(dist_temp)f;0\n" % {'name_temp': name_temp, "dist_temp": dist_temp}
        file.write(str_temp)
file.close()

#Here we save into a variable the list returned by the function searchHeterogen.
heterogens_coordinates = PCFB.searchHeterogen(heterogens, "FE")

#We save into a list the distance of all aminoacids to the FE atom.
distance2Fe = []
for coordinates in coordinates:
    x_temp = coordinates[0]
    y_temp = coordinates[1]
    z_temp = coordinates[2]
    distance = m.sqrt((heterogens_coordinates[0][0] - x_temp)**2 + (heterogens_coordinates[0][1] - y_temp)**2 + (heterogens_coordinates[0][2] - z_temp)**2)
    distance2Fe.append(distance)

#We create a tuple with the distance and its index and we sort them by distance.
dist_index = zip(distance2Fe, range(len(distance2Fe)))
sorted_list = sorted(dist_index, key=lambda v: v[0])

#For the 5 closest residues to the FE atom we print 1. the residue number
#2. its distance to the Fe atom and 3. its 3 code letter name.
for i in range(0, 5):
    print("The residue %(residue_number)f is %(residue_dist)f ångströms from the Fe atom and it's a %(residue_name)s.\n" % {'residue_number': (sorted_list[i][1] + 1), "residue_dist": sorted_list[i][0], 'residue_name': residues[sorted_list[i][1]]})
