import pandas as pd
from math import ceil, sqrt
import numpy as np

LB2KG = 0.453592
FT2M = 0.3048
IN2M = 0.0254

def codeFriendlyName(str, delimiter='_'):
    return delimiter.join(str.split(' (')[0].split(' ')).lower() # Delimiter is character between words ex. "_"

def getRocketSections(xlsx):
    rocket_sections = {}
    rocket_inputs = pd.read_excel(xlsx, sheet_name='Vehicle Sections') # Opens and reads a sheet called "Vehicle Sections"

    for section in rocket_inputs.iloc: # Iterates through each row (rocket section) of the sheet
        sec_name = codeFriendlyName(section['Name'])
        # if sec_name == 'nosecone': continue # Nosecone is ejected in recovery
        # print(sec_name) # TEST
        section_dict = section.drop('Name').to_dict() # Each rocket section will have the format { mass: xx, length: xx, ... }

        rocket_sections[sec_name] = {codeFriendlyName(k):v for k,v in section_dict.items() if v==v} # Creates dictionary for each rocket section, eliminates dictionary items wihout values (empty cells in spreadsheet)
        # Convert units to metric
        if 'mass' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['mass'] *= LB2KG
        if 'prop_mass' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['mass'] += 0.08 * rocket_sections[sec_name]['prop_mass'] * LB2KG # Need to know what to change 0.08 to and how much prop is left at max q
        if 'length' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['length'] *= FT2M
        if 'thickness' in rocket_sections[sec_name]:
            rocket_sections[sec_name]['thickness'] *= IN2M

    return rocket_sections

def getPointMasses(xlsx):
    point_mass_inputs = pd.read_excel(xlsx, sheet_name='Point Masses', header=0, usecols=[1,2])
    point_masses = []
    for point_mass in point_mass_inputs.iloc:
        mass, xcoord = list(point_mass)
        point_masses.append({'mass':mass * LB2KG, 'x_coordinate':xcoord * FT2M}) # Only point mass is batteries
    # print(point_masses) # TEST
    return point_masses

def getMassModel(rocket_dict, point_masses, dy):
    element_length = dy # Set element length of 5 mm
    mass_model = []
    model_index = 0
    pointMassIndices = [ceil(pm['x_coordinate']/element_length) for pm in point_masses]
    # print(pointMassIndices) # TEST
    for section in reversed(rocket_dict.values()):
        secMass = section.get('mass')
        secLength = section.get('length')
        if not secMass: continue # Skips fins
        elementsPerSection = ceil(secLength / element_length)
        massPerElement = secMass / elementsPerSection
        for i in range(elementsPerSection):
            if model_index in pointMassIndices:
                point_mass_mass = point_masses[pointMassIndices.index(model_index)]['mass']
                mass_model.append(massPerElement + point_mass_mass)
            else:
                mass_model.append(massPerElement)
            model_index += 1
    # print(mass_model) # TEST
    return np.array(mass_model)

def calculateTotalMass(rocket_dict, point_masses):
    totalPointMass = sum([pm['mass'] for pm in point_masses]) # Sum point masses
    totalSectionMass = sum([rocket_section.get('mass', 0) for rocket_section in rocket_dict.values()]) # Sum section masses
    return totalPointMass + totalSectionMass

def calculateTotalLength(rocket_dict, point_masses):
    totalLength = sum([rocket_section.get('length', 0) for rocket_section in rocket_dict.values()]) # Sum section lengths
    return totalLength

def calculateCG(rocket_dict, point_masses, totalMass):
    moment = 0
    length = 0
    for section in reversed(rocket_dict.values()):
        secMass = section.get('mass', 0)
        secLength = section.get('mass', 0)
        xCoord = length + secLength / 2 # Each section can be represented as a point mass
        # Multiply masses by x coordinate to get moment for each section
        moment += secMass * xCoord
        length += secLength # Go to next section
    
    pointMass = np.array([pm['mass'] for pm in point_masses])
    pointCoord = np.array([pm['x_coordinate'] for pm in point_masses])
    # Multiply masses by x coordinate to get moment for each point mass
    pointMoments = pointMass * pointCoord
    moment += np.sum(pointMoments)

    cg = moment / totalMass # Center of gravity calculation
    return cg

# def getFinSD(rocket_dict, rocket_diameter):

# def getFinCP(noseconeToFin, fin_obj, length, sweep_length):




sfd_inputs = pd.ExcelFile('sfd_inputs.xlsx', engine='openpyxl')
rocket_dict = getRocketSections(sfd_inputs)
point_masses = getPointMasses(sfd_inputs)

def returnMassModel(dy):
    mass_model = getMassModel(rocket_dict, point_masses, dy)
    return mass_model

# print(rocket_dict) # TEST
# print(point_masses) # TEST
# print(returnMassModel(0.005)) # TEST
