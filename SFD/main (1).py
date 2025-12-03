'''

Welcome to GDB Online.
GDB online is an online compiler and debugger tool for C, C++, Python, Java, PHP, Ruby, Perl,
C#, OCaml, VB, Swift, Pascal, Fortran, Haskell, Objective-C, Assembly, HTML, CSS, JS, SQLite, Prolog.
Code, Compile, Run and Debug online from anywhere in world.

'''
import math
''' Fruity Chutes specs '''

rocket_mass = 69
gravity = 9.81
air_density = 1.81
drag_coefficent = 2.2
canopy_area = 24.6677824063
max_height = 3500

weight_force = gravity*rocket_mass # W=mg calc 

terminal_velocity = math.sqrt((2*rocket_mass*gravity)/(canopy_area*drag_coefficent*air_density)) # solving for velocity setting weight and Drag equal

drag_force = drag_coefficent*air_density*canopy_area*((terminal_velocity**2)/2) # Formula from NASA website

decent_time = max_height/terminal_velocity

print ('Termial Velocity: ', terminal_velocity, 'm/s')
print ("Decent Time: ", decent_time, 'seconds')
print ('Drag Force: ', drag_force, 'N')
print ('Weight Force: ', weight_force, 'N')


''' Sphereacutes specs '''
rocket_mass = 69
gravity = 9.81
air_density = 1.81
drag_coefficent1 = .75
canopy_area1 = 15.029593683
max_height = 3500

terminal_velocity1 = math.sqrt((2*rocket_mass*gravity)/(canopy_area1*drag_coefficent1*air_density))
decent_time1 = max_height/terminal_velocity1

print('Terminal Velocity: ', terminal_velocity1, 'm/s')
print('Decent Time: ',decent_time1, 'seconds')




