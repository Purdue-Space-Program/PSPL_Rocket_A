def calculate_thermal_strain(alpha, deltaT):
    return alpha*deltaT*0.000001

def calculate_strain_rate(time, strain):
    return time/strain

def main():
    coefficient_of_thermal_expansion = 18.0
    change_in_temperature = 1400.0
    burn_time = 2
    thermal_strain = calculate_thermal_strain(coefficient_of_thermal_expansion,change_in_temperature)
    print(thermal_strain)
    stress_rate = calculate_strain_rate(burn_time, thermal_strain)
    print(stress_rate)

main()
