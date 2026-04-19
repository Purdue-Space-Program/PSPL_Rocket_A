import numpy as np
from dataclasses import dataclass
from tabulate import tabulate

import sys
import os
from pathlib import Path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.join(os.path.dirname(__file__), "..", ".."))))

import constants as c
import vehicle_parameters_functions
import vehicle_parameters
import vehicle_main
import print_filter

from pint import UnitRegistry
u = UnitRegistry()


# Written by David Gustafsson for performing analysis on tank bulkhead O-ring seals
# Resources:
# Parker O-Ring Handbook: https://www.parker.com/content/dam/Parker-com/Literature/O-Ring-Division-Literature/ORD-5700.pdf

# Terminology:
# O-Ring grooves are the machined groove in the male side of a seal (the side that carries the o-ring)
# O-Ring glands are the entire area that the o-ring exists in (so including the gap between bore diameter and the plug diameter)


def Calculate_Circle_Area(diameter):
    circle_radius = (diameter/2)
    circle_area = np.pi*(circle_radius**2)
    return (circle_area)

def Perform_O_Ring_Seals_Analysis(parameters):
    
    def Calculate_O_Ring_Squeeze(o_ring_seal_type, gland_outer_diameter, gland_inner_diameter, o_ring_width):
        match o_ring_seal_type:
            case "radial male":
                gland_depth = (gland_outer_diameter - gland_inner_diameter) / 2
                squeeze = (o_ring_width - gland_depth) / o_ring_width # percentage squeeze
            
            case "face":
                raise ValueError("this o_ring_seal_type is not programmed yet")
            case _:
                raise ValueError("unsupported o_ring_seal_type")
    
        print(f"squeeze: {squeeze * 100:.2f} %")
        
        return(squeeze)
    
    def Calculate_Radial_Seal_Diametral_Clearance(gland_outer_diameter, groove_outer_diameter):
        diametral_clearance = gland_outer_diameter - groove_outer_diameter
        
        print(f"diametral_clearance: {diametral_clearance * c.M2IN:.4f} inches")
        return(diametral_clearance)
    
    
    def Calculate_O_Ring_Fill(o_ring_seal_type, gland_outer_diameter, gland_inner_diameter, gland_width, o_ring_width):
        match o_ring_seal_type:
            case "radial male":
                gland_depth = (gland_outer_diameter - gland_inner_diameter) / 2            
            case "face":
                raise ValueError("this o_ring_seal_type is not programmed yet")
            case _:
                raise ValueError("unsupported o_ring_seal_type")
    
        gland_cross_sectional_area = gland_depth * gland_width
        o_ring_cross_sectional_area = Calculate_Circle_Area(diameter = o_ring_width)
        
        o_ring_fill = o_ring_cross_sectional_area/gland_cross_sectional_area
        print(f"o_ring_fill: {o_ring_fill * 100:.2f} %")
        
        return(o_ring_fill)
    
    def Calculate_O_Ring_Stretch(free_inner_diameter, fixed_inner_diameter):
        # empirical formulas from page 54 of parker o-ring handbook
        
        diameter_stretch_percent = ((fixed_inner_diameter - free_inner_diameter) / free_inner_diameter) * 100
        
        X = diameter_stretch_percent
        
        # i put a negative in front of the original formula because the change is a decrease
        if ((X >= 0) and (X <= 3)):
            o_ring_width_change_percent = -(-0.005 + 1.19*X - 0.19*(X**2) - 0.001*(X**3) + 0.008*(X**4))

        elif ((X > 3) and (X <= 25)):
            o_ring_width_change_percent = -(0.56 + 0.59*X - 0.0046*(X**2))
        else:
            raise ValueError("fuck")
        
        # the empirical formula will give a non-zero o_ring_width_change_percent for a stretch of zero
        o_ring_width_change_percent = min(0.0, o_ring_width_change_percent)
        
        return(o_ring_width_change_percent)
        
    @dataclass
    class Material:
        name: str
        CTE: float
    
    @dataclass
    class ORing:
        original_temperature_width: float # referred to as the Cross-Section in the diagram of page 89 of Parker O-Ring Handbook
        original_temperature_free_outer_diameter: float
        original_temperature_free_inner_diameter: float
        # squeeze??????????
        original_temperature: float # temperature at which the dimensions are given
        current_temperature: float # temperature at which the dimensions are given
        material: Material
        
        @property
        def temperature_strain(self):
            temperature_change = self.current_temperature - self.original_temperature
            temperature_strain = temperature_change * self.material.CTE
            return(temperature_strain)
        
        @property
        def free_outer_diameter(self):
            return self.original_temperature_free_outer_diameter * (1 + ((self.current_temperature - self.original_temperature) * self.material.CTE))
        
        @property
        def free_inner_diameter(self):
            return self.original_temperature_free_inner_diameter * (1 + ((self.current_temperature - self.original_temperature) * self.material.CTE))
    
    @dataclass
    class ORingGland:
        o_ring_seal_type: str
        original_temperature_gland_outer_diameter: float # referred to as Bore Diameter or A dia. in the diagram on page 89 of Parker O-Ring Handbook
        original_temperature_gland_inner_diameter: float # referred to as Groove Diameter or B-1 dia. in the diagram on page 89 of Parker O-Ring Handbook
        original_temperature_groove_outer_diameter: float # referred to as Plug Diameter or C dia. in the diagram on page 89 of Parker O-Ring Handbook
        original_temperature_width: float # referred to as G in the diagram on page 89 of Parker O-Ring Handbook
        original_temperature: float # temperature at which the dimensions are given
        current_temperature: float # temperature at which the dimensions are given
        material: Material

        @property
        def temperature_strain(self):
            temperature_change = self.current_temperature - self.original_temperature
            temperature_strain = temperature_change * self.material.CTE
            return(temperature_strain)
        
        @property
        def gland_outer_diameter(self):
            current_temperature_gland_outer_diameter = self.original_temperature_gland_outer_diameter * (1 + self.temperature_strain)
            return(current_temperature_gland_outer_diameter)
        
        @property
        def gland_inner_diameter(self):
            current_temperature_gland_inner_diameter = self.original_temperature_gland_inner_diameter * (1 + self.temperature_strain)
            return(current_temperature_gland_inner_diameter)
        
        @property
        def groove_outer_diameter(self):
            current_temperature_groove_outer_diameter = self.original_temperature_groove_outer_diameter * (1 + self.temperature_strain)
            return(current_temperature_groove_outer_diameter)
        
        @property
        def current_temperature_width(self):
            current_temperature_width = self.original_temperature_width * (1 + self.temperature_strain)
            return(current_temperature_width)
        

    @dataclass
    class RadialORingSeal:
        o_ring: ORing
        gland: ORingGland

        # you need values from both the o-ring and the gland to calculate this because of squeeze
        @property
        def o_ring_width(self):
            
            current_temperature_o_ring_width = self.o_ring.original_temperature_width * (1 + self.o_ring.temperature_strain)
            
            o_ring_stretch_width_change_percent = Calculate_O_Ring_Stretch(free_inner_diameter = self.o_ring.free_inner_diameter, fixed_inner_diameter = self.gland.gland_inner_diameter)            
            current_temperature_and_squeezed_width = current_temperature_o_ring_width * (1 + (o_ring_stretch_width_change_percent/100))
            
            return(current_temperature_and_squeezed_width)


        def Calculate_O_Ring_Squeeze(self):
             return(Calculate_O_Ring_Squeeze(
                                        o_ring_seal_type = self.gland.o_ring_seal_type, 
                                        gland_outer_diameter = self.gland.gland_outer_diameter, 
                                        gland_inner_diameter = self.gland.gland_inner_diameter, 
                                        o_ring_width = self.o_ring_width,
                                     ))
             
        def Calculate_Radial_Seal_Diametral_Clearance(self):
            return(Calculate_Radial_Seal_Diametral_Clearance(
                                                        gland_outer_diameter = self.gland.gland_outer_diameter,
                                                        groove_outer_diameter = self.gland.groove_outer_diameter,
                                                     ))
            
        def Calculate_O_Ring_Fill(self):
            return(Calculate_O_Ring_Fill(
                                    o_ring_seal_type = self.gland.o_ring_seal_type, 
                                    gland_outer_diameter = self.gland.gland_outer_diameter, 
                                    gland_inner_diameter = self.gland.gland_inner_diameter, 
                                    gland_width = self.gland.current_temperature_width,
                                    o_ring_width = self.o_ring_width,
                                 ))
            
        def Analyze(self):
            squeeze = self.Calculate_O_Ring_Squeeze()
            diametral_clearance = self.Calculate_Radial_Seal_Diametral_Clearance()
            fill = self.Calculate_O_Ring_Fill()
            o_ring_temperature = self.o_ring.current_temperature
            gland_temperature = self.gland.current_temperature
            
            output = ["LOx", o_ring_temperature, squeeze*100, diametral_clearance*c.M2IN, fill*100]
            return(output)
    
    
    PTFE = Material(
                    name = "PTFE",
                    CTE = c.PTFE_CTE,
                    )
    
    Aluminum = Material(
                    name = "Aluminum",
                    CTE = c.AlUMINUM_CTE,
                    )
    
    
    o_ring_254 = ORing(original_temperature_width = 0.139 * c.IN2M,
                        original_temperature_free_outer_diameter = 5.762 * c.IN2M,
                        original_temperature_free_inner_diameter = 5.484 * c.IN2M,
                        original_temperature = c.T_AMBIENT,
                        current_temperature = c.T_AMBIENT,
                        material = PTFE,
                        )
    
    bulkhead_o_ring_gland = ORingGland( 
                                        o_ring_seal_type = "radial male",
                                        original_temperature_gland_outer_diameter = 5.750 * c.IN2M,
                                        original_temperature_gland_inner_diameter = 5.528 * c.IN2M,
                                        original_temperature_groove_outer_diameter = 5.747 * c.IN2M,
                                        original_temperature_width = 0.187 * c.IN2M,
                                        original_temperature = c.T_AMBIENT,
                                        current_temperature = c.T_AMBIENT,
                                        material = Aluminum,
                                       )
    
    oxidizer_tank_radial_o_ring_seal = RadialORingSeal(
                                            o_ring = o_ring_254,
                                            gland = bulkhead_o_ring_gland,)
    
    # fuel_tank_o_ring_seal = RadialORingSeal()

 
    print(f"\nOxidizer Tank Radial O-ring Seal:")
    
    # room temperature installation
    output_table_data = []
    output_table_headers = ["Seal Name", "O-Ring Temperature [°K]", "Squeeze [%]", "Diametral Clearance [in]", "Fill [%]"]
    
    output = oxidizer_tank_radial_o_ring_seal.Analyze()
    output_table_data.append(output)

    # when fully in cryogenic fluid
    nitrogen_saturation_temperature = 77 # [°K]
    oxidizer_tank_radial_o_ring_seal.o_ring.current_temperature = nitrogen_saturation_temperature
    oxidizer_tank_radial_o_ring_seal.gland.current_temperature = nitrogen_saturation_temperature
    
    output = oxidizer_tank_radial_o_ring_seal.Analyze()
    output_table_data.append(output)
    
    
    print(tabulate(output_table_data, headers=output_table_headers, colglobalalign=("center"), tablefmt="grid"))
    # print(tabulate(output_table_data, headers=output_table_headers, tablefmt="pipe"))
    
    
def main(parameters):
    parameters = Perform_O_Ring_Seals_Analysis(parameters)
    return(parameters)

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    main(parameters)
