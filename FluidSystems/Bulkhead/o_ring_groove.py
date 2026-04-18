import numpy as np
from dataclasses import dataclass

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

def Perform_O_Ring_Seal_Analysis(parameters):
    
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
        def current_temperature_width(self):
            current_temperature_width = self.original_temperature_width * (1 + self.temperature_strain)
            return(current_temperature_width)
        
        @property
        def current_temperature_free_outer_diameter(self):
            return self.original_temperature_free_outer_diameter * (1 + ((self.current_temperature - self.original_temperature) * self.material.CTE))
        
        @property
        def current_temperature_free_inner_diameter(self):
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
        def current_temperature_gland_outer_diameter(self):
            current_temperature_gland_outer_diameter = self.original_temperature_gland_outer_diameter * (1 + self.temperature_strain)
            return(current_temperature_gland_outer_diameter)
        
        @property
        def current_temperature_gland_inner_diameter(self):
            current_temperature_gland_inner_diameter = self.original_temperature_gland_inner_diameter * (1 + self.temperature_strain)
            return(current_temperature_gland_inner_diameter)
        
        @property
        def current_temperature_groove_outer_diameter(self):
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

        def Calculate_O_Ring_Squeeze(self):
             return(Calculate_O_Ring_Squeeze(
                                        o_ring_seal_type = self.gland.o_ring_seal_type, 
                                        gland_outer_diameter = self.gland.current_temperature_gland_outer_diameter, 
                                        gland_inner_diameter = self.gland.current_temperature_gland_inner_diameter, 
                                        o_ring_width = self.o_ring.current_temperature_width,
                                     ))
             
        def Calculate_Radial_Seal_Diametral_Clearance(self):
            return(Calculate_Radial_Seal_Diametral_Clearance(
                                                        gland_outer_diameter = self.gland.current_temperature_gland_outer_diameter,
                                                        groove_outer_diameter = self.gland.current_temperature_groove_outer_diameter,
                                                     ))
            
        def Calculate_O_Ring_Fill(self):
            return(Calculate_O_Ring_Fill(
                                    o_ring_seal_type = self.gland.o_ring_seal_type, 
                                    gland_outer_diameter = self.gland.current_temperature_gland_outer_diameter, 
                                    gland_inner_diameter = self.gland.current_temperature_gland_inner_diameter, 
                                    gland_width = self.gland.current_temperature_width,
                                    o_ring_width = self.o_ring.current_temperature_width,
                                 ))    
    
    
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
    squeeze = oxidizer_tank_radial_o_ring_seal.Calculate_O_Ring_Squeeze()
    diametral_clearance = oxidizer_tank_radial_o_ring_seal.Calculate_Radial_Seal_Diametral_Clearance()
    fill = oxidizer_tank_radial_o_ring_seal.Calculate_O_Ring_Fill()
    o_ring_temperature = oxidizer_tank_radial_o_ring_seal.o_ring.current_temperature
    gland_temperature = oxidizer_tank_radial_o_ring_seal.gland.current_temperature
    
    print("\n")
    print(f"-------------------------------------------------------------------------------")
    print(f"| Seal | O-Ring Temperature [°K] | Squeeze [%] | Diametral Clearance [in] | Fill [%] |")
    print(f"-------------------------------------------------------------------------------")
    print(f"| LOx  |{o_ring_temperature:^25.2f}|{squeeze * 100:^13.2f}|{diametral_clearance * c.M2IN:^26.4f}|{fill * 100:^10.2f}|")
    
    
    # print(f"oxidizer_tank_radial_o_ring_seal.o_ring.temperature: {oxidizer_tank_radial_o_ring_seal.o_ring.original_temperature}°K")
    # nitrogen_saturation_temperature = 77 # [k]
    # oxidizer_tank_radial_o_ring_seal.o_ring.original_temperature = nitrogen_saturation_temperature
    # oxidizer_tank_radial_o_ring_seal.o_ring_gland.original_temperature = nitrogen_saturation_temperature
    
    # print(f"\noxidizer_tank_radial_o_ring_seal.o_ring.temperature: {oxidizer_tank_radial_o_ring_seal.o_ring.original_temperature}°K")
    # oxidizer_tank_radial_o_ring_seal.Calculate_O_Ring_Squeeze()
    # oxidizer_tank_radial_o_ring_seal.Calculate_Radial_Seal_Diametral_Clearance()
    
    
def main(parameters):
    parameters = Perform_O_Ring_Seal_Analysis(parameters)
    return(parameters)

if __name__ == "__main__":
    parameters, wet_mass_distribution, dry_mass_distribution = vehicle_parameters.main()
    main(parameters)
