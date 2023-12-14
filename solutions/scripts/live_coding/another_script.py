import argparse
import math

class Circle:

    def __init__(self, radius):
        self.radius = radius

    def compute_area(self):
        return self.radius * self.radius * math.pi
    
    def compute_circumference(self):
        return self.radius * 2 * math.pi
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simple Circle Calculator",
                                     epilog="Developed by Benjamin Rudski, Oct 2023.")
    
    parser.add_argument("radius", type=float, help="Radius of circle (floating point number).")
    parser.add_argument("-a", "--no-area", action="store_false", help="Indicate whether to compute area.")
    parser.add_argument("-c", "--circumference", action="store_true", help="Indicate whether to compute circumference.")

    args = parser.parse_args()

    circle_radius = args.radius

    circle = Circle(radius=circle_radius)

    should_calculate_area = args.no_area
    should_calculate_circumference = args.circumference

    if should_calculate_area:
        area = circle.compute_area()
        print(f"Circle has area {area}")

    if should_calculate_circumference:
        circumference = circle.compute_circumference()
        print(f"Circle has circumference {circumference}")