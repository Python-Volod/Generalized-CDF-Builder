###########################################################################################################
# Generalized Cumulative Distribution Function builder for devices with discrete number of random outcomes#
###########################################################################################################
import numpy as np
import matplotlib.pyplot as plt 
from decimal import *


# Main variable inputs, except ps present lower
precision = 25 # amount of decimal places precision
q = Decimal(3) # amount of possible outcomes, int > 0
depth_gl = 4 # length of the associated outcome strings, int > 0
inp_mode = "manual" # the way of inputing probabilities associated with outcomes, either manual (hard coded into this file) or iterative for iterative terminal input

# Applying precision 
getcontext().prec = precision
int_q = int(q) # Nessecary for iterative construction


# Main class, representing a portion of a CDF graph associated with a given possible outcome string
class Line_():
    def __init__(self, final_x, ass_p, initial_x=0, depth = 1, prev_acc_ass_p = 1):
        # X coordinates associated with the line
        self.initial_x = Decimal(initial_x)
        self.final_x = Decimal(final_x)
        # Y coordinates associated with the line
        self.initial_y = Decimal(0)
        self.final_y = Decimal(0)
        # Other descriptive qualities of the line, like slope that depends upon the depth at which the line is created
        self.depth = depth
        self.ass_p = Decimal(ass_p)
        self.acc_ass_p = Decimal(prev_acc_ass_p * self.ass_p)

        if self.depth != 2:
            self.slope = Decimal(self.acc_ass_p/(1/(q**(depth-1))))
        elif self.depth == 2:
            self.slope = Decimal(self.ass_p**(self.depth-1)/(1/(q**(depth-1))))
        else:
            self.slope = 1 

        if self.slope == 0: # If the slope of the line is zero, all subsequent divisions of the line are obsolete
            self.iterable = False
        else: 
            self.iterable = True
    
    # Setters and getters 
    def set_y(self, previous = 0):
        self.initial_y = previous
        self.final_y = ((self.final_x - self.initial_x) * self.slope) + self.initial_y

    # Wraper
    def __str__(self):
        return f"Line segment x1:{self.initial_x} x2:{self.final_x}, slope {self.slope}, depth {self.depth}" 
    def __repr__(self):
        return f"Line segment x1:{self.initial_x} x2:{self.final_x}, slope {self.slope}, depth {self.depth}" 
    
    # Comparison operators 

    def __eq__(self, other):
        return self.final_x == other.final_x
    def __lt__(self, other):
        return self.final_x < other.final_x
    def __bt__(self, other):
        return self.final_x > other.final_x
    
    # Hashing

    def __hash__(self):
        return hash(self.initial_x)
    
    # Generating x-points on a given line that would be used to create new lines at the next depth 
    def __defining_points_for_iteration(self):
        base = (self.final_x - self.initial_x) / q
        defining_points = [] 
        # Generate defining points based on the base value
        for i in range(int_q+1):  # Generate q points
            defining_points.append(self.initial_x + (i * base))

        return defining_points
    
    # Actually splits the line into subsequent q Line objects
    def split(self):
        if self.iterable == True:
            resulting_lines = set([])
            def_points = self.__defining_points_for_iteration()
            for i in range(0, len(def_points)-1):
                resulting_lines.add(Line_(def_points[i+1], ps[i%int_q], def_points[i], self.depth + 1, self.acc_ass_p))
            return resulting_lines
        else:
            return set([self])
        
    # Calculates dx integral of the Line object
    def calculate_area(self):
        area_triangle = ((self.final_x - self.initial_x) * (self.final_y - self.initial_y))/2
        area_rectangle = ((self.initial_y) * (self.final_x - self.initial_x))
        return (area_rectangle + area_triangle)
    # Calculates lemgth of the Line segment
    def calculate_arclength(self):
        return Decimal((self.final_x - self.initial_x)**2 + (self.final_y - self.initial_y)**2).sqrt()
    

# Sets probabilities for all q-s
def generate_ps():
    global ps
    ps = []
    if inp_mode == "iterative":
        count = 1
        while count <= q:
            p = input(f"Please enter the P{count} value in a/b form: ").split("/")
            try:
                p = Decimal(int(p[0])/int(p[1]))
            except:
                p = 0.0
            if p > 1:
                raise OverflowError
            ps.append(p)
            count += 1
    else:
        ps = [Decimal(0.5), Decimal(0), Decimal(0.5)] # Adjust manually
    return ps

# Taking two line segments, determines if they have an intersection on the specified range, using Newtonian method
def do_they_intersect(line_segment, range_, ground=(1, 0)):
    m1, b1 = line_segment  # line_segment: y = m1 * x + b1
    m2, b2 = ground        # ground line: y = m2 * x + b2
    x1, x2 = range_        # x range [x1, x2]
    
    # Calculate y-values at the range boundaries for both lines
    y1_line = Decimal(m1 * x1 + b1)
    y2_line = Decimal(m1 * x2 + b1)
    y1_ground = Decimal(m2 * x1 + b2)
    y2_ground = Decimal(m2 * x2 + b2)
    
    intersect = ((y1_line != y1_ground) and (y1_line - y1_ground) * (y2_line - y2_ground) <= 0)
    if intersect == False: intersect = (y2_line==y2_ground)
    if (y1_line==y1_ground) and (y2_ground==y2_line): intersect = "Infinite or all"
    
    return intersect

# Main functional part, takes an interval from 0-1 and splits it depending on probabilities  
def generate_all_defining_lines():
    initial_line = Line_(1, 1, 0, 1)
    graph_segments = set()
    graph_segments.add(initial_line)

    for i in range(depth_gl-1):
        new_segments = set([])
        delete_segments = set([])

        for segment in graph_segments:
            if segment.slope != 0:
                delete_segments.add(segment)
            new_segments.update(segment.split())
        for segment in delete_segments:
            graph_segments.remove(segment)

        graph_segments.update(new_segments)
    return list(graph_segments)

def calculate_l1(segment, b, sol= None, x=None):
    l1 = 0
    if sol != None:
        upper_integral = (sol[x] ** 2)/2 - (segment.initial_x ** 2)/2
        lower_integral = ((sol[x] ** 2) * segment.slope/2 + b * sol[x]) - ((segment.initial_x ** 2) * segment.slope/2 + b * segment.initial_x) 
        l1 += abs(upper_integral - lower_integral)
        upper_integral = (segment.final_x ** 2)/2 - (sol[x] ** 2)/2
        lower_integral = ((segment.final_x ** 2) * segment.slope/2 + b * segment.final_x) - ((sol[x] ** 2) * segment.slope/2 + b * sol[x])  
        l1 += abs(upper_integral - lower_integral)

    else: 
        upper_integral = (segment.final_x ** 2)/2 - (segment.initial_x ** 2)/2
        lower_integral = lower_integral = ((segment.final_x ** 2) * segment.slope/2 + b * segment.final_x) - ((segment.initial_x ** 2) * segment.slope/2 + b * segment.initial_x) 
        l1 += abs(upper_integral - lower_integral)
    return l1

# Generates all the nessecary data for a given depth
def calculate_for_depth():
    generate_ps()
    graph_segments = generate_all_defining_lines()
    graph_segments.sort()
    l1 = 0
    arclength = 0
    integral = 0 

    xs = [0]
    ys = [0]
    for i in range(0, len(graph_segments)):
        if i != 0:
            xs.append(graph_segments[i].final_x)
            graph_segments[i].set_y(ys[i])
            ys.append(graph_segments[i].final_y)
        else:
            xs.append(graph_segments[i].final_x)
            graph_segments[i].set_y()
            ys.append(graph_segments[i].final_y)
    
    for segment in graph_segments:
        arclength += segment.calculate_arclength()
        integral += segment.calculate_area()


    intersections = 2 # Initial and terminal points
    solutions = set()
    i = 0
    for segment in graph_segments:
        b = Decimal(segment.initial_y - segment.slope * segment.initial_x)
        intersection = do_they_intersect((segment.slope, b), (segment.initial_x, segment.final_x))
        if intersection == "Infinite or all":
            intersections = "Infinite or all"
            break
        elif intersection == True and i!=0 and i!=len(graph_segments)-1:
            intersections += 1
        else:
            pass
        i += 1   
    return xs, ys, integral, arclength, intersections, l1



def main():
    xs, ys, integral, arclength, intersections, l1 = calculate_for_depth()
    fig, ax = plt.subplots()
    fig.canvas.manager.set_window_title('Generalized CDF builder with L1 colmogorov-smirnov test')


    
    # Plot the function and y=x reference line
    ax.plot(xs, ys, "c", label="Cantor-like Function")
    ax.plot([0, 1], [0, 1], "r", label="y=x")
    # ax.set_title("CDF for Q = 3; P = (P1 = 0.5, P2 = 0, P3 = 0.5); Iteration: 13")

    
    # Create the stats text to display
    stats_text = (
        f"Integral: {integral}\n"
        f"Arclength: {arclength}\n"
        f"Intersections: {intersections}\n"
        "Area between the two: {Under development}"
    )
    
    ax.text(0.95, 0.05, stats_text, transform=ax.transAxes, verticalalignment='bottom', 
            horizontalalignment='right', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))

    # Show legend and plot
    ax.legend()

    plt.show()

main()
