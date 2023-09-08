import numpy as np

class bounding_box:

    def __init__(self,points):
        
        self.lower_left = points[0]
        self.upper_left = points[1]
        self.lower_right = points[2]
        self.upper_right = points[3]

    def get_upper_right_point(self):
        return self.upper_right
    
    def get_upper_left_point(self):
        return self.upper_left

    def get_lower_right_point(self):
        return self.lower_right

    def get_lower_left_point(self):
        return self.lower_left
    
        
class point:

    def __init__(self,latitude,longitude):
        self.latitude = latitude
        self.longitude = longitude
    
    def get_latitude(self):
        return self.latitude

    def get_longitude(self):
        return self.longitude
    
