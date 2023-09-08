import numpy as np

class bounding_box:

    def __init__(self,points):
        
        self.lowerleft = points[0]
        self.upperleft = points[1]
        self.lowerright = points[2]
        self.upperright = points[3]

    def getUpperRightPoint(self):
        return self.upperright
    
    def getUpperLeftPoint(self):
        return self.upperleft

    def getLowerRightPoint(self):
        return self.lowerright

    def getLowerLeftPoint(self):
        return self.lowerleft
    
        
class point:

    def __init__(self,latitude,longitude):
        self.latitude = latitude
        self.longitude = longitude
    
    def getLatitude(self):
        return self.latitude

    def getLongitude(self):
        return self.longitude
    
