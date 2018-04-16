#!/usr/bin/env python3 
# Name: Andrew Zarzar (azarzar) 
# Group Members: None

'''
Program uses the Triad class to perform calculations that returns the bond length between point1 and point2, and 
the bond length between point 2 and point3. The bond angle between point1-point2-point3 is then returned. Triad methods are used
for all calculations and main() is used for proper formating of inputs. 

Assumptions: There are only 3 points.
'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))

"""Houses the semantics of the program. Sorts through input, and performs calculations using Triad class"""
def main():
    data = input("Please enter 3 atomic coordinates on the same line: \n")
    data = data.split("=") #splits input at "=" to start the seperation of points
    data = "".join(data) # data is joined back into a string after the split
    data = data.replace(" ","") #data removes all white space from the string.
    
    #indexes the location of the first point data using the locations of parenthesis
    p1start = data.index("(")
    p1end = data.index(")")
    #indexes the location of the second point data using the locations of parenthesis and first point data
    p2start = data.index("(", p1end)
    p2end = data.index(")", p2start)
    #indexes the location of the third point data using the locations of parenthesis and second point data
    p3start = data.index("(", p2end)
    p3end = data.index(")", p3start)
    #The names of the elements corresponding to the points are found and stored using the relative location of the parenthesis.
    Element1 = data[0:p1start]
    Element2 = data[p1end+1:p2start]
    Element3 = data[p2end+1:p3start]
    #stores the values of coordinates without including parenthesis
    coord1 = data[p1start+1:p1end]
    coord2 = data[p2start+1:p2end]
    coord3 = data[p3start+1:p3end]
    #splits string of each coordinate at commas, to make a list for each of the coordinates
    Lcoord1 = coord1.split(',')
    Lcoord2 = coord2.split(',')
    Lcoord3 = coord3.split(',')
    # casts the type of of every value in the 3 lists from a string into a float
    for i in range(len(Lcoord1)):
        Lcoord1[i] = float(Lcoord1[i])
    for i in range(len(Lcoord2)):
        Lcoord2[i] = float(Lcoord2[i])
    for i in range(len(Lcoord1)):
        Lcoord3[i] = float(Lcoord3[i])
        
    """#debugging test
    print (data)
    print (Element1)
    print (Element2)
    print (Element3)
    print(coord1)
    print(coord2)
    print(coord3)
    print(Lcoord1)
    print(Lcoord2)
    print(Lcoord3)
    """
    #the float values are put into the triad class
    coords = Triad(Lcoord1,Lcoord2,Lcoord3)
    # triad methods are used to calculated bond lengths between PQ and QR. Bond angle is also calcuated
    bondLengthqp = coords.dPQ()
    bondLengthqr = coords.dQR()
    bondAngle = math.degrees(coords.angleQ())
    print("")
    #Results are printed using correct formating
    print (Element2 +"-"+Element1 + " bond length = "+ "{0:0.2f}".format(bondLengthqp))
    print (Element2 +"-"+Element3 + " bond length = "+ "{0:0.2f}".format(bondLengthqr))
    print (Element1 +"-"+Element2+"-"+Element3 + " bond angle = "+ "{0:0.1f}".format(bondAngle))
    
main()