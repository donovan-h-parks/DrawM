###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'


def polygon_centroid(pts):
    """Calculate centroid for polygon."""
    
    signed_area = 0
    c_x = 0
    c_y = 0
    for i in xrange(0, len(pts)):
        x_i = pts[i][0]
        y_i = pts[i][1]
        if i+1 == len(pts):
            x_ii = pts[0][0]
            y_ii = pts[0][1]
        else:
            x_ii = pts[i+1][0]
            y_ii = pts[i+1][1]
        
        cross = (x_i*y_ii - x_ii*y_i)
        signed_area += cross
        c_x += (x_i + x_ii)*cross
        c_y += (y_i + y_ii)*cross

    signed_area = 0.5*signed_area
    c_x = c_x / (6*signed_area)
    c_y = c_y / (6*signed_area)
    
    return c_x, c_y
    
def length(v):
    return math.sqrt(v[0]**2+v[1]**2)
    
def dot_product(v,w):
   return v[0]*w[0]+v[1]*w[1]
   
def determinant(v,w):
   return v[0]*w[1]-v[1]*w[0]
   
def inner_angle(v,w):
   cosx=dot_product(v,w)/(length(v)*length(w))
   rad=math.acos(cosx) # in radians
   return rad*180/math.pi # returns degrees
   
def angle_clockwise(A, B):
    inner=inner_angle(A,B)
    det = determinant(A,B)
    if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
        return inner
    else: # if the det > 0 then A is immediately clockwise of B
        return 360-inner
        
def project_pt_to_line(line1, line2, pt):
    m = (line2[1] - line1[1]) / (line2[0] - line1[0])
    b = line1[1] - (m * line1[0])

    x = (m * pt[1] + pt[0] - m * b) / (m * m + 1)
    y = (m * m * pt[1] + m * pt[0] + b) / (m * m + 1)

    return (x, y)
