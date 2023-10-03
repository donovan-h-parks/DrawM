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


def color_str(r, g, b):
    return "rgb(%d,%d,%d)" % (int(r+0.5), int(g+0.5), int(b+0.5))
    

def rgb_from_str(color_str):
    rgb = color_str[4:-1]
    return [int(c) for c in rgb.split(',')]
    
    
def donut(dwg, x, y, inner_radius, outer_radius, color, opacity=1.0, group=None, id=None):
    """Render a donut."""
    
    path = dwg.path("M%f,%f" % (x,y), id=id)
    path.fill(color=color, opacity=opacity, rule='evenodd')
    #path.stroke(color=color, width=1)
    
    path.push("m%f,%f" % (0, -outer_radius))
    path.push_arc(target=(1, 0), 
                            rotation=0, 
                            r=outer_radius,
                            large_arc=True,
                            angle_dir='-',
                            absolute=False)
    path.push("m%f,%f" % (0, outer_radius-inner_radius))
    path.push_arc(target=(-1, 0), 
                            rotation=0, 
                            r=inner_radius,
                            large_arc=True,
                            angle_dir='+',
                            absolute=False)

    if group:
        group.add(path)
    else:
        # add to 'root' group
        dwg.add(t)
    
    
def render_label(dwg, x, y, 
                    angle, 
                    label, 
                    font_size, 
                    color,
                    middle_x=False,
                    middle_y=False,                    
                    group=None, 
                    id_prefix=None):
    """Render label."""
    
    if label is None:
        return
    
    # make sure angle is between -180 and 180
    if angle > 180:
        angle = angle - 360

    # determine rendering angle and direction of text
    direction = 'ltr'   
    text_anchor = 'start'
    if angle < -90 or angle > 90:
        text_anchor = 'end'
        direction = 'rtl'

        if angle < -90:
            angle = 180 + angle
        elif angle > 90:
            angle = angle - 180
            
    if middle_x:
        text_anchor = 'middle'
       
    y_label = y
    if middle_y:
        y_label += 0.45*font_size
                   
    id_label = label.replace(' ', '_')
    if id_prefix:
        id_label = '%s_%s' % (id_prefix, id_label)

    t = dwg.text(label, 
                    x=[(x)], 
                    y=[(y_label)],
                    font_size="%fpt" % font_size,
                    text_anchor=text_anchor,
                    direction=direction,
                    fill=color,
                    id=id_label)

    t.rotate(angle, (x, y))

    if group:
        group.add(t)
    else:
        # add to 'root' group
        dwg.add(t)