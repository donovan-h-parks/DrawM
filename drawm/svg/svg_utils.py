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
    
    
def donut(dwg, x, y, inner_radius, outer_radius, color, opacity=1.0, group=None, id=None):
    """Render a donut."""
    
    path = dwg.path("M%d,%d" % (x,y), id=id)
    path.fill(color=color, opacity=opacity, rule='evenodd')
    path.stroke(color=color, width=1)
    
    path.push("m%d,%d" % (0, -outer_radius))
    path.push_arc(target=(1, 0), 
                            rotation=0, 
                            r=outer_radius,
                            large_arc=True,
                            angle_dir='-',
                            absolute=False)
    path.push("m%d,%d" % (0, outer_radius-inner_radius))
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
    
def render_label(dwg, x, y, angle, label, font_size, color, group=None):
    """Render label."""
    
    if label is None:
        return
    
    # make sure angle is between -180 and 180
    if angle > 180:
        angle = angle - 360

    # determine rendering angle and direction of text
    text_anchor = 'start'
    direction = 'ltr'        
    if angle < -90 or angle > 90:
        text_anchor = 'end'
        direction = 'rtl'

        if angle < -90:
            angle = 180 + angle
        elif angle > 90:
            angle = angle - 180

    t = dwg.text(label, 
                    x=[(x)], 
                    y=[(y + 0.4*font_size)], 
                    font_size=font_size,
                    text_anchor=text_anchor,
                    direction=direction,
                    fill=color,
                    id=label.replace(' ', '_'))

    t.rotate(angle, (x, y))
    
    if group:
        group.add(t)
    else:
        # add to 'root' group
        dwg.add(t)