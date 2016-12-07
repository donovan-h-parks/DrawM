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
    
    
def render_label(dwg, x, y, angle, label, font_size, color, group=None):
    """Render label."""
    
    # make sure angle is between -180 and 180
    if angle > 180:
        angle = angle - 360

    # determine text offset for best visual quality
    if -45 < angle <= 0:
        x_offset_dir = 1
        y_offset_dir = 1
    elif -90 < angle <= -45:
        x_offset_dir = 1
        y_offset_dir = -1
    elif -135 < angle <= -90:
        x_offset_dir = -1
        y_offset_dir = -1
    elif -180 < angle <= -135:
        x_offset_dir = -1
        y_offset_dir = 1
    elif 135 < angle <= 180:
        x_offset_dir = -1
        y_offset_dir = 1
    elif 90 < angle <= 135:
        x_offset_dir = 1
        y_offset_dir = 1
    elif 45 < angle <= 90:
        x_offset_dir = -1
        y_offset_dir = 1
    elif 0 < angle <= 45:
        x_offset_dir = 1
        y_offset_dir = 1

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

    label_x = x + x_offset_dir*0.5*font_size
    label_y = y + y_offset_dir*0.35*font_size
    t = dwg.text(label, 
                    x=[(label_x)], 
                    y=[(label_y)], 
                    font_size=font_size,
                    text_anchor=text_anchor,
                    direction=direction,
                    #dominant_baseline="central",
                    #dy=['5.3em'],
                    fill=color,
                    id=label)

    t.rotate(angle, (label_x, label_y))
    
    if group:
        group.add(t)
    else:
        # add to 'root' group
        dwg.add(t)