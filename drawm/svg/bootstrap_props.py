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

import os
import sys
import logging

import svgwrite

from drawm.tree.newick_utils import parse_label


class BootstrapProps:
    """Visual attributed for support values."""

    def __init__(self, config_file, dwg, tree_height, inch, font_size):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.tree_height = tree_height
        self.inch = inch
        self.font_size = font_size
        
        self.node_radius = 0.005 * self.tree_height
         
        self.show_bootstraps = False
        self.show_bootstrap_labels = False
        self.discrete_cm = []
        self.continuous_cm = []
        
        with open(config_file) as f:
            prop_type = f.readline().strip()
            
            if prop_type != 'BOOTSTRAP':
                self.logger.error("[BootstrapProps] Unexpected property type '%s' in %s." % (prop_type, config_file))
                sys.exit()
            
            for line in f:
                if line[0] == '#' or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                attribute, values = fields[0], fields[1:]
                if attribute == 'show_bootstraps':
                    self.show_bootstraps = (values[0] == 'True')
                elif attribute == 'show_bootstrap_labels':
                    self.show_bootstrap_labels = (values[0] == 'True')
                elif attribute == 'discrete_cm':
                    support_value, color = values
                    self.discrete_cm.append([float(support_value), color])
                elif attribute == 'continuous_cm':
                    support_value, color = values
                    self.continuous_cm.append([float(support_value), color])
                else:
                    self.logger.warning('[BootstrapProps] Unexpected attribute: %s' % attribute)
                    
    def _linear_interpolate(self, value, min_value, max_value, min_color, max_color):
        """Linear interpolate color."""
        
        if value < min_value or value > max_value:
            return None
            
        index = float(value - min_value) / (max_value - min_value)
        
        r_min, g_min, b_min = rgb_from_str(min_color)
        r_max, g_max, b_max = rgb_from_str(max_color)
        
        r = r_min + index*(r_max-r_min)
        g = g_min + index*(g_max-g_min)
        b = b_min + index*(b_max-b_min)
        
        return color_str(r,g,b)
        
    def render_legend(self):
        """Render legend for support values."""
        
        if not self.show_bootstraps:
            return
        
        if self.discrete_cm:
            legend_radius = 2*self.node_radius
            legend_step = 2*2*legend_radius
            
            legendX = 0.1*self.inch
            legendY = self.dwg.canvas_height - self.inch - legend_step*len(self.discrete_cm)

            bs_legend_group = svgwrite.container.Group(id='bootstrap_legend')
            self.dwg.add(bs_legend_group)
            for item_index, (support, color) in enumerate(self.discrete_cm):
                c = self.dwg.circle(center=(legendX, legendY), 
                                    r=legend_radius,
                                    id='bootstrap_legend_symbol_%d' % item_index)
                c.fill(color=color)
                c.stroke(color='black')
                bs_legend_group.add(c)
                
                t = self.dwg.text('>%g%%' % support, 
                                    x=[(legendX + 1.5*legend_radius)], 
                                    y=[(legendY + 0.35*self.font_size)], 
                                    font_size=self.font_size,
                                    fill='black',
                                    id='bootstrap_legend_label_%d' % item_index)
                bs_legend_group.add(t)
                    
                legendY += legend_step

    def render(self, tree):
        """Render node based on bootstrap support."""
        
        if not self.show_bootstraps:
            return
            
        bs_node_group = svgwrite.container.Group(id='bootstrap_node')
        self.dwg.add(bs_node_group)
        
        for node in tree.postorder_node_iter():
            if node.is_collapsed:
                continue
                
            support, _taxon, _auxiliary_info = parse_label(node.label)
            if not support:
                continue
            
            # get color
            color = None
            if self.discrete_cm:
                for threshold, color in self.discrete_cm:
                    if support > threshold:
                        color = color
                        break
                        
            if self.continuous_cm:
                max_value, max_color = self.continuous_cm[0]
                min_value, min_color = self.continuous_cm[1]
                color = self._linear_interpolate(support, min_value, max_value, min_color, max_color)
            
            # render node
            if color:
                c = self.dwg.circle(center=(node.x, node.y), 
                                    r=self.node_radius,
                                    id='support_%g' % support)
                c.stroke(color='black')
                c.fill(color=color)
                c = bs_node_group.add(c)
            