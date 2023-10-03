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
import math
import logging

import svgwrite

from drawm.tree.newick_utils import parse_label
from drawm.svg.svg_utils import render_label, color_str, rgb_from_str


class BootstrapProps:
    """Visual attributed for support values."""

    def __init__(self, config_file, dwg, inch):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.inch = inch

        self.show_bootstraps = False
        
        self.show_bootstrap_labels = False
        self.min_bootstrap_label = None
        self.font_size = 10
        self.font_color = 'rgb(0,0,0)'
        
        self.discrete_cm = []
        self.continuous_cm = []
        
        if not config_file:
            return # use default values
        
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
                elif attribute == 'min_bootstrap_label':    
                    self.min_bootstrap_label = float(values[0])
                elif attribute == 'font_size':
                    self.font_size = float(values[0]) * (self.inch/90.0)
                elif attribute == 'font_color':
                    self.font_color = values[0]
                elif attribute == 'discrete_cm':
                    support_value, color, node_radius = values
                    self.discrete_cm.append([float(support_value), color, float(node_radius) * self.inch])
                elif attribute == 'continuous_cm':
                    support_value, color, node_radius  = values
                    self.continuous_cm.append([float(support_value), color, float(node_radius) * self.inch])
                else:
                    self.logger.warning('[BootstrapProps] Unexpected attribute: %s' % attribute)
                    
    def _interpolate_color(self, value, min_value, max_value, min_color, max_color):
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
        
    def _linear_interpolate(self, value, min_value, max_value, min_output, max_output):
        """Linear interpolate."""
        
        if value < min_value or value > max_value:
            return None
            
        index = float(value - min_value) / (max_value - min_value)
        output = min_output + index*(max_output-min_output)

        return output
        
    def render_legend(self, tree):
        """Render legend for support values."""
        
        if not self.show_bootstraps:
            return
        
        color_map = None
        if self.discrete_cm:
            color_map = self.discrete_cm
        elif self.continuous_cm:
            color_map = self.continuous_cm
            
        if color_map:
            legend_height = 0
            max_radius = 0
            for item_index, (_support, _color, node_radius) in enumerate(color_map):
                legend_radius = node_radius
                if legend_radius > max_radius:
                    max_radius = legend_radius
                legend_step = max(2*legend_radius, 1.25*(self.font_size)) # 1.25 pixels/pt
                legend_height += legend_step
            
            legendX = 0.1*self.inch
            legendY = self.dwg.canvas_height - self.inch - legend_height

            bs_legend_group = svgwrite.container.Group(id='support_legend')
            self.dwg.add(bs_legend_group)
            for item_index, (support, color, node_radius) in enumerate(color_map):
                legend_radius = node_radius
                legend_step = max(2*legend_radius, 1.25*(self.font_size)) # 1.25 pixels/pt
            
                c = self.dwg.circle(center=(legendX + max_radius, legendY), 
                                    r=legend_radius,
                                    id='support_legend_symbol_%d' % item_index)
                c.fill(color=color)
                c.stroke(color='black')
                bs_legend_group.add(c)

                render_label(self.dwg, 
                                legendX + 2*max_radius + 0.02*self.inch, 
                                legendY, 
                                0, 
                                '>%g%%' % support, 
                                self.font_size, 
                                self.font_color,
                                middle_y=True,
                                group=bs_legend_group,
                                id_prefix='support_legend_text')
                    
                legendY += legend_step

    def render(self, tree):
        """Render node based on bootstrap support."""
        
        if not self.show_bootstraps:
            return
            
        bs_node_group = svgwrite.container.Group(id='support_nodes')
        self.dwg.add(bs_node_group)
        
        bs_text_group = svgwrite.container.Group(id='support_labels')
        self.dwg.add(bs_text_group)
        
        for node in tree.postorder_node_iter():
            if node.is_collapsed:
                continue
                
            support, _taxon, _auxiliary_info = parse_label(node.label)
            if not support:
                continue
            
            # get color
            color = None
            if self.discrete_cm:
                for threshold, color, node_radius in self.discrete_cm:
                    if support > threshold:
                        break
                        
            if self.continuous_cm:
                max_value, max_color, max_radius = self.continuous_cm[0]
                min_value, min_color, min_radius = self.continuous_cm[1]
                color = self._interpolate_color(support, min_value, max_value, min_color, max_color)
                node_radius = self._linear_interpolate(support, min_value, max_value, min_radius, max_radius)
            
            # render node
            if color:
                node_x, node_y = node.x, node.y
                if node.is_collapsed_root:
                    node_x = node.x - node.x_dir*node_radius
                    node_y = node.y - node.y_dir*node_radius
                    
                c = self.dwg.circle(center=(node_x, node_y), 
                                    r=node_radius,
                                    id='support_%g' % support)
                c.stroke(color='black')
                c.fill(color=color)
                bs_node_group.add(c)
                
            if self.show_bootstrap_labels and support > self.min_bootstrap_label:
                label_x = node.x + 0.5*(node.corner_x - node.x)
                label_y = node.y + 0.5*(node.corner_y - node.y)
                
                offset = 0.5*tree.branch_width + 0.02*self.inch
                if node.angle < 90 or node.angle > 270:
                    label_x += (node.y_dir * offset)
                    label_y -= (node.x_dir * offset)
                else:
                    label_x -= (node.y_dir * offset)
                    label_y += (node.x_dir * offset)
                        
                render_label(self.dwg, 
                                label_x, 
                                label_y,
                                node.angle, 
                                '%g' % support, 
                                self.font_size, 
                                self.font_color,
                                middle_x=True,
                                group=bs_text_group,
                                id_prefix='support_text')
