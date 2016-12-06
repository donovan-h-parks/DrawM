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
import math

from drawm.svg.svg_utils import render_label
from drawm.tree.newick_utils import parse_label


class LabelProps:
    """Visual attributed for internal and leaf labels."""

    def __init__(self, config_file, dwg, tree_height, inch, font_size):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.tree_height = tree_height
        self.inch = inch
        self.font_size = font_size
        
        self.node_radius = 0.005 * self.tree_height
        
        self.show_leaf_labels = False
        self.show_internal_labels = False
        
        self.internal_font_size = None
        self.internal_font_color = None
        
        self.leaf_font_size = None
        self.leaf_font_color = None

        with open(config_file) as f:
            prop_type = f.readline().strip()
            
            if prop_type != 'LABELS':
                self.logger.error("[LabelProps] Unexpected property type '%s' in %s." % (prop_type, config_file))
                sys.exit()
            
            for line in f:
                if line[0] == '#' or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                attribute, values = fields[0], fields[1:]
                if attribute == 'show_leaf_labels':
                    self.show_leaf_labels = (values[0] == 'True')
                elif attribute == 'show_internal_labels':
                    self.show_internal_labels = (values[0] == 'True')
                elif attribute == 'internal_font_size':
                    self.internal_font_size = int(values[0])
                elif attribute == 'internal_font_color':
                    self.internal_font_color = values[0]
                elif attribute == 'leaf_font_size':
                    self.leaf_font_size = int(values[0])
                elif attribute == 'leaf_font_color':
                    self.leaf_font_color = values[0]
     
    def _render_internal_labels(self, tree):
        """Render internal labels."""
    
        font_size = 0.5*self.font_size
        
        for node in tree.preorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)
            if taxon and taxon in ['p__Patescibacteria', 'p__Firmicutes', 'p__Actinobacteria', 'p__Proteobacteria', 'p__Bacteroidetes']:
                dx = node.x - node.corner_x
                dy = node.y - node.corner_y
                angle = math.atan2(dy, dx)
                angle_deg = math.degrees(angle)
                render_label(self.dwg, 
                                node.x, 
                                node.y, 
                                angle_deg, 
                                taxon, 
                                self.internal_font_size, 
                                self.internal_font_color)

    def _render_leaf_labels(self, tree):
        """"Render labels for extant taxa."""
        
        font_size = 0.5*self.font_size
        
        for i, leaf in enumerate(tree.leaf_node_iter()):
            dx = leaf.x - leaf.corner_x
            dy = leaf.y - leaf.corner_y
            angle = math.atan2(dy, dx)
            angle_deg = math.degrees(angle)
            render_label(self.dwg, 
                            leaf.x, 
                            leaf.y, 
                            angle_deg, 
                            leaf.taxon.label, 
                            self.leaf_font_size, 
                            self.leaf_font_color)
                    
    def render(self, tree):
        """Render labels."""
        
        if self.show_leaf_labels:
            self.logger.info("Rendering leaf labels.")
            self._render_leaf_labels(tree)
        
        if self.show_internal_labels:
            self.logger.info("Rendering internal labels.")
            self._render_internal_labels(tree)

