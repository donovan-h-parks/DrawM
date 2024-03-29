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

import svgwrite

from drawm.svg.svg_utils import render_label
from drawm.tree.newick_utils import parse_label


class LabelProps:
    """Visual attributed for internal and leaf labels."""

    def __init__(self, config_file, dwg, inch):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.inch = inch
        
        self.show_leaf_labels = False
        self.show_internal_labels = False
        
        self.internal_font_size = 10
        self.internal_font_color = 'rgb(0,0,0)'
        self.internal_sample_rate = 1
        
        self.leaf_font_size = 10
        self.leaf_font_color = 'rgb(0,0,0)'
        self.leaf_sample_rate = 1
        
        if not config_file:
            return # use default values

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
                    self.internal_font_size = float(values[0]) * (self.inch/90.0)
                elif attribute == 'internal_font_color':
                    self.internal_font_color = values[0]
                elif attribute == 'internal_sample_rate':
                    self.internal_sample_rate = int(values[0])
                elif attribute == 'leaf_font_size':
                    self.leaf_font_size = float(values[0]) * (self.inch/90.0)
                elif attribute == 'leaf_font_color':
                    self.leaf_font_color = values[0]
                elif attribute == 'leaf_sample_rate':
                    self.leaf_sample_rate = int(values[0])
                else:
                    self.logger.warning('[LabelProps] Unexpected attribute: %s' % attribute)
     
    def _render_internal_labels(self, tree):
        """Render internal labels."""
    
        label_group = svgwrite.container.Group(id='internal_node_labels')
        self.dwg.add(label_group)
        
        node_count = -1
        for node in tree.preorder_node_iter(lambda n: not n.is_leaf()):
            if node.is_collapsed:
                continue
                
            node_count += 1
            if node_count % self.internal_sample_rate:
                continue
                
            _support, taxon, _auxiliary_info = parse_label(node.label)
            offset = 0.02*self.inch
            render_label(self.dwg,
                            node.x + offset*node.x_dir, 
                            node.y + offset*node.y_dir, 
                            node.angle, 
                            taxon, 
                            self.internal_font_size, 
                            self.internal_font_color,
                            middle_y=True,
                            group=label_group,
                            id_prefix='internal_label')

    def _render_leaf_labels(self, tree):
        """"Render labels for extant taxa."""
        
        label_group = svgwrite.container.Group(id='leaf_node_labels')
        self.dwg.add(label_group)
        
        node_count = -1
        for leaf in tree.leaf_node_iter():
            if leaf.is_collapsed:
                continue

            node_count += 1
            if node_count % self.leaf_sample_rate:
                continue
            
            offset = 0.02*self.inch
            render_label(self.dwg, 
                            leaf.x + offset*leaf.x_dir, 
                            leaf.y + offset*leaf.y_dir, 
                            leaf.angle, 
                            leaf.taxon.label, 
                            self.leaf_font_size, 
                            self.leaf_font_color,
                            middle_y=True,
                            group=label_group,
                            id_prefix='leaf_label')
                    
    def render(self, tree):
        """Render labels."""
        
        if self.show_leaf_labels:
            self.logger.info("Rendering leaf labels.")
            self._render_leaf_labels(tree)
        
        if self.show_internal_labels:
            self.logger.info("Rendering internal labels.")
            self._render_internal_labels(tree)

