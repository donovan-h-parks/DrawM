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
from collections import defaultdict

from drawm.svg.svg_utils import render_label
from drawm.tree.newick_utils import parse_label


class LineageProps:
    """Visual attributed for named lineages."""

    def __init__(self, config_file, dwg, tree_height, inch, font_size):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.tree_height = tree_height
        self.inch = inch
        self.font_size = font_size
        
        self.node_radius = 0.005 * self.tree_height
        
        self.show_lineages = False
        self.display_method = None
        self.contour_width = 1
        self.font_size = None
        self.font_color = None
        self.lineage_map = {}

        with open(config_file) as f:
            prop_type = f.readline().strip()
            
            if prop_type != 'LINEAGES':
                self.logger.error("[LineageProps] Unexpected property type '%s' in %s." % (prop_type, config_file))
                sys.exit()
            
            for line in f:
                if line[0] == '#' or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                attribute, values = fields[0], fields[1:]
                if attribute == 'show_lineages':
                    self.show_lineages = (values[0] == 'True')
                elif attribute == 'display_method':
                    display_method = values[0]
                    assert(display_method in ['outlines', 'arc_labels'])
                    self.display_method = values[0]
                elif attribute == 'contour_width':
                    self.contour_width = float(values[0])
                elif attribute == 'font_size':
                    self.font_size = float(values[0])
                elif attribute == 'font_color':
                    self.font_color = values[0]
                elif attribute == 'lineage':
                    lineage_name, color, alpha, stroke_width = values
                    self.lineage_map[lineage_name] = (color, alpha, int(stroke_width))
       
    def render(self, tree):
        """Render lineages."""
        
        if not self.show_lineages:
            return
        
        if self.display_method == 'outlines':
            self._render_outlines(tree)
        elif self.display_method == 'arc_labels':
            self._render_arc_labels(tree)
            
    def _render_outlines(self, tree):
        """Color named lineages."""
        
        for node in tree.preorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)
            if taxon in self.lineage_map:
                color, alpha, stroke_width = self.lineage_map[taxon]
                
                path = self.dwg.path()
                path.fill(color=color, opacity=alpha)
                path.stroke(color=color, width=stroke_width)

                # start at current node
                path.push("M%d,%d" % (node.x, node.y))
                
                # descend 'right' of lineage
                right_branch = node.child_node_iter().next()
                while True:
                    # draw arc
                    angle_dir = '+'
                    if (right_branch.parent_node.angle - right_branch.angle) % 360 <= 180:
                        # child node is further clockwise than parent so must
                        # draw angle counter-clockwise (i.e., negative) direction
                        angle_dir = '-'
                        
                    path.push_arc(target=(right_branch.corner_x, right_branch.corner_y), 
                                rotation=0, 
                                r=right_branch.parent_node.rel_depth,
                                large_arc=False,
                                angle_dir=angle_dir,
                                absolute=True)
                    
                    if right_branch.is_leaf():
                        break
                     
                    path.push("L%d,%d" % (right_branch.x, right_branch.y))
                    right_branch = right_branch.child_node_iter().next()

                # move across children    
                last_leaf = None
                for leaf in node.leaf_iter():
                    path.push("L%d,%d" % (leaf.x, leaf.y))
                    last_leaf = leaf
                
                # ascend 'left' of lineage
                left_branch = last_leaf
                while left_branch != node:
                    path.push("L%d,%d" % (left_branch.corner_x, left_branch.corner_y))
                    
                    # draw arc
                    angle_dir = '+'
                    if (left_branch.angle - left_branch.parent_node.angle) % 360 <= 180:
                        # child node is further clockwise than parent so must
                        # draw angle counter-clockwise (i.e., negative) direction
                        angle_dir = '-'
                        
                    path.push_arc(target=(left_branch.parent_node.x, left_branch.parent_node.y), 
                                rotation=0, 
                                r=left_branch.parent_node.rel_depth,
                                large_arc=False,
                                angle_dir=angle_dir,
                                absolute=True)

                    left_branch = left_branch.parent_node

                self.dwg.add(path)
                
    def _render_arc_labels(self, tree):
        """Color named lineages."""
        
        rel_depth = 1.00
        label_depth = defaultdict(int)
        for node in tree.postorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)
            
            if node.is_leaf():
                label_depth[node.id] = 00
                continue
            
            max_child_label_depth = 0
            for c in node.child_node_iter():
                if label_depth[c.id] > max_child_label_depth:
                    max_child_label_depth = label_depth[c.id]
            
            if taxon not in self.lineage_map:
                label_depth[node.id] = max_child_label_depth
            else:
                # find deepest leaf node in lineage
                deepest_rel_depth = 0
                for leaf in node.leaf_iter():
                    if leaf.rel_depth > deepest_rel_depth:
                        deepest_rel_depth = leaf.rel_depth
                
                # draw arc
                color, alpha, stroke_width = self.lineage_map[taxon]
    
                label_depth[node.id] = max_child_label_depth+1
                
                leaves = [leaf for leaf in node.leaf_iter()]
                start_leaf = leaves[0]
                end_leaf = leaves[-1]
                
                # draw arc to parent
                angle_dir = '+'
                large_arc = False
                if (start_leaf.angle - end_leaf.angle) % 360 <= 180:
                    # start leaf is further clockwise than end leaf so must
                    # draw angle counter-clockwise (i.e., negative) direction
                    #angle_dir = '-'
                    large_arc=True
 
                depth = deepest_rel_depth * (1.0 + 0.05*label_depth[node.id])

                start_angle_rad = math.radians(start_leaf.angle)
                start_x = depth * math.cos(start_angle_rad) + 0.5*self.dwg.canvas_width
                start_y = depth * math.sin(start_angle_rad) + 0.5*self.dwg.canvas_height
                
                end_angle_rad = math.radians(end_leaf.angle)
                end_x = depth * math.cos(end_angle_rad) + 0.5*self.dwg.canvas_width
                end_y = depth * math.sin(end_angle_rad) + 0.5*self.dwg.canvas_height
                    
                p = self.dwg.path('m%d,%d' % (start_x, start_y))
                p.push_arc(target=(end_x, end_y), 
                            rotation=0, 
                            r=depth,
                            large_arc=large_arc,
                            angle_dir=angle_dir,
                            absolute=True)
                p.fill(color='none')
                p.stroke(color=color, width=stroke_width)
                self.dwg.add(p)
                
                label_angle = (0.5*(start_leaf.angle + end_leaf.angle)) % 360
                label_angle_rad = math.radians(label_angle)
                label_x = depth * math.cos(label_angle_rad) + 0.5*self.dwg.canvas_width + stroke_width
                label_y = depth * math.sin(label_angle_rad) + 0.5*self.dwg.canvas_height + stroke_width
                render_label(self.dwg, 
                                label_x, 
                                label_y, 
                                label_angle, 
                                taxon, 
                                self.font_size, 
                                self.font_color)
   