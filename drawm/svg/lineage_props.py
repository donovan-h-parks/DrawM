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

import svgwrite

from drawm.svg.geometry import unit_vector
from drawm.svg.svg_utils import render_label
from drawm.tree.newick_utils import parse_label
from drawm.tree.tree_utils import find_node


class LineageProps:
    """Visual attributed for named lineages."""

    def __init__(self, config_file, dwg, inch):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.inch = inch

        self.show_lineages = False
        self.display_method = None
        self.display_depth = 'MAX'
        self.font_size = None
        self.font_color = None
        self.lineage_map = {}
        
        if not config_file:
            return # use default values

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
                    assert(display_method in ['OUTLINE_LINEAGE', 'ARC_LABELS'])
                    self.display_method = values[0]
                elif attribute == 'display_depth':
                    display_depth = values[0]
                    assert(display_depth in ['TIGHT', 'MAX'])
                    self.display_depth = values[0]    
                elif attribute == 'font_size':
                    self.font_size = float(values[0]) * (self.inch/90.0)
                elif attribute == 'font_color':
                    self.font_color = values[0]
                elif attribute == 'lineage':
                    lineage_name, lineage_label, color, alpha, stroke_width = values
                    self.lineage_map[lineage_name] = (lineage_name, lineage_label, color, alpha, int(stroke_width))
                else:
                    self.logger.warning('[LineageProps] Unexpected attribute: %s' % attribute)
       
    def render(self, tree):
        """Render lineages."""
        
        if not self.show_lineages:
            return
            
        self.logger.info('Rendering lineages.')
            
        self._lineage_nodes(tree)
        
        lineage_group = svgwrite.container.Group(id='lineage')
        self.dwg.add(lineage_group)
        
        lineage_text_group = svgwrite.container.Group(id='lineage_text')
        self.dwg.add(lineage_text_group)
        
        if self.display_method == 'OUTLINE_LINEAGE':
            self._render_outlines(tree, lineage_group, lineage_text_group)
        elif self.display_method == 'ARC_LABELS':
            self._render_arc_labels(tree, lineage_group, lineage_text_group)
            
    def _lineage_nodes(self, tree):
        """Identify nodes in tree associated with lineages to render."""
        
        new_lineage_map = {}
        for lineage_name, data in self.lineage_map.iteritems():
            node = find_node(tree, lineage_name)

            if node:
                new_lineage_map[node] = data
            else:
                self.logger.warning('Failed to identify node with label: %s.' % lineage_name)
            
        self.lineage_map = new_lineage_map
            
    def _outline_circular(self, 
                            node, 
                            taxon, 
                            color, 
                            alpha, 
                            stroke_width, 
                            lineage_group, 
                            lineage_text_group):
        """Outline lineage in circular tree."""
        
        path = self.dwg.path(id='lineage_%s' % taxon.replace(' ', '_'))
        path.fill(color=color, opacity=alpha)
        path.stroke(color=color, width=stroke_width)

        # start at current node
        path.push("M%f,%f" % (node.x, node.y))
        
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
             
            path.push("L%f,%f" % (right_branch.x, right_branch.y))
            right_branch = right_branch.child_node_iter().next()

        # move across children    
        last_leaf = None
        for leaf in node.leaf_iter():
            path.push("L%f,%f" % (leaf.x, leaf.y))
            last_leaf = leaf
        
        # ascend 'left' of lineage
        left_branch = last_leaf
        while left_branch != node:
            path.push("L%f,%f" % (left_branch.corner_x, left_branch.corner_y))
            
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

        lineage_group.add(path)
                
    def _outline_rectangular(self,
                                tree,
                                node, 
                                taxon, 
                                color, 
                                alpha, 
                                stroke_width, 
                                lineage_group, 
                                lineage_text_group):
        """Outline lineage in circular tree."""
        
        # get top and bottom leaf nodes
        leaves = [leaf for leaf in node.postorder_iter(lambda n: n.is_leaf())]
        start_leaf = leaves[0]
        end_leaf = leaves[-1]

        # determine start and end of rectangle covering lineage
        start_x = node.x
        end_x = 0
        if self.display_depth == 'MAX':
            end_x = tree.width + tree.start_x
        elif self.display_depth == 'TIGHT':
            for leaf in node.leaf_iter():
                if leaf.x > end_x:
                    end_x = leaf.x
                
        start_y = start_leaf.y
        end_y = end_leaf.y
        
        rect = self.dwg.rect(insert=(start_x, start_y),
                                size=(abs(end_x-start_x), abs(end_y-start_y)),
                                id='lineage_%s' % taxon.replace(' ', '_'))
        rect.fill(color=color, opacity=alpha)
        rect.stroke(color=color, width=stroke_width)
        
        lineage_group.add(rect)

        # render label
        label_x = end_x + 0.05*self.inch
        label_y = 0.5*(start_y + end_y)

        render_label(self.dwg, 
                        label_x, 
                        label_y, 
                        0, 
                        taxon, 
                        self.font_size, 
                        self.font_color,
                        middle_y=True,
                        group=lineage_text_group,
                        id_prefix='lineage_text')
            
    def _render_outlines(self, tree, lineage_group, lineage_text_group):
        """Color named lineages."""
                
        for node in tree.preorder_node_iter(lambda n: not n.is_leaf()):
            support, taxon, auxiliary_info = parse_label(node.label)
            if node in self.lineage_map:
                # make sure lineage isn't collapsed
                if node.is_collapsed or node.is_collapsed_root:
                    continue
                    
                lineage_name, lineage_label, color, alpha, stroke_width = self.lineage_map[node]
                
                if tree.display_method == 'CIRCULAR':
                    self._outline_circular(node, 
                                            lineage_label, 
                                            color, 
                                            alpha, 
                                            stroke_width, 
                                            lineage_group, 
                                            lineage_text_group)
                elif tree.display_method == 'RECTANGULAR':
                    self._outline_rectangular(tree,
                                                node, 
                                                lineage_label, 
                                                color, 
                                                alpha, 
                                                stroke_width, 
                                                lineage_group, 
                                                lineage_text_group)
                
    def _render_arc_labels(self, tree, lineage_group, lineage_text_group):
        """Color named lineages."""

        rel_depth = 1.00
        label_depth = defaultdict(int)
        labels = []
        for node in tree.postorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)
            
            if node.is_leaf():
                label_depth[node.id] = 0
                continue
            
            max_child_label_depth = 0
            for c in node.child_node_iter():
                if label_depth[c.id] > max_child_label_depth:
                    max_child_label_depth = label_depth[c.id]
            
            if node not in self.lineage_map:
                label_depth[node.id] = max_child_label_depth
            else:
                # make sure lineage isn't collapsed
                if node.is_collapsed or node.is_collapsed_root:
                    continue
                    
                # find deepest leaf node in lineage
                deepest_rel_depth = 0
                deepest_x = 0
                for leaf in node.leaf_iter():
                    if leaf.rel_depth > deepest_rel_depth:
                        deepest_rel_depth = leaf.rel_depth
                        deepest_x = leaf.x
                
                # draw arc
                lineage_name, lineage_label, color, alpha, stroke_width = self.lineage_map[node]
    
                label_depth[node.id] = max_child_label_depth+1
                
                leaves = [leaf for leaf in node.postorder_iter(lambda n: n.is_leaf())]
                start_leaf = leaves[0]
                end_leaf = leaves[-1]
                
                if tree.display_method == 'CIRCULAR':
                    # draw arc to parent
                    angle_dir = '+'
                    large_arc = False
                    if (start_leaf.angle - end_leaf.angle) % 360 <= 180:
                        # start leaf is further clockwise than end leaf so must
                        # draw angle counter-clockwise (i.e., negative) direction
                        #angle_dir = '-'
                        large_arc=True
     
                    if self.display_depth == 'MAX':
                        depth = tree.height
                    elif self.display_depth == 'TIGHT':
                        depth = deepest_rel_depth
                    
                    depth += 0.05*self.inch*label_depth[node.id]

                    start_x = depth * start_leaf.x_dir + 0.5*self.dwg.canvas_width
                    start_y = depth * start_leaf.y_dir + 0.5*self.dwg.canvas_height

                    end_x = depth * end_leaf.x_dir + 0.5*self.dwg.canvas_width
                    end_y = depth * end_leaf.y_dir + 0.5*self.dwg.canvas_height
 
                    p = self.dwg.path('M%f,%f' % (start_x, start_y), 
                                        id='lineage_%s' % lineage_label.replace(' ', '_'))
                    p.push_arc(target=(end_x, end_y), 
                                rotation=0, 
                                r=depth,
                                large_arc=large_arc,
                                angle_dir=angle_dir,
                                absolute=True)
   
                    x = 0.5 * (start_leaf.x_dir + end_leaf.x_dir)
                    y = 0.5 * (start_leaf.y_dir + end_leaf.y_dir)
                    x_dir, y_dir = unit_vector((x,y))
                    if large_arc:
                        y_dir *= -1
                        x_dir *= -1
                    
                    label_x = (depth + 0.5*stroke_width + 0.02*self.inch) * x_dir + 0.5*self.dwg.canvas_width
                    label_y = (depth + 0.5*stroke_width + 0.02*self.inch) * y_dir + 0.5*self.dwg.canvas_height
                    label_angle = math.degrees(math.atan2(y_dir, x_dir))
                elif tree.display_method == 'RECTANGULAR':
                    if self.display_depth == 'MAX':
                        depth = tree.width + tree.start_x
                    elif self.display_depth == 'TIGHT':
                        depth = deepest_x
                    
                    depth += 0.05*self.inch*label_depth[node.id]
                    
                    p = self.dwg.line(start=(depth, start_leaf.y), 
                                        end=(depth, end_leaf.y),
                                        id='lineage_%s' % lineage_label.replace(' ', '_'))
                                        
                    label_x = depth + 0.5*stroke_width
                    label_y = 0.5*(start_leaf.y + end_leaf.y)
                    label_angle = 0
                    
                p.fill(color='none')
                p.stroke(color=color, opacity=alpha, width=stroke_width)
                lineage_group.add(p)
                
                labels.append((lineage_label, label_x, label_y, label_angle))
            
        # render arc labels after arcs to ensure text is on top
        for lineage_label, label_x, label_y, label_angle in labels:
            render_label(self.dwg, 
                            label_x, 
                            label_y, 
                            label_angle, 
                            lineage_label, 
                            self.font_size, 
                            self.font_color,
                            middle_y=True,
                            group=lineage_text_group,
                            id_prefix='lineage_text')
       