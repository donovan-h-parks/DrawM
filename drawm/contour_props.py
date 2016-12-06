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


class ContourProps:
    """Visual attributed for contours."""

    def __init__(self, config_file, dwg, tree_height, inch, font_size):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.tree_height = tree_height
        self.inch = inch
        self.font_size = font_size
        
        self.node_radius = 0.005 * self.tree_height
        
        self.show_contours = False
        self.contour_width = None
        self.contour_file = None
        self.contour_index = None
        self.contour_leaf_value = None
        self.contour_cm = []
        
        with open(config_file) as f:
            prop_type = f.readline().strip()
            
            if prop_type != 'CONTOURS':
                self.logger.error("[ContourProps] Unexpected property type '%s' in %s." % (prop_type, config_file))
                sys.exit()
            
            for line in f:
                if line[0] == '#' or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                attribute, values = fields[0], fields[1:]
                if attribute == 'show_contours':
                    self.show_contours = (values[0] == 'True')
                elif attribute == 'contour_file':
                    self.contour_file = values[0]
                elif attribute == 'contour_index':
                    self.contour_index = int(values[0])
                elif attribute == 'contour_width':
                    self.contour_width = float(values[0])
                elif attribute == 'contour_leaf_value':
                    self.contour_leaf_value = float(values[0])
                elif attribute == 'contour_cm':
                    outer_threshold, inner_threshold, color = values
                    self.contour_cm.append((float(outer_threshold), float(inner_threshold), color))
                    
    def decorate(self, tree):
        """Decorate tree with contour information."""
        
        if not self.show_contours:
            return
        
        contours = {}
        for line in open(self.contour_file):
            if line[0] == '#':
                continue

            line_split = line.strip().split('\t')
            tree_id = line_split[0]
            print 'tree_id', tree_id
            value = line_split[self.contour_index]
            num_leaves = int(line_split[-1])
            if value != 'none':
                contours[int(tree_id)] = (float(value), num_leaves)
            else:
                contours[int(tree_id)] = (-1, -1)
            
        for node in tree.preorder_node_iter():
            if node.is_leaf():
                node.contour = self.contour_leaf_value
            else:
                value, num_leaves = contours[node.id] 
                
                if num_leaves != -1:   
                    node.contour = value
                else:
                    node.contour = None
        
    def render_legend(self):
        """Render legend."""
        
        if not self.show_contours:
            return
            
        legend_radius = 2*self.node_radius
        legend_step = 2*2*legend_radius
        
        legendX = 0.1*self.inch
        legendY = 0.1*self.inch

        for outer_threshold, inner_threshold, color in self.contour_cm:
            c = self.dwg.circle(center=(legendX, legendY), r=legend_radius)
            c.fill(color=color)
            c.stroke(color='black')
            self.dwg.add(c)

            t = self.dwg.text('%d to %d%%' % (inner_threshold, outer_threshold), 
                                x=[(legendX + 1.5*legend_radius)], 
                                y=[(legendY + 0.35*self.font_size)], 
                                font_size=self.font_size,
                                fill='black')
            self.dwg.add(t)
                
            legendY += legend_step
            
    def _contour_pts(self, tree, contour_threshold, draw_threshold=None):
        """Get points defining contour."""
        
        if not draw_threshold:
            draw_threshold = contour_threshold
        
        pts = []
        nodes = []
        stack = [tree.seed_node]
        while stack:
            node = stack.pop()
            
            if node.contour == None:
                # nodes below this will also have undefined
                # contour values
                continue
             
            # check if node meets contouring criterion
            if node.contour < contour_threshold or node == tree.seed_node:
                # node is below threshold so add children
                for c in node.child_node_iter():
                    stack.append(c)
                continue
            else:
                # first node in lineage above threhold so
                # draw contour along this branch
                bl = node.contour
                pbl = node.parent_node.contour
                
                index = 0
                if bl != pbl:
                    index = float(bl - draw_threshold) / (bl - pbl)

                x = node.x + index * (node.corner_x - node.x)
                y = node.y + index * (node.corner_y - node.y)
                    
                pts.append((x, y))
                nodes.append(node)
                
        return pts, nodes
        
    def render_contour(self, tree):
        """Render contour."""
        
        if not self.show_contours:
            return
        
        for outer_threshold, inner_threshold, color in self.contour_cm: 
            outer_pts, outer_nodes = self._contour_pts(tree, outer_threshold)
            cheaty_inner_pts, cheaty_inner_nodes = self._contour_pts(tree, outer_threshold, inner_threshold)
            inner_pts, inner_nodes = self._contour_pts(tree, inner_threshold)

            # draw outer contour
            path = self.dwg.path("M%d,%d" % outer_pts[0])
            path.fill(color=color, opacity=0.5, rule='evenodd')
            path.stroke(color=color, width=self.contour_width)

            for pt in outer_pts[1:]:
                path.push("L%d,%d" % pt)
                   
            # draw inner ring of contour
            path.push("L%d,%d" % inner_pts[-1])
            for pt in inner_pts[::-1]:
                path.push("L%d,%d" % pt)
                
            # this improves the visual quality for trees that
            # only have deep nodes near the end of the tree
            angle_rad = math.radians(outer_nodes[0].angle)
            visual_x = inner_nodes[0].rel_depth * math.cos(angle_rad) + 0.5*self.dwg.canvas_width
            visual_y = inner_nodes[0].rel_depth * math.sin(angle_rad) + 0.5*self.dwg.canvas_height
            path.push("L%d,%d" % (visual_x,visual_y))

            # connect the inner and outer contours
            path.push("L%d,%d" % outer_pts[0])
                
            self.dwg.add(path)
