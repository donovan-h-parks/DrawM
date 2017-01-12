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

from drawm.svg.svg_utils import donut


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
                elif attribute == 'contour_method':
                    self.contour_method = values[0]
                    assert (self.contour_method in ['CONCENTRIC', 'BY_FILE'])
                elif self.contour_method == 'BY_FILE' and attribute == 'contour_file':
                    self.contour_file = values[0]
                elif attribute == 'contour_width':
                    self.contour_width = float(values[0])
                elif attribute == 'contour_cm':
                    outer_threshold, inner_threshold, color, alpha, label = values
                    self.contour_cm.append((float(outer_threshold), float(inner_threshold), color, float(alpha), label))
                else:
                    self.logger.warning('[ContourProps] Unexpected attribute: %s' % attribute)
                    
    def decorate(self, tree):
        """Decorate tree with contour information."""
        
        if not self.show_contours:
            return
            
        if self.contour_method != 'BY_FILE':
            return
        
        contours = {}
        for line in open(self.contour_file):
            if line[0] == '#':
                continue

            line_split = line.strip().split('\t')
            taxa = line_split[0]
            contour_value = float(line_split[1])
            
            if '|' in taxa:
                node = tree.mrca(taxon_labels=taxa.split('|'))
            else:
                # specified leaf node
                node = tree.find_node_with_taxon_label(taxa)
                            
            node.contour = contour_value
        
    def render_legend(self):
        """Render legend."""
        
        if not self.show_contours:
            return
            
        legend_radius = 2*self.node_radius
        legend_step = 2*2*legend_radius
        
        legendX = 0.1*self.inch
        legendY = 0.1*self.inch
        
        legend_group = svgwrite.container.Group(id='contour_legend')
        self.dwg.add(legend_group)

        for item_index, (outer_threshold, inner_threshold, color, alpha, label) in enumerate(self.contour_cm):
            c = self.dwg.circle(center=(legendX, legendY), 
                                r=legend_radius,
                                id='contour_legend_symbol_%d' % item_index)
            c.fill(color=color, opacity=alpha)
            c.stroke(color='black')
            legend_group.add(c)

            t = self.dwg.text('%s: %g to %g' % (label, inner_threshold, outer_threshold), 
                                x=[(legendX + 1.5*legend_radius)], 
                                y=[(legendY + 0.35*self.font_size)], 
                                font_size=self.font_size,
                                fill='black',
                                id='contour_legend_label_%d' % item_index)
            legend_group.add(t)
                
            legendY += legend_step
            
    def _contour_pts(self, tree, contour_threshold, draw_threshold=None):
        """Get points defining contour."""
        
        if not draw_threshold:
            draw_threshold = contour_threshold
            
        # determine direction of contour
        descending = tree.seed_node.contour > tree.leaf_node_iter().next().contour
        
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
            if (node == tree.seed_node 
                or (descending and node.contour > contour_threshold)
                or (not descending and node.contour < contour_threshold)):
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
        
    def _contour_by_file(self, tree):
        """Draw contour based on node information in file."""
        
        contour_group = svgwrite.container.Group(id='contour')
        self.dwg.add(contour_group)
        
        for index, (outer_threshold, inner_threshold, color, alpha, label) in enumerate(self.contour_cm): 
            outer_pts, outer_nodes = self._contour_pts(tree, outer_threshold)
            cheaty_inner_pts, cheaty_inner_nodes = self._contour_pts(tree, outer_threshold, inner_threshold)
            inner_pts, inner_nodes = self._contour_pts(tree, inner_threshold)

            # draw outer contour
            path = self.dwg.path("M%d,%d" % outer_pts[0], id='contour_%d' % index)
            path.fill(color=color, opacity=alpha, rule='evenodd')
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
                
            contour_group.add(path)
        
    def _contour_concentric(self, tree):
        """Draw concentric contours using branch length."""
        
        contour_group = svgwrite.container.Group(id='contour')
        self.dwg.add(contour_group)
        
        for index, (outer_threshold, inner_threshold, color, alpha, label) in enumerate(self.contour_cm):
            donut(self.dwg, 
                    tree.seed_node.x, tree.seed_node.x, 
                    inner_threshold*self.tree_height, outer_threshold*self.tree_height, 
                    color, opacity=alpha, 
                    group=contour_group, id='contour_%d' % index)
            if False:
                c = self.dwg.circle(center=(tree.seed_node.x, tree.seed_node.x), 
                                    r=outer_threshold*self.tree_height,
                                    id='contour_%d' % index)
                c.fill('none')
                c.stroke(color=color, 
                            width=self.contour_width)
                contour_group.add(c)  
        
    def render_contour(self, tree):
        """Render contour."""
        
        if not self.show_contours:
            return
            
        contour_group = svgwrite.container.Group(id='contour')
        self.dwg.add(contour_group)
        
        if self.contour_method == 'BY_FILE':
            self._contour_by_file(tree)
        elif self.contour_method == 'CONCENTRIC':
            self._contour_concentric(tree)

