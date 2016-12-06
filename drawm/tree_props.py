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
import random

from tree_utils import dist_to_ancestor
from newick_utils import parse_label

from biolib.taxonomy import Taxonomy

import dendropy


class TreeProps:
    """Visual attributed for displaying tree."""

    def __init__(self, config_file, dwg, inch):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')
        
        self.dwg = dwg
        self.inch = inch
        
        self.tree_x = 0.5 * self.dwg.canvas_width
        self.tree_y = 0.5 * self.dwg.canvas_height

        self.show_tree = False
        self.display_method = None
        self.ladderize = None
        self.branch_transformation = None
        
        self.width = None
        self.height = None
        
        self.rotation = None
        self.arc = None
        
        self.branch_width = None
        self.node_radius = None
        
        self.show_scale_bar = False
        self.scale_bar_width = None
        self.scale_font_size = None
        self.show_scale_bar_contours = False
        self.scale_bar_contour_width = None
        
        self.prune_by_taxon = None

        with open(config_file) as f:
            prop_type = f.readline().strip()
            
            if prop_type != 'TREE':
                self.logger.error("[TreeProps] Unexpected property type '%s' in %s." % (prop_type, config_file))
                sys.exit()
            
            for line in f:
                if line[0] == '#' or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                attribute, values = fields[0], fields[1:]
                if attribute == 'show_tree':
                    self.show_tree = (values[0] == 'True')
                elif attribute == 'display_method':
                    self.display_method = values[0]
                elif attribute == 'ladderize':
                    self.ladderize = values[0]
                elif attribute == 'branch_transformation':
                    self.branch_transformation = values[0]
                elif attribute == 'width':
                    self.width = float(values[0])*self.dwg.canvas_width
                elif attribute == 'height':
                    self.height = float(values[0])*self.dwg.canvas_height
                elif attribute == 'rotation':
                    self.rotation = float(values[0])
                elif attribute == 'arc':
                    self.arc = float(values[0])
                elif attribute == 'branch_width':
                    self.branch_width = int(values[0])
                elif attribute == 'node_radius':
                    self.node_radius = float(values[0])*self.width
                elif attribute == 'show_scale_bar':
                    self.show_scale_bar = (values[0] == 'True')
                elif attribute == 'scale_bar_width':
                    self.scale_bar_width = float(values[0])
                elif attribute == 'scale_font_size':
                    self.scale_font_size = float(values[0])
                elif attribute == 'show_scale_bar_contours':
                    self.show_scale_bar_contours = (values[0] == 'True')
                elif attribute == 'scale_bar_contour_width':
                    self.scale_bar_contour_width = float(values[0])
                elif attribute == 'prune_by_taxon':
                    self.prune_by_taxon = (values[0], int(values[1]))
                    
        # modify values for specific visual attributes
        if self.display_method == 'circular':
            self.width = 0.5*self.width
            self.height = self.width
            
    def _cladogram(self, tree):
        """Transform branch lengths to form a cladogram."""
        
        # get depth of all nodes
        tree.seed_node.depth = 0
        for n in tree.preorder_node_iter():
            if n is not tree.seed_node:
                n.depth = n.parent_node.depth+1
                
        # get deepest leaf node for each node
        for n in tree.preorder_node_iter():
            if n.is_leaf():
                n.deepest_leaf = 0
            else:
                deepest_leaf = 0
                for leaf in n.leaf_iter():
                    if leaf.depth > deepest_leaf:
                        deepest_leaf = leaf.depth
                        
                n.deepest_leaf = deepest_leaf - n.depth
        
        # rescale branches
        for n in tree.preorder_node_iter():
            if n is not tree.seed_node:
                n.edge.length = n.parent_node.deepest_leaf - n.deepest_leaf
                
    def _prune(self, tree):
        """Prune tree."""
        
        rank_label, num_genomes_to_retain = self.prune_by_taxon
        self.logger.info('Pruning tree to %d taxa at each %s.' % (num_genomes_to_retain, rank_label))
        
        taxa_to_retain = set()
        for node in tree.preorder_node_iter():
            _support, taxon, _auxiliary_info = parse_label(node.label)
            
            if taxon:
                most_specific_taxon = taxon.split(';')[-1].strip()
                most_specific_rank_prefix = Taxonomy.rank_prefixes.index(most_specific_taxon[0:3])
                most_specific_rank_label = Taxonomy.rank_labels[most_specific_rank_prefix]
                
                if most_specific_rank_label == rank_label:
                    taxa = [leaf.taxon for leaf in node.leaf_iter()]
                    selected = random.sample(taxa, min(len(taxa), num_genomes_to_retain))
                    taxa_to_retain.update(selected)
                    
        tree.retain_taxa(taxa_to_retain)
                    
        self.logger.info('Pruned tree to %d taxa.' % len(taxa_to_retain))

    def read_tree(self, input_tree_file):
        """Read tree from file."""
        
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree_file,
                                            schema='newick',
                                            rooting='force-rooted',
                                            preserve_underscores=True)
                                            
        # give each node an ID that is recoverable in different programs
        num_taxa = 0
        for tree_id, node in enumerate(tree.preorder_node_iter()):
            node.id = tree_id
            if node.is_leaf():
                num_taxa += 1
            
        self.logger.info('Tree contains %d taxa.' % num_taxa)
            
        # check if tree needs to be pruned
        if self.prune_by_taxon:
            self._prune(tree)
                    
        # ladderize tree as requested
        if self.ladderize == 'top':
            tree.ladderize(ascending=False)
        elif self.ladderize == 'bottom':
            tree.ladderize(ascending=True)
            
        # transform branches
        if self.branch_transformation == 'cladogram':
            self._cladogram(tree)
            
        return tree
                        
    def _circular_layout(self, tree):
        """Calculate position of nodes in circular tree layout."""
        
        self.logger.info('Performing circular layout.')
        
        # calculate number of leaf nodes below each node
        self.deepest_node = 0
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.num_leaves = 1
                dist_to_root = dist_to_ancestor(node, tree.seed_node)
                if dist_to_root > self.deepest_node :
                    self.deepest_node  = dist_to_root
            else:
                node.num_leaves = sum([c.num_leaves for c in node.child_node_iter()])
                
        self.logger.info('Total number of leaves: %d' % tree.seed_node.num_leaves)
        self.logger.info('Deepest leaf node: %.2f' % self.deepest_node)
        
        # calculate position of each node in x,y plane
        self.logger.info('Calculating position of each node.')
        cur_leaf_angle = self.rotation
        angle_step_size = self.arc / tree.seed_node.num_leaves
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                # leaf nodes are layout at equal angles around the rooted
                angle_rad = math.radians(cur_leaf_angle)
                angle = cur_leaf_angle
                cur_leaf_angle = (cur_leaf_angle + angle_step_size) % 360
            else:      
                # internal nodes are placed at angle between children
                angles = []
                for c in node.child_node_iter():
                    angles.append(c.angle)
                
                angle_diff = abs(angles[0] - angles[1])
                if angle_diff < 180:
                    angle = (sum(angles) / 2.0) % 360
                else:
                    angle = ((sum(angles) + 360.0) / 2.0) % 360
                     
                angle_rad = math.radians(angle)

            depth = dist_to_ancestor(node, tree.seed_node)
            rel_depth = (depth / self.deepest_node ) * self.width
                
            x = rel_depth * math.cos(angle_rad) + self.tree_x
            y = rel_depth * math.sin(angle_rad) + self.tree_y
            
            # save position information for node
            node.angle = angle
            node.x = x
            node.y = y
            node.rel_depth = rel_depth
            
        # layout corners
        for node in tree.postorder_node_iter():
            if node == tree.seed_node:
                node.corner_x = node.x - 1
                node.corner_y = node.y
                continue
                
            corner_angle_rad = math.radians(node.angle)
            corner_x = node.parent_node.rel_depth * math.cos(corner_angle_rad) + self.tree_x
            corner_y = node.parent_node.rel_depth * math.sin(corner_angle_rad) + self.tree_y
            
            node.corner_x = corner_x
            node.corner_y = corner_y
            
    def _rectangular_layout(self, tree):
        """Calculate position of nodes in rectangular tree layout."""
        
        self.logger.info('Performing rectangular layout.')

        # calculate number of leaf nodes below each node
        self.deepest_node = 0
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.num_leaves = 1
                dist_to_root = dist_to_ancestor(node, tree.seed_node)
                if dist_to_root > self.deepest_node :
                    self.deepest_node  = dist_to_root
            else:
                node.num_leaves = sum([c.num_leaves for c in node.child_node_iter()])
                
        self.logger.info('Total number of leaves: %d' % tree.seed_node.num_leaves)
        self.logger.info('Deepest leaf node: %.2f' % self.deepest_node)
        
        # calculate position of each node in x,y plane
        self.logger.info('Calculating position of each node.')
        border_x = 0.5*(self.dwg.canvas_width - self.width)
        border_y = 0.5*(self.dwg.canvas_height - self.height)
        y_step = float(self.height) / (tree.seed_node.num_leaves - 1.0)
        y_pos = border_y
        for node in tree.postorder_node_iter():
            depth = dist_to_ancestor(node, tree.seed_node)
            node.rel_depth = (depth / self.deepest_node ) * self.width    
            node.x = node.rel_depth + border_x
            
            if node.is_leaf():
                node.y = y_pos
                y_pos += y_step
            else:
                # put node at midpoint of its children
                y_children = []
                for c in node.child_node_iter():
                    y_children.append(c.y)
                node.y = np_mean(y_children)
                
        # layout corners  
        for node in tree.postorder_node_iter():
            if node == tree.seed_node:
                node.corner_x = node.x - 0.001*self.width
                node.corner_y = node.y
                continue
                
            node.corner_x = node.parent_node.x
            node.corner_y = node.y
            
    def layout(self, tree):
        """Layout tree."""
        
        if self.display_method == 'circular':
            self._circular_layout(tree)
        elif self.display_method == 'rectangular':
            self._rectangular_layout(tree)
            
    def _render_circular(self, tree, bootstrap_props):
        """Render circular tree."""
        
        # draw all nodes as small circles and branches as lines
        for node in tree.postorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)

            if node.parent_node:
                # draw line to corner leading to parent
                b = self.dwg.line(start=(node.x, node.y), 
                                    end=(node.corner_x, node.corner_y), 
                                    fill='black', 
                                    stroke_width=self.branch_width)
                b.stroke(color='black')
                self.dwg.add(b)
                
                # draw arc to parent
                angle_dir = '+'
                if (node.angle - node.parent_node.angle) % 360 <= 180:
                    # child node is further clockwise than parent so must
                    # draw angle counter-clockwise (i.e., negative) direction
                    angle_dir = '-'
                    
                p = self.dwg.path('m%d,%d' % (node.corner_x, node.corner_y))
                p.push_arc(target=(node.parent_node.x, node.parent_node.y), 
                            rotation=0, 
                            r=node.parent_node.rel_depth,
                            large_arc=False,
                            angle_dir=angle_dir,
                            absolute=True)
                p.fill(color='none')
                p.stroke(color='black')
                self.dwg.add(p)
            else:
                # take special care of root
                pass
                
            bootstrap_props.render(node.x, node.y, support)
            
    def _render_rectangular(self, tree, bootstrap_props):
        """Render rectangular tree."""
        
        # draw all nodes as small circles and branches as lines
        for node in tree.postorder_node_iter():
            support, taxon, auxiliary_info = parse_label(node.label)

            if node.parent_node:
                # draw line from node to corner
                b = self.dwg.line(start=(node.x, node.y), 
                                    end=(node.corner_x, node.corner_y), 
                                    fill='none')
                b.stroke(color='black', width=self.branch_width)
                self.dwg.add(b)
                
                # draw line from corner to parent
                b = self.dwg.line(start=(node.corner_x, node.corner_y), 
                                    end=(node.parent_node.x, node.parent_node.y), 
                                    fill='none')
                b.stroke(color='black', width=self.branch_width)
                self.dwg.add(b)
                
            else:
                # take special care of root
                pass
                
            bootstrap_props.render(node.x, node.y, support)

    def render(self, tree, bootstrap_props):
        """Render tree in x,y plane."""
        
        if not self.show_tree:
            return
            
        self.logger.info('Rendering tree.')
            
        if self.display_method == 'circular':
            self._render_circular(tree, bootstrap_props)
        elif self.display_method == 'rectangular':
            self._render_rectangular(tree, bootstrap_props)
            
    def _scale_width_bl(self, deepest_node):
        """Calculate width of scale in terms of branch length."""
        
        target_scale_bl = deepest_node / 10.0
        scalebar_width_bl = round(target_scale_bl, -int(math.floor(math.log10(abs(target_scale_bl)))))
        
        return scalebar_width_bl
            
    def render_scale_bar(self):
        """Render scale bar."""
        
        if not self.show_scale_bar:
            return
        
        scalebarX = 0.1*self.inch
        scalebarY = self.dwg.canvas_height-0.1*self.inch
        
        # scale bar should be ~10% of the longest branch,
        # but rounded to first significant figure
        scalebar_width_bl = self._scale_width_bl(self.deepest_node)
        scalebar_width = (scalebar_width_bl / self.deepest_node) * self.height

        # draw scale bar
        bar = self.dwg.line(start=(scalebarX, scalebarY), 
                                    end=(scalebarX+scalebar_width, scalebarY), 
                                    fill='black', 
                                    stroke_width=self.scale_bar_width)
        bar.stroke(color='black')
        self.dwg.add(bar)
        
        left_tick = self.dwg.line(start=(scalebarX, scalebarY-0.05*scalebar_width), 
                                    end=(scalebarX, scalebarY+0.05*scalebar_width), 
                                    fill='black', 
                                    stroke_width=self.scale_bar_width)
        left_tick.stroke(color='black')
        self.dwg.add(left_tick)
        
        right_tick = self.dwg.line(start=(scalebarX+scalebar_width, scalebarY-0.05*scalebar_width), 
                                    end=(scalebarX+scalebar_width, scalebarY+0.05*scalebar_width), 
                                    fill='black', 
                                    stroke_width=self.scale_bar_width)
        right_tick.stroke(color='black')
        self.dwg.add(right_tick)
        
        # draw scale bar label
        if scalebar_width_bl > 1:
            scalebar_width_label = '%d' % scalebar_width_bl
        else:
            scalebar_width_label = '%.1g' % scalebar_width_bl
        
        label = self.dwg.text(scalebar_width_label, 
                                x=[(scalebarX+0.5*scalebar_width)], 
                                y=[(scalebarY-0.05*scalebar_width)], 
                                font_size=self.scale_font_size,
                                text_anchor="middle",
                                fill='black')
        self.dwg.add(label)
        
    def render_scale_lines(self):
        """Render scale bar."""
        
        if not self.show_scale_bar_contours:
            return

        colors = [color_str(64, 64, 64), 'blue']
        
        # scale bar should be ~10% of the longest branch,
        # but rounded to first significant figure
        scalebar_width_bl = self._scale_width_bl(self.deepest_node)
        scalebar_width = (scalebar_width_bl / self.deepest_node) * self.height

        radius = scalebar_width
        max_radius = self.height
        color_index = 0
        while radius < max_radius:
            c = self.dwg.circle(center=(self.tree_x, self.tree_y), r=radius)
            c.fill('none')
            c.stroke(color=colors[color_index], width=self.scale_bar_contour_width)
            #c.dasharray([0.01*self.inch, 0.01*self.inch])
            self.dwg.add(c)  
            
            radius += scalebar_width
            color_index = (color_index + 1) % 2