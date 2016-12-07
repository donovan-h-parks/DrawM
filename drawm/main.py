#!/usr/bin/env python

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

import os
import sys
import math
import logging
import argparse

import dendropy
import svgwrite

from biolib.common import check_file_exists, check_dir_exists
from biolib.logger import logger_setup

from drawm.svg.tree_props import TreeProps
from drawm.svg.label_props import LabelProps
from drawm.svg.bootstrap_props import BootstrapProps
from drawm.svg.lineage_props import LineageProps
from drawm.svg.contour_props import ContourProps

from drawm.tree.reroot import Reroot
from drawm.tree.prune import Prune
from drawm.tree.subtree import Subtree


"""
To Do:

- branch transform: equal
- node order: increasing or decreasing

- show bootstaps as dots
-- color map (need legend)
-- show actual values

- collapse lineages
-- list in file or config file

- labels for extent taxa
-- align at tips or at same depth (with optional dashes)
- coloring of named lineages
- option to not use branch lenths (i.e., cladogram)
- branch thickness, color, line style
- label font, style, size, etc.

- rectangular phylogram
-- with optional leaf sorting (branches with fewest leaves 'on top': ladder-left, ladder-right)

"""

class OptionsParser():
    """Parse input commands."""
    
    def __init__(self):
        """Initialization"""
        
        self.logger = logging.getLogger('timestamp')
        self.reporter = logging.getLogger('no_timestamp')
          
    def draw(self, options):
        """Create SVG image of tree in Newick format."""
        
        check_file_exists(options.input_tree)
        check_dir_exists(options.config_dir)
        
        tree_name = os.path.split(os.path.basename(options.input_tree))[0]

        # setup SVG file
        canvas_width = options.width * options.dpi
        canvas_height = options.height * options.dpi
        font_size = int(8 * (float(dpi)/90) + 0.5)
        
        self.logger.info('Setting up SVG file.')
        svg_output = output_prefix + '.svg'
        dwg = svgwrite.Drawing(filename=svg_output, 
                                    size=(canvas_width, canvas_height))
        dwg.set_desc(title='DrawM rendering of %d' % tree_name, desc=tree_name)
        dwg.canvas_width = canvas_width
        dwg.canvas_height = canvas_height
        
        background = "rgb(" + str(255) + "," + str(255) + "," + str(255) + ")"
        self.dwg.add(dwg.rect(insert=(0, 0), 
                        size=('100%', '100%'), 
                        rx=None, ry=None, 
                        fill=background, 
                        stroke='none'))
                        
        # read configuration files
        self.logger.info('Reading configuration files.')
        tree_props = TreeProps(os.path.join(config_dir, 'tree.cfg'),
                                dwg,
                                options.dpi)
                                            
        bootstrap_props = BootstrapProps(os.path.join(config_dir, 'bootstrap.cfg'),
                                            dwg,
                                            tree_props.height,
                                            options.dpi,
                                            font_size)
                                            
        contour_props = ContourProps(os.path.join(config_dir, 'aai_contours.cfg'),
                                        dwg,
                                        tree_props.height,
                                        options.dpi,
                                        font_size)
                                            
        lineage_props = LineageProps(os.path.join(config_dir, 'lineages.cfg'),
                                        dwg,
                                        tree_props.height,
                                        options.dpi,
                                        font_size)
                                            
        label_props = LabelProps(os.path.join(config_dir, 'labels.cfg'),
                                        dwg,
                                        tree_props.height,
                                        options.dpi,
                                        font_size)
        
        # read and ladderize tree
        tree = tree_props.read_tree(input_tree)
        tree_props.layout(tree)
                                           
        contour_props.decorate(tree)
                                                    
        bootstrap_props.render_legend()
        
        rd = RelativeDistance()
        rd.decorate_rel_dist(tree)
        
        phylorank_props.render_contour(tree)
        phylorank_props.render_legend()
        
        contour_props.render_contour(tree)
        contour_props.render_legend()
        
        lineage_props.render(tree)
   
        tree_props.render(tree, bootstrap_props)
        tree_props.render_scale_bar()
        tree_props.render_scale_lines()

        #decorate_mean_branch_length(tree)
        #self._render_mean_tip_bl(tree, 0.7)
        #self._render_mean_tip_bl_test(tree, 0.7)
        
        label_props.render(tree)
        
        self.logger.info('Saving SVG image.')
        dwg.save()
        
        self.logger.info('Saving PNG image.')
        png_output = output_prefix + '.png'
        os.system('"C:\Program Files\Inkscape\inkscape.exe" -e %s -d %d -z -w %d -h %d %s' % (png_output,
                                                                                                self.dpi,
                                                                                                draw_tree.canvas_width, 
                                                                                                draw_tree.canvas_height, 
                                                                                                svg_output))
    def reroot(self, options):
        """Reroot tree."""
        
        check_file_exists(options.input_tree)
        
        if options.midpoint and options.outgroup:
            self.logger.error("The 'midpoint' and 'outgroup' arguments cannot be used together.")
            sys.exit()
        elif not options.midpoint and not options.outgroup:
            self.logger.error("Either the 'midpoint' or 'outgroup' argument must be specified.")
            sys.exit()
            
        reroot = Reroot()
        reroot.run(options.input_tree, options.midpoint, options.outgroup, options.output_tree)
        
    def prune(self, options):
        """Prune tree."""
        
        check_file_exists(options.input_tree)
        check_file_exists(options.taxa_to_retain)
        
        prune = Prune()
        prune.run(options.input_tree,
                    options.taxa_to_retain,
                    options.output_tree)
        
    def subtree(self, options):
        """Extract subtree."""
        
        check_file_exists(options.input_tree)
        
        subtree = Subtree()
        subtree.run(options.input_tree,
                    options.output_tree,
                    options.subtree_taxa)
           
    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if(options.subparser_name == 'draw'):
            self.draw(options)
        elif(options.subparser_name == 'reroot'):
            self.reroot(options)
        elif(options.subparser_name == 'prune'):
            self.prune(options)
        elif(options.subparser_name == 'subtree'):
            self.subtree(options)
  
        else:
            self.logger.error('Unknown DrawM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
  
