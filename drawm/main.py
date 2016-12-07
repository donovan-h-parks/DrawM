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

from biolib.common import check_file_exists, check_dir_exists
from biolib.logger import logger_setup

from drawm.draw_tree import DrawTree

from drawm.tree.reroot import Reroot
from drawm.tree.prune import Prune
from drawm.tree.subtree import Subtree
from drawm.tree.rename import Rename


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
        check_file_exists(options.config_file)
        
        draw_tree = DrawTree()
        draw_tree.render(options.input_tree,
                            options.config_file,
                            options.width,
                            options.height,
                            options.dpi,
                            options.output_prefix)
        
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
                    
    def rename(self, options):
        """Rename taxa."""
        
        check_file_exists(options.input_tree)
        
        rename = Rename()
        rename.run(options.input_tree,
                    options.taxa_label_map,
                    options.output_tree)
                    
    def convert(self, options):
        """Convert format of tree."""
        
        check_file_exists(options.input_tree)
        
        self.logger.info('Reading input tree.')
        tree = dendropy.Tree.get_from_path(options.input_tree, 
                                            schema=options.input_format,
                                            preserve_underscores=True)
        
        self.logger.info('Writing output tree.')        
        tree.write_to_path(options.output_tree, 
                            schema=options.output_format, 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
        
           
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
        elif(options.subparser_name == 'rename'):
            self.rename(options)
        elif(options.subparser_name == 'convert'):
            self.convert(options)
        else:
            self.logger.error('Unknown DrawM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
  
