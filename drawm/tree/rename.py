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

import sys
import logging

import dendropy

from newick_utils import parse_label


class Rename(object):
    """Rename taxa in tree."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

    def run(self, input_tree, taxa_label_map, output_tree):
        """Rename taxa in tree.
        
        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        taxa_label_map : str
          File specifying taxa to rename.
        output_tree : str
          Name of file for rerooted tree.
        """
        
        self.logger.info('Reading input tree.')
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        
        # read mapping
        self.logger.info('Reading taxa map.')
        taxon_map = {}
        for line in open(taxa_label_map):
            if line[0] == '#':
                continue
                
            line_split = line.strip().split('\t')
            if len(line_split) == 2:
                taxon_map[line_split[0]] = line_split[1]
                
        # rename taxa
        modified_taxa = 0
        for node in tree.preorder_node_iter():
            if node.is_leaf():
                _support, taxon, _aux_info = parse_label(node.taxon.label)
            else:
                _support, taxon, _aux_info = parse_label(node.label)
                
            new_taxon = taxon_map.get(taxon, None)
            if new_taxon:
                modified_taxa += 1
                taxon_map.pop(taxon)
                
                if node.is_leaf():
                    node.taxon.label = node.taxon.label.replace(taxon, new_taxon)
                else:
                    node.label = node.label.replace(taxon, new_taxon)
                    
        self.logger.info('Replaced %d taxon labels.' % modified_taxa)
        if len(taxon_map) != 0:
            self.logger.warning('Not all taxa were identified in tree: %s' % ','.join(taxon_map))
  
        # write out tree
        self.logger.info('Writing output tree.')  
        tree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
                    