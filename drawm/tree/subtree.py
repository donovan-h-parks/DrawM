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


class Subtree(object):
    """Extract subtree from tree."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

    def run(self, input_tree, output_tree, subtree_taxa):
        """Reroot tree.
        
        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        output_tree : str
          Name of file for rerooted tree.
        subtree_taxa : list
          Taxa specifying subtree to extract.
        """
        
        self.logger.info('Reading input tree.')
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
                                            
        # find taxa among nodes
        subtree_taxa_in_tree = set()
        subtree_taxa = set(subtree_taxa)
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                support, taxon, _auxiliary_info = parse_label(node.taxon.label)
                if taxon in subtree_taxa:
                    subtree_taxa_in_tree.add(node.taxon)
                    subtree_taxa.remove(taxon)
            else:
                support, taxon, _auxiliary_info = parse_label(node.label)
                if taxon in subtree_taxa:
                    for leaf in node.leaf_iter():
                        subtree_taxa_in_tree.add(leaf.taxon)
                    subtree_taxa.remove(taxon)
                        
            # check if all outgroup taxa have been identified          
            if not subtree_taxa:
                break
                
        # identify MRCA of taxa in subtree
        mrca = tree.mrca(taxa=subtree_taxa_in_tree)
        mrca_leaves = [leaf.taxon for leaf in mrca.leaf_iter()]
        self.logger.info('Identified %d taxa in subtree.' % len(mrca_leaves))
        subtree = tree.extract_tree_with_taxa(mrca_leaves)
  
        # write out results
        self.logger.info('Writing output tree.')  
        subtree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
                    