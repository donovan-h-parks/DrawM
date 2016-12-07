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


class Reroot(object):
    """Reroot tree."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

    def root_with_outgroup(self, tree, outgroup):
        """Reroot the tree using the given outgroup.

        Parameters
        ----------
        tree : Tree
          Tree to reroot. 
        outgroup : iterable
          Labels of taxa in outgroup.
        """

        # find taxa among nodes
        outgroup_in_tree = set()
        outgroup = set(outgroup)
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                support, taxon, _auxiliary_info = parse_label(node.taxon.label)
                if taxon in outgroup:
                    outgroup_in_tree.add(node.taxon)
                    outgroup.remove(taxon)
            else:
                support, taxon, _auxiliary_info = parse_label(node.label)
                if taxon in outgroup:
                    for leaf in node.leaf_iter():
                        outgroup_in_tree.add(leaf.taxon)
                    outgroup.remove(taxon)
                        
            # check if all outgroup taxa have been identified          
            if not outgroup:
                break

        self.logger.info('Identified %d outgroup taxa in the tree.' % len(outgroup_in_tree))

        if len(outgroup_in_tree) == 0:
            self.logger.warning('No outgroup taxa identified in the tree.')
            self.logger.warning('Tree was not rerooted.')
            sys.exit()

        mrca = tree.mrca(taxa=outgroup_in_tree)

        if len(mrca.leaf_nodes()) != len(outgroup_in_tree):
            self.logger.info('Outgroup is not monophyletic. Tree will be rerooted at the MRCA of the outgroup.')
            self.logger.info('The outgroup consisted of %d taxa, while the MRCA has %d leaf nodes.' % (len(outgroup_in_tree), len(mrca.leaf_nodes())))
            if len(mrca.leaf_nodes()) == len(tree.leaf_nodes()):
                self.logger.warning('The MRCA spans all taxa in the tree.')
                self.logger.warning('This indicating the selected outgroup is likely polyphyletic in the current tree.')
                self.logger.warning('Polyphyletic outgroups are not suitable for rooting. Try another outgroup.')
        else:
            self.logger.info('Outgroup is monophyletic.')

        if mrca.edge_length is None:
            self.logger.info('Tree is already rooted on this outgroup.')
        else:
            self.logger.info('Rerooting tree.')
            tree.reroot_at_edge(mrca.edge,
                                length1=0.5 * mrca.edge_length,
                                length2=0.5 * mrca.edge_length)
            return True
            
        return False

    def run(self, input_tree, midpoint, outgroup, output_tree):
        """Reroot tree.
        
        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        midpoint : boolean
          Flag indicating tree should be rerooted at midpoint.
        outgroup : list
          Taxa forming outgroup.
        output_tree : str
          Name of file for rerooted tree.
        """
        
        self.logger.info('Reading input tree.')
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        
        if midpoint:
            self.logger.info('Rerooting tree at midpoint.')
            rooted = True
            tree.reroot_at_midpoint()
        else:
            self.logger.info('Rerooting tree with outgroup.')
            rooted = self.root_with_outgroup(tree, outgroup)
            
        if rooted:
            self.logger.info('Writing output tree.')  
            tree.write_to_path(output_tree, 
                                schema='newick', 
                                suppress_rooting=True, 
                                unquoted_underscores=True)
                        