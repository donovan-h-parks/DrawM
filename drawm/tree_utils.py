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


def dist_to_ancestor(child, ancestor):
    """Calculate distance from child node to ancestor node."""
    
    d = 0.0
    cur_node = child
    while cur_node != ancestor:
        d += cur_node.edge.length
        cur_node = cur_node.parent_node
    
    return d

    
def decorate_mean_branch_length(tree):
    """Decorate tree with mean branch length to tips (mean_tip_bl)."""
    
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            node.mean_tip_bl = 0.0
            continue
            
        # check if node meets mean branch length criterion
        dists_to_tips = []
        for t in node.leaf_iter():
            dists_to_tips.append(dist_to_ancestor(t, node))
            
        node.mean_tip_bl = np_mean(dists_to_tips)