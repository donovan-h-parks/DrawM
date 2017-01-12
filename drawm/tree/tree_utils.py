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


from drawm.tree.newick_utils import parse_label


def find_node(tree, label):
    """Find node in tree."""
    
    if '|' in label:
        try:
            node = tree.mrca(taxon_labels=label.split('|'))
        except:
            node = None
    else:
        node = tree.find_node_with_taxon_label(label)
        if not node:
            for n in tree.preorder_internal_node_iter():
                support, taxon, auxiliary_info = parse_label(n.label)
                if label == taxon:
                    node = n
                    break
                    
    return node
    

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
      

"""THESE NEED TO BE REMOVED OR SPECIFIED AS A VISUAL PROPERTY."""      
def _render_mean_tip_bl(self, tree, test_dist):
    """Render mean branch length to extent taxa (root to tip)."""

    contour_pts = []
    stack = [tree.seed_node]
    while stack:
        node = stack.pop()
                         
        # check if node meets mean branch length criterion
        if node.mean_tip_bl > test_dist:
            # node is above threshold so add children
            for c in node.child_node_iter():
                stack.append(c)
            continue
        else:
            # first node in lineage below threhold so
            # draw contour along this branch
            bl = node.mean_tip_bl
            pbl = node.parent_node.mean_tip_bl
            
            index = 0
            if bl != pbl:
                index = float(pbl - test_dist) / (pbl - bl)

            x = node.corner_x + index * (node.x - node.corner_x)
            y = node.corner_y + index * (node.y - node.corner_y)
                
            contour_pts.append((x, y))

    # draw contour
    path = self.dwg.path("M%d,%d" % contour_pts[0])
    path.fill(color='none', opacity=0.5, rule='evenodd')
    path.stroke(color='red', width=1)
    
    c = self.dwg.circle(center=contour_pts[0], r=self.node_radius)
    c.stroke(color='black')
    c.fill(color='red')
    c = self.dwg.add(c)
    
    # draw contour
    for pt in contour_pts[1:]:
        path.push("L%d,%d" % pt)
        
        c = self.dwg.circle(center=pt, r=self.node_radius)
        c.stroke(color='black')
        c.fill(color='red')
        c = self.dwg.add(c)
        
    self.dwg.add(path)
            
def _render_mean_tip_bl_test(self, tree, test_dist):
    """Render mean branch length to extent taxa."""

    contour_pts = []
    for node in tree.preorder_node_iter():
        if not node.parent_node:
            continue
            
        bl = node.mean_tip_bl
        pbl = node.parent_node.mean_tip_bl
        
        index = 0
        if bl <= test_dist <= pbl:
            if bl != pbl:
                index = float(pbl - test_dist) / (pbl - bl)

            x = node.corner_x + index * (node.x - node.corner_x)
            y = node.corner_y + index * (node.y - node.corner_y)
                
            contour_pts.append((x, y))

    # draw contour
    path = self.dwg.path("M%d,%d" % contour_pts[0])
    path.fill(color='none', opacity=0.5, rule='evenodd')
    path.stroke(color='green', width=1)
    
    c = self.dwg.circle(center=contour_pts[0], r=self.node_radius)
    c.stroke(color='black')
    c.fill(color='green')
    c = self.dwg.add(c)
    
    # draw contour
    for pt in contour_pts[1:]:
        path.push("L%d,%d" % pt)
        
        c = self.dwg.circle(center=pt, r=self.node_radius)
        c.stroke(color='black')
        c.fill(color='green')
        c = self.dwg.add(c)
        
    self.dwg.add(path)