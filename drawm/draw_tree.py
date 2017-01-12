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
import logging

import dendropy
import svgwrite

from drawm.svg.bootstrap_props import BootstrapProps
from drawm.svg.contour_props import ContourProps
from drawm.svg.label_props import LabelProps
from drawm.svg.lineage_props import LineageProps
from drawm.svg.tree_props import TreeProps


class DrawTree(object):
    """Create SVG image of tree in Newick format."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')
        
    def _read_config_file(self, config_file):
        """Read configuration information.
        
        Parameters
        ----------
        config_file : str
            File indicating location of all property files.
        """
        
        props = set(['bootstrap_props',
                        'collapse_props',
                        'contour_props',
                        'label_props',
                        'lineage_props',
                        'tree_props'])
                        
        config_file_dir = os.path.split(os.path.abspath(config_file))[0]
        
        prop_files = {}
        for line in open(config_file):
            if line[0] == '#' or not line.strip():
                continue
                
            prop, prop_file = [x.strip() for x in line.strip().split('=')]
            if prop in props:
                # get absolute path to property file
                prop_file = os.path.join(config_file_dir, prop_file)
                if not os.path.exists(prop_file):
                    self.logger.error('Could not find %s property file: %s' % (prop, prop_file))
                    sys.exit()
                prop_files[prop] = prop_file
            else:
               self.logger.warning('Unrecognized line in configuration file: %s' % line.strip())
               
        return prop_files

    def render(self, 
                input_tree, 
                config_file, 
                width, 
                height,
                dpi,
                output_prefix):
        """Render tree.
        
        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        config_file : str
          File specifying path to all property files.
        width : float
          Width of image.
        height : float
          Height of image.
        dpi : int
          Resolution of image (dots per inch).
        output_prefix : str
          Prefix for output files.
        """
        
        tree_name = os.path.split(os.path.basename(input_tree))[0]

        # setup SVG file
        canvas_width = width * dpi
        canvas_height = height * dpi
        font_size = int(8 * (float(dpi)/90) + 0.5)
        
        self.logger.info('Setting up SVG file.')
        svg_output = output_prefix + '.svg'
        dwg = svgwrite.Drawing(filename=svg_output, 
                                    size=(canvas_width, canvas_height))
        dwg.set_desc(title='DrawM rendering of %s' % tree_name, desc=tree_name)
        dwg.canvas_width = canvas_width
        dwg.canvas_height = canvas_height
        
        background = "rgb(" + str(255) + "," + str(255) + "," + str(255) + ")"
        dwg.add(dwg.rect(insert=(0, 0), 
                        size=('100%', '100%'), 
                        rx=None, ry=None, 
                        fill=background, 
                        stroke='none'))
                        
        # read configuration files
        self.logger.info('Reading configuration files.')
        prop_files = self._read_config_file(config_file)
        
        tree_props = TreeProps(prop_files['tree_props'], 
                                prop_files['collapse_props'],
                                dwg,
                                dpi)
                                            
        bootstrap_props = BootstrapProps(prop_files['bootstrap_props'],
                                            dwg,
                                            tree_props.height,
                                            dpi,
                                            font_size)
                                            
        contour_props = ContourProps(prop_files['contour_props'],
                                        dwg,
                                        tree_props.height,
                                        dpi,
                                        font_size)
                                            
        lineage_props = LineageProps(prop_files['lineage_props'],
                                        dwg,
                                        dpi,
                                        font_size)
                                            
        label_props = LabelProps(prop_files['label_props'],
                                        dwg,
                                        tree_props.height,
                                        dpi,
                                        font_size)
        
        # read and ladderize tree
        tree = tree_props.read_tree(input_tree)
        tree_props.layout(tree)
                                           
        contour_props.decorate(tree)
                                                           
        contour_props.render_contour(tree)
        contour_props.render_legend()
        
        lineage_props.render(tree)
   
        tree_props.render(tree)
        tree_props.render_scale_bar()
        tree_props.render_scale_lines()
        
        bootstrap_props.render(tree)
        bootstrap_props.render_legend()
        
        label_props.render(tree)
        
        self.logger.info('Saving SVG image.')
        dwg.save()
        
        self.logger.info('Saving PNG image.')
        png_output = output_prefix + '.png'
        os.system('inkscape -e %s -d %d -z -w %d -h %d %s' % (png_output,
                                                                dpi,
                                                                dwg.canvas_width, 
                                                                dwg.canvas_height, 
                                                                svg_output))