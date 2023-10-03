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
__copyright__ = 'Copyright 2017'
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
from collections import defaultdict

import svgwrite

from drawm.svg.svg_utils import donut, render_label


class SymbolProps:
    """Visual attributed for symbols to display in columns after extent taxa."""

    def __init__(self, config_file, dwg, inch):
        """Read attributes from config file."""
        
        self.logger = logging.getLogger('timestamp')

        self.dwg = dwg
        self.inch = inch

        self.show_symbols = False
        self.symbols = {}
        self.symbol_file = None
        
        if not config_file:
            return # use default values
        
        with open(config_file) as f:
            prop_type = f.readline().strip()
            
            if prop_type != 'SYMBOLS':
                self.logger.error("[SymbolProps] Unexpected property type '%s' in %s." % (prop_type, config_file))
                sys.exit()
            
            for line in f:
                if line[0] == '#' or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                attribute, values = fields[0], fields[1:]
                if attribute == 'show_symbols':
                    self.show_symbols = (values[0] == 'True')
                elif attribute == 'symbol':
                    label, column, shape, color, symbol_size = values
                    self.symbols[label] = (int(column), shape, color, float(symbol_size) * self.inch)
                elif attribute == 'symbol_file':
                    self.symbol_file = os.path.join(os.path.split(config_file)[0], 
                                                        values[0])
                else:
                    self.logger.warning('[SymbolProps] Unexpected attribute: %s' % attribute)
                            
    def _read(self, symbol_file):
        """Read file indicating symbols for each extent taxa."""
        
        extent_symbols = defaultdict(lambda: defaultdict(int))
        for line in open(symbol_file):
            line_split = line.strip().split('\t')
            for s in map(str.strip, line_split[1].split(',')):
                extent_symbols[line_split[0]][s] += 1
            
        return extent_symbols
        
    def _draw_column_lines(self, tree, symbols, symbol_offset, symbol_group):
        """Draw lines between symbol columns."""

        for symbol_label in symbols:
            column, _shape, _color, symbol_size = symbols[symbol_label]
            if column == 0:
                continue
                
            x = tree.width + tree.start_x + symbol_offset + 3*symbol_size*column - 1.5*symbol_size
            y_start = tree.start_y - symbol_size
            y_end = tree.height + tree.start_y + symbol_size
            
            p = self.dwg.line(start=(x, y_start), end=(x, y_end))
            p.fill(color='none')
            p.stroke(color='grey', opacity=0.5, width=1)
                    
            symbol_group.add(p)
    
    def render(self, tree):
        """Render symbols."""
        
        if not self.show_symbols:
            return
            
        self.logger.info('Rendering symbols.')
        
        # read symbols for each extant taxa
        extent_symbols = self._read(self.symbol_file)

        symbol_group = svgwrite.container.Group(id='symbols')
        self.dwg.add(symbol_group)
        
        symbol_offset = 20
        self._draw_column_lines(tree, self.symbols, symbol_offset, symbol_group)
        
        for leaf in tree.leaf_node_iter():
            extent_id = leaf.taxon.label
            if extent_id in extent_symbols:
                for symbol_label, count in extent_symbols[extent_id].iteritems():
                    column, shape, color, symbol_radius = self.symbols[symbol_label]

                    symbol_offset = 20 # TBD: this needs to fall after all labels???
                    x = tree.width + tree.start_x + symbol_offset + 3*symbol_radius*column
                    y = leaf.y
                    if shape == 'circle':
                        s = self.dwg.circle(center=(x, y), r=symbol_radius)
                    elif shape == 'square':
                        s = self.dwg.rect(insert=(x-symbol_radius, y-symbol_radius),
                                            size=(2*symbol_radius, 2*symbol_radius))
                    else:
                        self.logger.warning('Symbol shape %s is currently not supported.' % shape)
                
                    s.fill(color=color)
                    s.stroke(color='grey')
                    symbol_group.add(s)
                    
                    if count > 1:
                        render_label(self.dwg, 
                                        x, 
                                        y, 
                                        0, 
                                        str(count), 
                                        6, 
                                        'black',
                                        middle_x=True,
                                        middle_y=True,
                                        group=symbol_group,
                                        id_prefix='symbols_text')
