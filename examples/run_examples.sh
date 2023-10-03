#!/bin/sh

drawm draw example.tree ./simple_rectangular/simple_rectangular.cfg simple_rectangular
drawm draw example.tree ./simple_circular/simple_circular.cfg simple_circular
drawm draw example.tree ./lineage_rectangular/lineage_rectangular.cfg lineage_rectangular
drawm draw example.tree ./lineage_circular/lineage_circular.cfg lineage_circular
drawm draw example.tree ./collapse_rectangular/collapse_rectangular.cfg collapse_rectangular
drawm draw example.tree ./collapse_circular/collapse_circular.cfg collapse_circular
drawm draw example.tree ./contour_rectangular/contour_rectangular.cfg contour_rectangular
drawm draw example.tree ./contour_circular/contour_circular.cfg contour_circular
drawm draw archaea.tree ./archaea_collapse/archaea.cfg archaea_collapse
drawm draw archaea.tree ./archaea_contour/archaea.cfg archaea_contour
drawm draw archaea.tree ./archaea_lineages/archaea.cfg archaea_lineages --height 100
drawm draw archaea.tree ./archaea_example/archaea.cfg archaea_example
