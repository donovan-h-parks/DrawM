#!/usr/bin/env python

from setuptools import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'drawm', 'VERSION'))
    return versionFile.read().strip()

if __name__ == '__main__':

    dirName = os.path.dirname(__file__)
    if dirName and os.getcwd() != dirName:
        os.chdir(dirName)

    setup(
        name='drawm',
        version=version(),
        author='Donovan Parks',
        author_email='donovan.parks@gmail.com',
        packages=['drawm', 'drawm.svg', 'drawm.tree'],
        scripts=['bin/drawm'],
        package_data={'drawm' : ['VERSION', './distributions/*.txt']},
        url='http://pypi.python.org/pypi/drawm/',
        license='GPL3',
        description='A toolkit for creating publication-quality images of phylogenetic trees.',
        install_requires=[
            "dendropy>=4.0.0",
            "biolib>=0.0.32",
            "svgwrite>1.1.19"],
    )
