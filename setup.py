import setuptools
from setuptools.command.install import install
from io import open
import os

install_requires = ["numpy", "matplotlib", "biom-format", "rpy2"]
setuptools.setup(
    name='lefse',
    version='1.1.2',
    author='Nicola Segata',
    author_email='nicola.segata@unitn.it',
    url='http://github.com/SegataLab/lefse/',
    scripts = ['lefse_plot_cladogram.py','lefse_format_input.py','lefse_plot_features.py','lefse_plot_res.py', 'lefse.py','lefse_run.py'],
    packages = ['lefsebiom'],
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    description='Hclust2 is a handy tool for plotting heat-maps with several useful options to produce high quality figures that can be used in publication',
    install_requires=install_requires
)
