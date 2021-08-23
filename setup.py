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
    packages = ['lefse', 'lefsebiom'],
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    entry_points={
        'console_scripts': [
            'lefse_format_input.py = lefse.lefse_format_input:format_input',
            'lefse_plot_cladogram.py = lefse.lefse_plot_cladogram:plot_cladogram',
            'lefse_plot_features.py = lefse.lefse_plot_features:plot_features',
            'lefse_plot_res.py = lefse.lefse_plot_res:plot_res',
            'lefse_run.py  = lefse.lefse_run:lefse_run',
            'lefse2circlader.py = lefse.lefse2circlader:lefse2circlader',
            'qiime2lefse.py = lefse.qiime2lefse:qiime2lefse'
        ]
    },
    description="""
    LEfSe determines the features (organisms, clades, operational taxonomic units, genes, or functions)
    most likely to explain differences between classes by coupling standard tests for statistical significance with additional tests encoding
    biological consistency and effect relevance.""",
    install_requires=install_requires
)
