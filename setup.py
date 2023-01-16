from setuptools import setup
from setuptools import find_packages

setup( name = 'STORM',
version = '1.0.0',
description='STORM: A software package restores the spatial structure from single-cell transcriptome with spectral graph sparsification.',
url='https://github.com/deepomicslab/STORM',
author='Jingwan WANG',
author_email='wanwang6-c@my.cityu.edu.hk',
license='MIT',
packages=find_packages(),
install_requires=['umap-learn',
'loess','tasklogger','time',
'random2','scipy']
)
