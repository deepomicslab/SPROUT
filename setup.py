from setuptools import setup
from setuptools import find_packages

setup( name = 'SPROUT',
version = '2.0.0',
description='SPROUT: spectral sparsification helps restore the spatial structure at single-cell resolution.',
url='https://github.com/deepomicslab/SPROUT/',
author='Jingwan WANG',
author_email='wanwang6-c@my.cityu.edu.hk',
license='MIT',
packages=find_packages(),
install_requires=['umap-learn',
'loess','tasklogger','time',
'random2','scipy']
)
