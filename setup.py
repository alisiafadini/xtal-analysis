#coding: utf8

"""

Setup for xtal_analysis

"""

from glob import glob

try:
	from setuptools import setup

except ImportError:
	from disutils.core import setup

    
setup(name='xtal_analysis',
      author='Alisia Fadini',
      packages=['xtal_analysis'],
      package_dir={'xtal_analysis': 'xtal_analysis'})
