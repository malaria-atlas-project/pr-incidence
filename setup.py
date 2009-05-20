# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from setuptools import setup
from numpy.distutils.misc_util import Configuration
import os
config = Configuration('pr_incidence',parent_package=None,top_path=None)

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(  version="0.1",
            description="The Malaria Atlas Project's PR-incidence code.",
            author="Anand Patil", 
            author_email="map@map.ox.ac.uk",
            url="www.map.ox.ac.uk",
            packages=['pr_incidence'],
            license="Creative Commons BY-NC-SA",
            **(config.todict()))
