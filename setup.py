# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys

from taxcheck import __version__ as _version
here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

	name="taxcheck"
	version = _version

	if sys.version_info.major != 3:
		raise EnvironmentError("""{toolname} is a python module that requires python3, and is not compatible with python2.""".format(toolname=name))

	setup(
		name=name,
		version=version,
		description=description,
		long_description=long_description,
		# url="https://github.com/cschu/gff_quantifier",
		author="Christian Schudoma",
		author_email="christian.schudoma@embl.de",
		license="MIT",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Scientific Engineering :: Bio/Informatics",
			"License :: OSI Approved :: MIT License",
			"Operating System :: POSIX :: Linux",
			"Programming Language :: Python :: 3.7"
		],
		zip_safe=False,
		keywords="pathseq postprocess taxonomy check",
		packages=find_packages(exclude=["test"]),
		install_requires=[
			"ete3",
			"six",
			"numpy",
			"biopython",
		],
		entry_points={
			"console_scripts": [
				"taxcheck=taxcheck.__main__:main",
			],
		},
		package_data={},
		include_package_data=True,
		data_files=[],
	)
