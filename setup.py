import setuptools
from MiBFSSI.pipeline import __version__, __author__, __email__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MiBFSSI",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="Pipeline to automatically process BFSSI MiSeq output",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bfssi-forest-dussault/MiBFSSI/",
    packages=setuptools.find_packages(),
    scripts=['MiBFSSI/pipeline.py'],
    install_requires=[]
)
