import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mission_sim",
    version="0.0.1",
    author="Andrew Inglis & Steven Christe",
    author_email="steven.christe@nasa.gov",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
)
