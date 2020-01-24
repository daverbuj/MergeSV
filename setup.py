import setuptools

with open("README.md", "r") as fh:
 long_description = fh.read()

setuptools.setup(
    name="mergesv",
    version="1.4",
    author="Dan Averbuj",
    author_email="dan.averbuj@gmail.com",
    description="Merge overlapping SVs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/daverbuj/MergeSV",
    requires=['csv', 'argparse', 'operator'],
    package_dir={'mergesv': 'mergesv/'},
    packages=setuptools.find_packages(),
    scripts=['mergesv/mergesv'],
    python_requires='>=3.6',
)