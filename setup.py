from setuptools import setup, find_packages

setup(
    name="pronet",
    version="0.1.1",
    author="Sarah",
    author_email="alsayedsarah01@gmail.com",
    description="A Python package for visualizing protein interaction networks from STRING database.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/pronet/",
    packages=find_packages(),
    install_requires=[
        "requests",
        "networkx",
        "matplotlib"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
