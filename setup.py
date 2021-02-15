import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="bitqt",
    version="0.0.1",
    description="A Graph-Based Approach to the  Quality Threshold Clustering of Molecular Dynamics",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/LQCT/BitQT.git,
    author="Roy Gonzalez-Aleman",
    author_email="roy_gonzalez@fq.uh.cu",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    packages=["bitqt"],
    include_package_data=True,
    install_requires=["feedparser", "html2text"],
    entry_points={
        "console_scripts": [
            "realpython=reader.__main__:main",
        ]
    },
)
