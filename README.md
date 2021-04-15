# BitQT
> A Graph-Based Approach to the  Quality Threshold Clustering of Molecular Dynamics

BitQT is a Python command-line interface (CLI) conceived to speed up
the Heyer's Quality Threshold (QT) clustering [1] of long Molecular Dynamics.
The package implements a heuristic approach to [this exact variant
of QT](https://doi.org/10.1021/acs.jcim.9b00558).

## BitQT Home Page


BitQT´s latest documentation including usage examples, tutorials, benchmarks, etc. is available [here](https://bitqt.readthedocs.io).  


## Installation

There are some easy-to-install dependencies you must have before running BitQT. MDTraj (mandatory) will perform the heavy RMSD calculations,
while VMD (optional) will help with visualization tasks. The rest of the dependencies (listed below) will be automatically
managed by BitQT.


#### 1. **MDTraj**

It is recommended that you install __MDTraj__ using conda.

`conda install -c conda-forge mdtraj`

#### 2. **BitQT**

+ __Via **pip**__


After successfully installing __MDTraj__ you can easily install BitQT and the rest of its dependencies using pip.

```pip install bitqt```


+ __Via **GitHub**__

```
git clone https://github.com/LQCT/bitqt
cd bitqt
python setup.py install
```
Then, you should be able to see BitQT help by typing in a console:

`bitqt -h`


#### 3. **VMD** and **VMD clustering plugin** (optional)

BitQT clusters can be visualized by loading a **.log**  file in VMD via a clustering plugin.
Please see the [VMD visualization tutorial](https://bitqt.readthedocs.io/en/latest/tutorial.html#visualizing-clusters-in-vmd) in the BitQT documentation web page.

Official site for VMD download and installation can be `found here <https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_.

Instructions on how to install the clustering plugin of VMD are `available here <https://github.com/luisico/clustering>`_.


## Basic Usage
You can display the BitQT sel-explanatory  basic usage by typing  ` bitclust -h` in the command line.

```
$ bitclust -h

usage: bitqt -traj trajectory [options]

BitQT: A Graph-based Approach to the Quality Threshold Clustering of Molecular
Dynamics

optional arguments:
  -h, --help           show this help message and exit

Trajectory options:
  -traj trajectory     Path to trajectory file [required]
  -top topology        Path to the topology file
  -first first_frame   First frame to analyze (counting from 0) [default: 0]
  -last last_frame     Last frame to analyze (counting from 0) [default: last
                       frame]
  -stride stride       Stride of frames to analyze [default: 1]
  -sel selection       Atom selection (MDTraj syntax) [default: all]

Clustering options:
  -cutoff k            RMSD cutoff [default: 2]
  -min_clust_size m    Minimum size of returned clusters [default: 2]
  -nclust n            Number of clusters to retrieve [default: 2]

Output options:
  -odir bitQT_outputs  Output directory to store analysis [default:
                       bitQT_outputs]
```

## Citation (work in-press)

If you make use of BitQT in your scientific work, [cite it ;)]()

## Release History

* 0.0.1
    * First Release (academic publication)

## Licence

**BitQT** is licensed under GNU General Public License v3.0.

## Reference

[1] Heyer, L. J.; Kruglyak, S.; Yooseph, S. Exploring Expression Data Identification and Analysis of Coexpressed Genes. Genome Res. 1999, 9 (11), 1106–1115.


