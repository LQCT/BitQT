Technical Reference
===================


Basic help
----------
BitQT help is displayed in the console when typing **bitqt -h** ::

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

  As simple as that ;)


Arguments in Details
--------------------

``-traj (str):`` This is the only argument that is **always** required. Valid
extensions for trajectories are ``.dcd``, ``.dtr``, ``.hdf5``, ``.xyz``, ``.binpos``,
``.netcdf``, ``.prmtop``, ``.lh5``, ``.pdb``, ``.trr``, ``.xtc``, ``.xml``,
``.arc``, ``.lammpstrj`` and ``.hoomdxml``.

``-top (str):`` If trajectory format (automatically inferred from file extension)
includes topological information this argument is not required. Otherwise, user
must pass a path to a topology file. Valid topology extensions are  ``.pdb``,
``.pdb.gz``, ``.h5``, ``.lh5``, ``.prmtop``, ``.parm7``, ``.prm7``, ``.psf``,
``.mol2``, ``.hoomdxml``, ``.gro``, ``.arc`` and ``.hdf5``.

``-first (int, default=0):`` First frame to analyze (starting count from 0)

``-last (int, default=last):`` Last frame to analyze (starting count from 0).
Last frame is internally detected.

``-stride (int, default=1):`` Stride of frames to analyze. You might want to use
this argument to reduce the trajectory size when performing exploratory analysis.

``-sel (str, default='all'):`` Atom selection. BitQT inherits ``MDtraj``
syntax selection which is very flexible. Common cases are listed at the
`Syntax Selection`_ section. 

  
``-cutoff (int, default=2):`` RMSD cutoff for similarity measures given in Angstroms
(1 A = 0.1 nm).

``-min_clust_size (int, default=2):`` Minimum number of frames inside returned clusters.
0 is not a meaningful value and 1 implies an unclustered frame (no other frame is
similar to it). Greater values of this parameter may speed up the algorithm with
loss of uniformity in retrieved clusters.

``-nclust (int, default=all):`` Maximum number of calculated clusters. Change the default
for a better performance whenever you only need to inspect the first clusters.

``-ref (int, default=0):`` Reference frame to align trajectory.

``-odir (str, default="./bitQT_outputs"):`` Output directory to store analysis.
BitQT checks for outdir existence to avoid overwriting it.


Syntax Selection 
----------------

BitQT inherits atom selection syntax from **MDTraj** which is similar to that
in VMD. We reproduce below some of the **MDTraj** examples. Note that in BitQT
all keywords (or their synonyms) string are passed directly to ``-sel`` argument.
For more details on possible syntax, please refer to
`MDTraj original documentation <http://mdtraj.org/1.9.4/atom_selection.html>`_.

BitQT recognizes the following keywords.

=============    ========================   =========      ================================================================
Keyword          Synonyms                   Type           Description
-------------    ------------------------   ---------      ----------------------------------------------------------------
``all``          ``everything``             ``bool``       Matches everything
``none``         ``nothing``                ``bool``       Matches nothing
``backbone``     ``is_backbone``            ``bool``       Whether atom is in the backbone of a protein residue
``sidechain``    ``is_sidechain``           ``bool``       Whether atom is in the sidechain of a protein residue
``protein``      ``is_protein``             ``bool``       Whether atom is part of a protein residue
``water``        ``is_water``, ``waters``   ``bool``       Whether atom is part of a water residue
``name``                                    ``str``        Atom name
``index``                                   ``int``        Atom index (0-based)
``type``         ``element``, ``symbol``    ``str``        1 or 2-letter chemical symbols from the periodic table
``mass``                                    ``float``      Element atomic mass (daltons)
``residue``      ``resSeq``                 ``int``        Residue Sequence record (generally 1-based, but depends on topology)
``resid``        ``resi``                   ``int``        Residue index (0-based)
``resname``      ``resn``                   ``str``        Residue name
``rescode``      ``code``, ``resc```        ``str``        1-letter residue code
``chainid``                                 ``int``        Chain index (0-based)
=============    ========================   =========      ================================================================


Operators
+++++++++

Standard boolean operations (``and``, ``or``, and ``not``) as well as their
C-style aliases (``&&``, ``||``, ``!``) are supported. The expected logical
operators (``<``, ``<=``, ``==``, ``!=``, ``>=``, ``>``) are also available, as
along with their FORTRAN-style synonyms (``lt``, ``le``, ``eq``, ``ne``,
``ge``, ``gt``).

Range queries
+++++++++++++

Range queries are also supported. The range condition is an expression of
the form ``<expression> <low> to <high>``, which resolves to ``<low> <=
<expression> <= <high>``.  For example ::

    # The following queries are equivalent
    -sel "resid 10 to 30"
    -sel "(10 <= resid) and (resid <= 30)"

