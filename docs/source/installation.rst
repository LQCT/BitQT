Installation
============
There are some easy-to-install dependencies you must have before running BitQT.
MDTraj (mandatory) will perform the heavy RMSD calculations while VMD (optional)
will be helpful for visualization tasks. The rest of dependencies will be automatically 
managed by BitQT.


**MDTraj**
----------

It is recommended that you install ``mdtraj`` using conda. ::

  $ conda install -c conda-forge mdtraj

You can install ``mdtraj`` with ``pip``, if you prefer. ::

  $ pip install mdtraj


**BitQT**
---------

Via **pip**:

After succesful installation of ``mdtraj`` you can easily proceed to
install BitQT and the rest of its key dependencies using pip. ::

  $ pip install bitqt


Via **GitHub**:


  $ git clone https://github.com/LQCT/bitqt
  $ cd bitqt
  $ python setup.py install

Then, you should be able to see BitQT help by typing in a console: ::

  $ bitqt -h



**VMD** and **VMD clustering plugin** (optional)
------------------------------------------------
BitQT clusters can be visualized by loading a *.log* file in VMD via a clustering plugin.
This is described in section :ref:`vmd_tutorial`.

Official site for VMD download and installation can be `found here <https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_.

Instructions on how to install the clustering plugin of VMD are `available here <https://github.com/luisico/clustering>`_.

