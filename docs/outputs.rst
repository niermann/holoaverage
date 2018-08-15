.. _sec-outputs:

Outputs
=======

Within the parameter file the output filename is specified by the :ref:`param-output_name` parameter. This output file
will always be a HDF5 file. The dataset names described below are prefixed by the value of the :ref:`param-output_prefix`
parameter (which is empty by default). As example, if the value of ``output_prefix`` would be ``alpha_``, the dataset
``data`` is saved as ``alpha_data`` in the output file. By using the prefix multiple outputs can be written to the same
HDF5 file.

The output file typically holds three datasets:

The dataset ``data`` contains the averaged reconstruction of the normalized and drift aligned object holograms. If the defocus value
for the individual holograms (as specified by the parameters :ref:`param-defocus_first` and :ref:`param-defocus_step`) is
not zero, the holograms are propagated to the zero defocus. The dataset ``empty`` contains the averaged
reconstruction of the empty holograms.

As HDF5 files do not support datasets of complex numbers, the datasets are stored as a compound data type of 8 bytes
length. The compound has two members, both 32 bit floats. The first member ``r`` at byte offset 0, contains the real
part of the complex numbers, the second member ``i`` at byte offset 4 the imaginary part. This convention is recognized
by HDF5 python library (`<http://www.h5py.org>`_) and the HDF5 Digital Micrograph plugin (
`<https://github.com/niermann/gms_plugin_hdf5>`_).

The dataset ``variance`` contains the estimate of the per-pixel variance of the object hologram series (drift aligned
and propagated to zero defocus), this is stored as 32 bit float.

The datasets ``data`` and ``variance`` are not present, when the object hologram reconstruction is disabled by omitting
the :ref:`param-object_names` parameter. The dataset ``empty`` is not present, when the empty hologram reconstruction
is disabled by omitting the :ref:`param-empty_names` parameter.

If the parameter :ref:`param-output_series` is set, additionally a group ``series`` is present in the output file,
which contains datasets ``000``, ``001``, ``002``, ... . These contain the reconstructions of the individual holograms
in the object series. The dataset ``000`` refers to the first hologram (as specified by the parameter
:ref:`param-object_names` and :ref:`param-object_first`). The consecutive numbers refer to the consecutive holograms
in the series. These are also stored as complex valued datasets.

If the parameter :ref:`param-output_aligned` is set, additional a group ``aligned_rois`` is present in the output file,
which contains datasets ``000``, ``001``, ``002``, etc. These contain the region of interest (parameter
:ref:`param-roi`) as tracked across the object series in the raw-alignment step. The datatype of the dataset is
typically the datatype used in the image files.

Additional the datasets/groups have attributes, which contain further metadata. Not all the
following attributes are present in all datasets. This list is incomplete.

=========================== =================== ================================================================================
Name                        Type                Description
=========================== =================== ================================================================================
holoaverage_version         String              Version number of holoaverage used for averaging.
holoaverage_param           String              Parameter string passed to holoaverage (JSON).
align_roi                   List of 4 Ints      Region of object holograms used for raw alignment (see :ref:`param-align_roi`)
binning                     List of 2 Ints      (X, Y) Detector binning of series
carrier(nm-1)               List of 2 Floats    Spatial frequency of the reconstructed side band in 1/nm
convergence                 List of Floats      Total squared residual for each iteration of averaging procedure
defocus(nm)                 List of N Floats    Defocus of the individual holograms (after alignment) in nm
detector                    String              Name of detector (according to image files)
dim_offset                  List of 2 Floats    (X, Y) offset (according to image files)
dim_scale                   List of 2 Floats    (X, Y) sampling
dim_unit                    List of 2 Strings   (X, Y) units for ``dim_scale`` and ``dim_offset``
error                       List of N Floats    Squared residual between individual reconstruction and average
factor                      List of Nx2 Floats  (real, imaginary) global amplitudes for individual holograms
microscope                  String              Name of microscope (according to image files)
raw_shift                   List of Nx2 Ints    (X, Y) shift of individual holograms in pixels after raw-alignment
reconstructionCutOff2(nm2)  Float               Squared cut-off frequency (in 1/nm^2)
reconstructionMaskType      ...                 Mask type as used for cutoff (see :ref:`param-filter_func`)
roi                         List of 4 Ints      Reconstructed region of object holograms (see :ref:`param-roi`)
shift(nm)                   List of Nx2 Floats  (X, Y) shift of individual holograms in nm (after fine-alignment)
tilt(1/nm)                  List of Nx2 Floats  (X, Y) tilt of individual holograms in 1/nm (after alignment)
voltage(kV)                 Float               Acceleration voltage in kV
=========================== =================== ================================================================================
