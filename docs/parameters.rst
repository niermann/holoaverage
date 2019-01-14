.. highlight:: javascript

.. _sec-parameters:

Parameters
==========

File format
-----------

Parameter files for *holoaverage* are plain text files with a JSON (`<http://www.json.org>`_) compatible syntax. Such files
can be edited by any text editor. The whole parameter set must be a JSON object. Beside pure JSON the parser also
allows Javascript-like comments (line and block comments), as well as trailing commas. See the :ref:`sec-tutorial` for
example files.

When a parameter file is provided, *holoaverage* expects it to be in UTF-8 encoding (see `Wikipedia article <https://en.wikipedia.org/wiki/UTF-8>`_
for details). For parameters from *stdin* it uses the standard Python encoding (you can use the `PYTHONIOENCODING
<https://docs.python.org/3/using/cmdline.html#envvar-PYTHONIOENCODING>`_ environment variable to override this behavior).

.. _sec-file_pathes:

File pathes / File formats
--------------------------

Some parameters are file pathes. All file pathes are relative to the directory, where the parameter file
is located. If the parameters are read from *stdin*, the pathes are relative to the current directory. This behavior
can be changed by setting the :ref:`param-path` parameter (see description of parameter for details).

Within a path forward slashes (``/``) as well as backslashes (``\\``) both are treated as path separators.
Internally all pathes are normalized to the platform's separator (e.g. backslashes on Windows platforms).

In JSON strings the backslash character is used as escape sequences with special meaning. For this reason, if you
wan't to write a single backslash (like in a Windows path separator) into a JSON string you always have to put in
a double backslash. For example the JSON string ``"this\\is\\my\\path"`` becomes the path ``this\is\my\path``.

Some file types allow to pass additional parameters to the file reader.
These parameters are separated from the file name by a question mark (``?``). The parameters itself are separated by an
ampersand (``&``). Parameters have a parameter name and a value, both separated by a equal sign (``=``).

.. versionchanged:: 1.1.2
    Separation of parameters by ampersand (instead of question marks), like in HTTP query strings.

For instance the path ``some_file.raw?xsize=1024&ysize=1024&dtype=int32`` is interpreted as file name ``some_file.raw``
with three parameters named ``xsize`` (value 1024), ``ysize`` (value 1024) and ``dtype`` (value ``int32``).



The file type is recognized from the extension. You can manually select the file type by passing the parameter ``type``
with an extension described below as value (e.g. ``file.hdf5?type=dm3`` is read as DM3-file even if its has the
extension of a different file type.

Supported file types are:

Digital Micrograph 3 Files
^^^^^^^^^^^^^^^^^^^^^^^^^^

:File Extension: ``dm3``
:Description: Sampling, Acceleration voltage and camera binning are read from file if possible.
:Parameters: None

Hierarchical Data Format 5
^^^^^^^^^^^^^^^^^^^^^^^^^^

:File Extension: ``hdf5``, ``h5``
:Description: Sampling, Acceleration voltage and camera binning must be given in parameters.
:Parameters: * *dataset* - The name of the dataset (required).

Raw binary files
^^^^^^^^^^^^^^^^

:File Extension: ``raw``
:Description: Sampling, Acceleration voltage and camera binning must be given in parameters.
:Parameters: * *xsize* - Width of the image in pixels(required).
    * *ysize* - Height of the image in pixels(required).
    * *dtype* - Datatype (required; numpy compatible description: ``int32``, ``uint16``, ``complex64``, ``float32``, ``F4``, ...)
    * *offset* - Number of bytes to skip at file head (optional, defaults to 0)
    * *swap_bytes* - Whether the bytes should be swapped (optional; 0 or 1; defaults to 0)

.. _sec-param_reference:

Parameter reference
-------------------

This reference gives a description of the parameters. The format field gives the expected type(s) of the parameter.
The actual types depend on whether the parameters are used in a JSON parameter file or as a Python object.
The corresponding types are:

    =========== ====================== ===============================
    Format      JSON parameter file    Python
    =========== ====================== ===============================
    None/Null   ``null``               ``None``
    Boolean     ``false`` or ``true``  ``False`` or ``True``
    Integer     Number type            Type ``int``
    String      String type            Type ``str``
    List        Array type             Type ``list``
    Dictionary  Object type            Type ``dict`` with ``str`` keys
    =========== ====================== ===============================

.. _param-adjust_defocus:

adjust_defocus
^^^^^^^^^^^^^^

:Parameter: ``adjust_defocus``
:Type: Optional (default is ``false``)
:Format: Boolean
:Description: Switch, which determines whether the object reconstructions should be aligned for defocus variations.
    Due to instabilities of the microscope's
    stage or lens currents the defocus between the individual exposures of the series might drift. When this switch
    is set to ``true``, the program tries to detect defocus deviations in the object hologram series.

.. _param-adjust_shift:

adjust_shift
^^^^^^^^^^^^

:Parameter: ``adjust_shift``
:Type: Optional (default is ``true``)
:Format: Boolean
:Description: Switch, which determines whether the object reconstructions should be aligned for specimen drift.
    When this switch is set to ``true``, the program tries to shift all object holograms to a common position during
    the averaging step. This "fine" alignment is performed independently from the "raw" alignment, which is controlled
    by the parameter :ref:`param-enable_raw_alignment`.

.. _param-adjust_tilt:

adjust_tilt
^^^^^^^^^^^

:Parameter: ``adjust_tilt``
:Type: Optional (default is ``false``)
:Format: Boolean
:Description: Switch, which determines whether the object reconstructions should be aligned for drift of the sideband
    position. Such a drift might occur when the voltage supply of the biprism is not stable. Usually this alignment is
    not needed.

.. _param-align_roi:

align_roi
^^^^^^^^^

:Parameter: ``align_roi``
:Type: Optional (by default region from parameter :ref:`param-roi` is taken)
:Format: List of four integers
:Unit: Pixels
:Description: ``[left, top, right, bottom]`` pixel coordinates of the region used for raw alignment of the object
    holograms. This region can be specified independently from the reconstruction region (as given by :ref:`param-roi`).

    If this parameter is not given the reconstruction region :ref:`param-roi` is also used for raw alignment.

    .. deprecated:: 1.1
        Setting this parameter to ``null`` disables the raw alignment. Set the parameter :ref:`param-enable_raw_alignment`
        to ``false`` instead.

.. _param-binning:

binning
^^^^^^^

:Parameter: ``binning``
:Type: Optional (taken from input files by default).
:Format: Integer
:Description: Binning used for recording of the holograms. This parameter affects, how the parameterization of the MTF
    (see :ref:`param-mtf`) is interpreted. If this parameter is not given, the binning is taken from the image files.
    If the image files provide no binning, it is assumed to be one.

.. _param-camera_distortions:

camera_distortions
^^^^^^^^^^^^^^^^^^

:Parameter: ``camera_distortions``
:Type: Optional
:Format: List of two Strings
:Description: Per pixel displacements due to camera distortions. The optics of the camera itself produce small
    displacements. This parameter contains two filenames. The first filename contains an array with the X-displacement
    of each pixel. The second filename contains the Y-displacements. The referenced arrays must have the same dimensions as the
    holograms. The displacements are given in units of pixels. These displacements are only used, if the parameter
    :ref:`param-synthesize_empty` is set.

.. _param-cut_off:

cut_off
^^^^^^^

:Parameter: ``cut_off``
:Type: Mandatory
:Format: Floating point number
:Unit: Reciprocal nanometer (1/nm)
:Description: This parameter defines in combination with the parameter :ref:`param-filter_func`, how the masking of the
    sideband in Fourier space is done. This is typically the radius of the mask used. The smaller this is chosen,
    the lower the resolution of the reconstructions will be. However, smaller values will spatially average the
    reconstructions more, thus decreasing the noise present in the holograms (at the cost of larger spatial correlations).
    The value specified by this parameter is also taken as cut-off frequency for the low pass used in the raw alignment
    step. For the raw alignment low pass, always a hard aperture (edge function) is taken.
    Please note, that if a wrong :ref:`param-sampling` is specified, the value of this parameter does not refer to the
    correct spatial frequency.
    Instead of this parameter the parameter :ref:`param-cut_off2` can be specified.

.. _param-cut_off2:

cut_off2
^^^^^^^^

:Parameter: ``cut_off2``
:Type: Alternative to (:ref:`param-cut_off`)
:Format: 2x2 matrix of floating point numbers (list of two lists of two floats)
:Unit: Reciprocal nanometer squared (1/nm2)
:Description:
    This parameter extents the functionality of the parameter :ref:`param-cut_off` for non-isotropic masking.
    For a general description of the overall parameter see :ref:`param-cut_off`. For a masking with radius `a` along
    the major axis with an angle of `alpha` to the x-axis and a radius of `b` along the minor axis, specify

    .. math::
        \begin{multline}
        R = \left[ \begin{array}{cc}
        \cos(\alpha) & \sin(\alpha) \\
        -\sin(\alpha) & \cos(\alpha) \\
        \end{array}\right] \\
        \mathrm{cut\_off2} = R^T \cdot \left[ \begin{array}{cc}
        a^2 & 0 \\
        0 & b^2 \\
        \end{array}\right] \cdot R
        \end{multline}

    If this parameter is specified, the parameter :ref:`param-cut_off` must not be present.

    Raw alignment still uses isotropic filtering with the geometric mean of both radii as radius.

    .. versionadded:: 1.1.4

.. _param-defocus_first:

defocus_first
^^^^^^^^^^^^^

:Parameter: ``defocus_first``
:Type: Optional (default is 0.0 nm)
:Format: Floating point number
:Unit: Nanometers
:Description: Defocus of first object hologram (hologram with index given by :ref:`param-object_first`).
    Negative focus values refer to underfocus. The reconstructed (averaged) object hologram is propagated to the
    Gaussian focus (i.e. defocus of zero) during reconstruction. No propagation of the reconstructed hologram is
    performed, when the defocus of an hologram is given as zero. The empty holograms are never propagated.
    Please note, that if the sampling of the holograms (see :ref:`param-sampling`) or the acceleration voltage (see
    :ref:`param-voltage`) are wrongly specified, the propagation will be performed wrongly. Also note, that if the
    defocus is specified wrongly, the holograms will be be propagated to a different focus than the Gaussian one.

.. _param-defocus_step:

defocus_step
^^^^^^^^^^^^^

:Parameter: ``defocus_step``
:Type: Optional (default is 0.0 nm)
:Format: Floating point number
:Unit: Nanometers
:Description: Step of defocus between consecutive object holograms in the series. This is intended for the
    case that the hologram series is also a focal series, where every hologram has a different defocus.
    Defaults to 0.0 nm (all object holograms were taken at same defocus).

.. _param-empty_exclude:

empty_exclude
^^^^^^^^^^^^^^

:Parameter: ``empty_exclude``
:Type: Optional (default is empty list)
:Format: List of integers
:Description: A list of empty hologram indices, which should **not** be used for averaging. See
    :ref:`param-object_exclude` for the rationale of this parameter. By default this list is empty and all empty
    holograms in the given range are used.

.. _param-empty_first:

empty_first
^^^^^^^^^^^^

:Parameter: ``empty_first``
:Type: Mandatory
:Format: Integer
:Description: Index of first hologram in the empty hologram series.

.. _param-empty_last:

empty_last
^^^^^^^^^^^

:Parameter: ``empty_last``
:Type: Mandatory
:Format: Integer
:Description: Index of last hologram (inclusive) in the empty hologram series.

.. _param-empty_names:

empty_names
^^^^^^^^^^^

:Parameter: ``empty_names``
:Type: Mandatory
:Format: String
:Description: File name of empty hologram series. See :ref:`param-object_names` for the description of the format of this
    parameter.

    If the parameter ``empty_names`` is not present in the parameter file, no empty hologram series will be
    reconstructed and averaged. In this case, the parameters :ref:`param-empty_first`, :ref:`param-empty_last` are not
    needed.

.. _param-empty_override:

empty_override
^^^^^^^^^^^^^^

:Parameter: ``empty_override``
:Type: Optional
:Format: String
:Description: File name of empty hologram used for normalization. If this parameter is present in the parameter
    files the empty hologram will be read from this file (see :ref:`sec-file_pathes` for format) and the parameters
    :ref:`param-empty_names`, :ref:`param-empty_first`, :ref:`param-empty_last`, and :ref:`param-empty_size` are
    ignored.

.. _param-empty_size:

empty_size
^^^^^^^^^^^

:Parameter: ``empty_size``
:Type: Optional (default is given by parameter :ref:`param-object_size`)
:Format: Integer
:Unit: Pixels
:Description: Size of the reconstructed empty hologram. See :ref:`param-object_size` for details concerning this
    parameter. For normalization of the reconstructed object holograms the reconstructed empty hologram is interpolated
    to the size of the object holograms (before its cropped to the :ref:`param-roi` region) by zero-padding.
    If parameter :ref:`param-empty_size` is missing, it is substituted by :ref:`param-object_size`.

.. _param-enable_raw_alignment:

enable_raw_alignment
^^^^^^^^^^^^^^^^^^^^^

:Parameter: ``enable_raw_alignment``
:Type: Optional (default is ``true``)
:Format: Boolean
:Description: Enables the raw alignment. If the raw alignment is disabled, the region of interest is taken from the
    same area in each hologram of the object hologram series. Otherwise, the region of interest is tracked across the
    series.

    .. versionadded:: 1.1

.. _param-filter_func:

filter_func
^^^^^^^^^^^

:Parameter: ``filter_func``
:Type: Optional (default is ``"EDGE"``)
:Format: see below
:Description: This parameter gives the function that will be used in combination with the parameter
    :ref:`param-cut_off` for masking the sideband in Fourier space. The format of this parameter is either
    a string describing the filter function, or a list with the function name as first element and further parameters
    in the remaining list.

    If ``filter_func`` is ``"EDGE"``, an edge function is used. This corresponds to a hard mask at the ``cut_off``
    spatial frequency. If the edge function is chosen, you might observe "ringing" artifacts in the reconstructions
    especially at the borders or at "hot pixels".

    If ``filter_func`` is ``"GAUSSIAN"``, a Gaussian function is used. The Gaussian is chosen such that a ``1/e``
    fall-off is reached at the ``cut_off`` spatial frequency.

    If ``filter_func`` is ``["BUTTERWORTH", order]``, a Butterworth function of the given order is used. This
    corresponds to a soft mask at the ``cut_off`` spatial frequency. The lower the order of the Butterworth function is,
    the softer this filter becomes.

    If this parameter is not given, the edge function is used.

.. _param-mtf:

mtf
^^^

:Parameter: ``mtf``
:Type: Optional
:Format: List
:Description: Parameterization of the camera MTF. The reconstruction are corrected for the effects of MTF (by
    dividing the Fourier transformed holograms by the MTF). See :ref:`sec-mtf` for details on the specification
    of this parameter. If this parameter is not given, no MTF correction is performed.

.. _param-object_exclude:

object_exclude
^^^^^^^^^^^^^^

:Parameter: ``object_exclude``
:Type: Optional (default is empty list)
:Format: List of integers
:Description: A list of object hologram indices, which should **not** be used for averaging. Usually all holograms
    with indices between :ref:`param-object_first` and :ref:`param-object_last` (inclusive) are used for averaging. Any indices
    occurring in this list are not used. For example with ``object_first`` of ``1``, ``object_last`` of ``5``, and
    ``object_exclude`` set to ``[3, 4]`` only object holograms with indices ``1``, ``2``, and ``5`` are used, since
    indices ``3`` and ``4`` were explicitly excluded. By default, this list is empty and all object holograms in the
    given range are used.

.. _param-object_first:

object_first
^^^^^^^^^^^^

:Parameter: ``object_first``
:Type: Mandatory
:Format: Integer
:Description: Index of first hologram in the object hologram series.

.. _param-object_last:

object_last
^^^^^^^^^^^

:Parameter: ``object_last``
:Type: Mandatory
:Format: Integer
:Description: Index of last hologram (inclusive) in the object hologram series.

.. _param-object_names:

object_names
^^^^^^^^^^^^

:Parameter: ``object_names``
:Type: Mandatory
:Format: String
:Description: File name of object hologram series. Typically a series hologram file names contain an increasing number.
    The number in this parameter is encoded with the *printf*-style format rules (`old-style formating in python
    <http://docs.python.org/3/library/stdtypes.html#old-string-formatting>`_). For instance simple numbers can be
    expressed as ``%d`` and become ``1``, ``2``, ``3``, etc. If you want to have zero padded three digit numbers use
    ``%03d``, which becomes ``001``, ``002``, ``003``, etc. Due to this formatting rules you have to write a double
    percent sign (i.e. ``%%``) if you want a single ``%`` in your filename.

    If the parameter ``object_names`` is not present in the parameter file, only the empty hologram series will be
    reconstructed and averaged. In this case, the parameters :ref:`param-object_first`, :ref:`param-object_last`,
    and :ref:`param-object_size` are not needed.

.. _param-object_size:

object_size
^^^^^^^^^^^

:Parameter: ``object_size``
:Type: Mandatory
:Format: Integer
:Unit: Pixels
:Description: Size of the reconstructed object hologram. Reconstructed holograms always have same size in width and
     height. This size in pixels is given by this parameter. The :ref:`param-roi` of the object holograms is scaled
     to this size during the reconstruction (by cropping in Fourier space). This parameter should be larger than the
     diameter of filter used during the reconstruction (see :ref:`param-cut_off` parameter). For performance
     reasons a number with low prime factors should be chosen, e.g. prefer ``384 = 3 * 2^7`` over ``383`` (prime).

.. _param-only_phase:

only_phase
^^^^^^^^^^

:Parameter: ``only_phase``
:Type: Optional  (default is ``false``)
:Format: Boolean
:Description: Switch, which determines how the object reconstructions are normalized. When this parameter is ``true``,
    the normalization is performed by dividing the individual reconstructed object holograms by the reconstructed
    (and averaged) empty hologram. This normalizes the object holograms in amplitude in phase. However, if the
    reconstructed empty hologram contains regions, where the amplitude is very small, the normalization will cause
    artifacts. Such cases typically occur when the interference region, does not cover the whole image.
    When this parameter is ``true``, only the phases of the reconstructed holograms are normalized.

.. _param-output_aligned:

output_aligned
^^^^^^^^^^^^^^

:Parameter: ``output_aligned``
:Type: Optional (default is ``false``)
:Format: Boolean
:Description: When set to ``true``, the region of interest of the individual object holograms (before
    reconstruction) are also stored in the output file.

.. _param-output_name:

output_name
^^^^^^^^^^^

:Parameter: ``output_name``
:Type: Mandatory
:Format: String
:Description: Name of the output file. The output(s) will be always stored in HDF5 format.

    .. versionchanged:: 1.1
        The parameter was renamed from ``output`` to ``output_name``.

.. _param-output_prefix:

output_prefix
^^^^^^^^^^^^^

:Parameter: ``output_prefix``
:Type: Optional (Defaults to empty string)
:Format: String
:Description: Prefix to dataset names in output file. By using the prefix multiple outputs can be written to the same
    HDF5 file. Especially forward slashes can be used in :ref:`param-output_prefix` to create the outputs in sub-groups.
    As example, if the value of ``output_prefix`` would be ``alpha_``, the dataset ``data`` is saved as ``alpha_data``
    in the output file.

    .. versionadded:: 1.1

.. _param-output_series:

output_series
^^^^^^^^^^^^^^

:Parameter: ``output_series``
:Type: Optional (default is ``false``)
:Format: Boolean
:Description: When set to ``true``, also the individual object hologram reconstructions are stored in the output file.
    The averaged hologram (and the variance estimation obtained during averaging) are always stored in the output file.
    The individual reconstructions of the empty hologram series are never stored.

.. _param-path:

path
^^^^

:Parameter: ``path``
:Type: Optional (default is none)
:Format: String
:Description: All (relative) file names are relative to this path. Absolute file names are still absolute.
    If this parameter itself is not an absolute path, the path is taken relative to the path
    of the parameter file (current directory, if the parameters are read from *stdin*). By default this path is left
    empty, which means all file names are relative to the parameter file path (or the current directory, when the
    parameters are read from *stdin*; see :ref:`sec-file_pathes`).

.. _param-roi:

roi
^^^

:Parameter: ``roi``
:Type: Optional (default is full image region)
:Format: List of four integers.
:Unit: Pixels
:Description: ``[left, top, right, bottom]`` pixel coordinates of the region of interest (ROI) in the first object
    hologram (as given by parameter :ref:`param-object_first`). The ROI is always a rectangular region. In the raw
    alignment step (:ref:`sec-overview`) of the hologram series the position of this ROI is aligned to the drift of the object,
    such that always the same object region is taken from each hologram.

    The *left* and *top* pixel positions given here refer to the top, left corner in this rectangular region
    (inclusive). The *right* and *bottom* positions refer to the bottom, right corner (exclusive), which means they
    refer the pixel coordinate adjacent to right (bottom) edge of the ROI.
    X coordinates are going from left to right, Y coordinates are going form top to bottom. For performance reasons,
    the size of the ROI, i.e. ``right - left`` and ``bottom - top``, should have only low prime-factors, e.g. prefer
    ``384 = 3 * 2^7`` over ``383`` (prime).

    If this value is not given, the whole object hologram region is taken as ROI.

.. _param-sampling:

sampling
^^^^^^^^

:Parameter: ``sampling``
:Type: Optional (taken from input files by default)
:Format: Floating point number
:Unit: Nanometer per pixel
:Description: Sampling of the object and empty holograms. The number given by this parameter corresponds to the size
    of a single pixel of the holograms. If this parameter is not given, the sampling from the image files is taken.
    Otherwise this parameter overrides the sampling given in the files.

    Please note that all holograms, independently of being part of the object or empty series must have the same
    sampling. Also only image files with samplings given in nanometer per pixel are supported. If the sampling recorded
    in the image files is wrong (or the file format does not provides this metadata), the ``sampling`` parameter must
    be set explicitly.

.. _param-sideband_pos:

sideband_pos
^^^^^^^^^^^^

:Parameter: ``sideband_pos``
:Type: Mandatory
:Format: List of two floating point numbers.
:Unit: Pixels
:Description: ``[X,Y]`` position of the sideband in the Fourier transformed image files. When the discrete Fourier
    transform of the holograms is calculated and the Fourier transform is shifted such that the Fourier space origin
    is in the center of the transformed images (like the numpy commands ``np.fft.fftshift(np.fft.fft2(image))`` would
    do), this parameter refers to the pixel position of the sideband to be reconstructed.

.. _param-synthesize_empty:

synthesize_empty
^^^^^^^^^^^^^^^^

:Parameter: ``synthesize_empty``
:Type: Optional (default is ``false``)
:Format: Boolean
:Description: When set to ``true``, the reconstructed object hologram series is normalized by a synthetic empty
    hologram instead of an experimental empty hologram. The synthesized empty hologram is calculated from the provided
    camera distortions. If ``synthesize_empty`` is set, the parameters :ref:`param-camera_distortions` and
    :ref:`param-empty_size` must be also given. If ``synthesize_empty`` is set, other emtpy holograms (provided either
    by :ref:`param-empty_names` or :ref:`param-empty_override`) are ignored.

.. _param-voltage:

voltage
^^^^^^^^

:Parameter: ``voltage``
:Type: Optional (taken from input files by default)
:Format: Floating point number
:Unit: Kilovolts
:Description: Acceleration voltage used during acquisition of the holograms. If this parameter is not given it is taken
    from the holograms files. This parameter must be given explicitly, if the acceleration voltage cannot be read
    from the hologram files.

.. _sec-mtf:

Modulation Transfer Function
----------------------------

The modulation transfer function (MTF) of the camera used for acquisition of the individual holograms is specified
in parameterized form.

In the following, it is assumed the MTF is a 2 dimensional function :math:`M(q_x, q_y)` of the
two dimensional spatial frequency :math:`(q_x, q_y)`. A spatial frequency of +/-0.5 gives the Nyquist frequency of the
detector. The MTF consists then of two parts, one due to the binning into pixels, and the other part due to the beam
broadening within the detector/scintillator.

.. math::
    M(q_x, q_y) = \mathrm{sinc}(q_x) \mathrm{sinc}(q_y) \sum_n f_n(q)

The effect of the binning is described by the two *sinc* functions, here defined as

.. math::
    \mathrm{sinc}(q) = \sin(\pi q) / (\pi q).

The beam broadening in the above parameterization is described by a sum over functions :math:`f_n(q)`, where

.. math::
    q = \sqrt{q_x^2 + q_y^2}.

These functions are specified in the parameter file as a list of terms, where each term describes one function
:math:`f_n(q)`. The terms itself are again lists, where the first element always is a string describing the kind of
function and the other elements are parameters to the function.

Possible terms are:

    * ``["CONSTANT", A]``

    .. math::
        f(q) = A

    * ``["GAUSSIAN, A, B]``

    .. math::
        f(q) = A \exp(-B q^2)

    * ``["LORENTZIAN", A, B]``

    .. math::
        f(q) = A / (B + q^2)

As example, if the MTF of the detector is given by:

    .. math::
        M(q_x, q_y) = \mathrm{sinc}(q_x) \mathrm{sinc}(q_y) \left[ 0.8 \exp(-0.03 q^2) + 0.2 \right]

the parameterization as specified by the :ref:`param-mtf` parameter is

::

    mtf = [["GAUSSIAN", 0.8, 0.03], ["CONSTANT", 0.2]]

