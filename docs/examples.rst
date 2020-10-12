Example parameter files
=======================

Only reconstruct empty hologram
-------------------------------

This will reconstruct only the empty hologram. This can be used in for instance
in subsequent reconstructions (see other examples)::

    {
        // Leave "object_names" unset

        // Usual empty hologram input
        "empty_names" : "empty.%d.dm3",
        "empty_first" : 1,
        "empty_last" : 10,

        // Size (in px) used for reconstruction of "empty" holograms. Required.
        "empty_size" : 1536,

        // X, Y Position of side band in FFT pixels (origin is in center). Required.
        "sideband_pos" : [1011, 1091],

        // Output file name (will be HDF5 file). Empty hologram will be in dataset "empty"
        "output_name" : "my_empty_reco.hdf5",

        // Mask type (see FilterFunction for details). Defaults to "EDGE"
        "filter_func" : ["BUTTERWORTH", 14],

        // Reconstruction cutoff in 1/nm
        "cut_off" : 3.0,

        // Parameterization for MTF
        // "mtf" : ...
    }


Use pre-reconstructed empty hologram for normalization
------------------------------------------------------

This example uses a prereconstructed empty hologram for normalization
(also see :ref:`sec-normalization`)::

    {
        // Object holograms (as usual). Here series of DM3 files
        "object_names" : "hologram.%d.dm3",     // Filename format
        "object_first" : 1,                     // Index of first
        "object_last" : 5,                      // Index of last (inclusive)

        // Use pre-reconstructed empty hologram for normalization, here dataset from HDF5 file.
        "empty_override" : "somefile.hdf5?dataset=empty",

        // Reconstruction parameters
        "object_size" : 1536,                   // Reconstruction size in px
        "sideband_pos" : [1011, 1091],          // X, Y Position of side band in FFT pixels (origin is in center).
        "cut_off" : 3.0,                        // Reconstruction cut off in 1/nm
        "filter_func" : ["BUTTERWORTH", 14],    // Mask type

        // Optional reconstruction region (L, T, R, B). Defaults to full region.
        //"roi" : [166, 388, 1701, 1923],

        // Output file name (will be HDF5 file). Required.
        "output_name" : "my_reco.hdf5",

        // Parameterization for MTF
        //"mtf" : ...
    }


Use predetermined camera distortions for normalization
------------------------------------------------------

This example creates an synthetic empty hologram for normalization. The
synthetic empty hologram is creacted from predetermined camera distortions
(also see :ref:`sec-normalization`)::

    {
        // Object holograms (as usual). Here series of DM3 files
        "object_names" : "hologram.%d.dm3",     // Filename format
        "object_first" : 1,                     // Index of first
        "object_last" : 5,                      // Index of last (inclusive)

        // Enable synthetic empty holograms
        "synthesize_empty": true,

        // Two datasets, same size as holograms, with displacements in px for each pixel
        // First in X direction, Second in Y direction
        "camera_distortions": ["camera.hdf5?dataset=dx", "camera.hdf5?dataset=dy"],

        // Reconstruction parameters (as usual)
        "object_size" : 1536,                   // Reconstruction size in px
        "sideband_pos" : [1011, 1091],          // X, Y Position of side band in FFT pixels (origin is in center).
        "cut_off" : 3.0,                        // Reconstruction cut off in 1/nm
        "filter_func" : ["BUTTERWORTH", 14],    // Mask type

        // Optional reconstruction region (L, T, R, B). Defaults to full region.
        //"roi" : [166, 388, 1701, 1923],

        // Output file name (will be HDF5 file). Required.
        "output_name" : "my_reco.hdf5",

        // Parameterization for MTF
        //"mtf" : ...
    }
