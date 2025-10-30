.. _sec-changelog:

Changelog
=========

Version 1.1.9
-------------

* Fix for timestamp parsing of Digital Micrograph files (thanks to Hüseyin Celik)
* Fix for deprecated numpy 2.0 functions  (thanks to Hüseyin Celik)

Version 1.1.8
-------------

* Dropped Python 2.X support
* Fix deprecated functions
* Added CI/CD support for GitLab

Version 1.1.7
-------------

* Output file attribute 'reconstructionCutOff2(nm2)' was set erroneously to isotropic cut off. Now squared anisotropic
  cut off is output as stated in the documentation.

Version 1.1.6
-------------

* Parameters 'object_names' and 'empty_names' now also allow file name lists instead of the printf-based syntax

Version 1.1.5 (unreleased)
--------------------------

* New parameter 'align_cut_off' allows specification of a dedicated raw alignment cut_off

Version 1.1.4
-------------

* Fixed bug, which caused the 'path' parameter to be ignored for some input files
* Fixed bug in configurable file importer, which occured when the deprecated syntas was used
* New parameter 'cut_off2' introduced, which allows nonisotropic reconstruction masks

Version 1.1.3
-------------

* Fixed a bug, which might have caused slight phase gradients over the image, if a ROI was used
* Display hint for use of PyFFTW
* Acceleration voltage for defocus propagation now taken from parameter file instead of image files, when present

Version 1.1.2
-------------

* Single hologram reconstructions do not require printf-format codes in input file paramters
* Include parameters in output file
* Include holoaverage version in output file
* Display global holograms amplitudes instead of intensities during reconstruction
* Made syntax of configurable file importer more URL-like.

Version 1.1.1
-------------

* Deprecation warnings are displayed by default
* Parameter 'output' was renamed to 'output_name'
* New parameter 'output_prefix' allows control over dataset names in output
* Added program overview to documentation
* Configurable image file importer
* Raw input image files
* New parameter 'enable_raw_align' to allow control over raw alignment
* Commandline option to print version number
* Improved handling of DM3 tag names encoding
