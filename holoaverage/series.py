# This file is part of holoaverage.
# Copyright (c) 2018 Tore Niermann
#
# holoaverage is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with holoaverage.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import numpy as np

from .grid import Grid
from .fft import empty_aligned

__all__ = ('DataSet', 'AbstractSeries', 'Series', 'LazyLoadingSeries')


class DataSet(object):
    """
    A dataset.
    """

    # general attributes.
    general_attributes = ["unit", "scale", "offset", "dim_unit", "dim_scale", "dim_offset",
                          "voltage(kV)", "thickness(nm)", "exposure(s)", "binning", "space",
                          "microscope", "detector", "timestamp"]

    def __init__(self, shape, dtype=float, factory=empty_aligned):
        """Create empty dataset of given shape and data type.
        Arguments:
            shape    Tuple of Numbers
            dtype    Datatype
            factory  A method with arguments (shape, dtype) that creates an empty array.
        """
        self._array = factory(shape, dtype)
        self._attrs = {}

    @property
    def array(self):
        """Returns numpy array of dataset. Might return None, if
        DataSet is not implemented as an numpy array."""
        return self._array

    @property
    def attrs(self):
        """Attributes of dataset."""
        return self._attrs

    @property
    def shape(self):
        """Returns shape of dataset."""
        return self._array.shape

    @property
    def dtype(self):
        """Returns datatype of dataset."""
        return self._array.dtype

    def copy(self):
        """Returns deep copy of dataset."""
        result = DataSet(self._array.shape, self._array.dtype)
        result._array[...] = self._array
        result._attrs.update(self._attrs)
        return result

    @staticmethod
    def fromArray(array):
        """Creates dataset from numpy array."""
        result = DataSet(array.shape, array.dtype)
        result._array[...] = array
        return result


class AbstractSeries(object):
    """A Series is an array of DataSets of equal shape and type.

    Implementers must implement __getitem__ and __setitem__.
    """

    def __init__(self, indexShape, shape, dtype=float):
        """Create empty series of given shape and data type.

        Arguments:
            indexShape : shape of array of datasets
            shape : Shape of datasets
            dtype : Type of datasets
        """
        indexShape = np.cast[int](np.atleast_1d(indexShape))
        shape = np.cast[int](np.atleast_1d(shape))
        self._shape = tuple(shape)
        self._size = np.prod(shape)
        self._indexShape = tuple(indexShape)
        self._indexSize = np.prod(indexShape)
        self._dtype = np.dtype(dtype)
        self._attrs = {}
        self._grid = None

    @property
    def dtype(self):
        """Each of the serie's DataSets will have this type."""
        return self._dtype

    @property
    def shape(self):
        """Each of the serie's DataSets will have this shape."""
        return self._shape

    @property
    def size(self):
        """Each of the serie's DataSets will have this size."""
        return self._size

    @property
    def attrs(self):
        """Attributes of dataset."""
        return self._attrs

    @property
    def indexShape(self):
        """Shape of the series."""
        return self._indexShape

    @property
    def indexSize(self):
        """Number of DataSets in series."""
        return self._indexSize

    def __len__(self):
        return self._indexShape[0]

    @property
    def grid(self):
        """Create for data. Created from first dataset if not present."""
        if self._grid is None:
            self._grid = Grid.fromDataSet(self)
        return self._grid

    @grid.setter
    def grid(self, value):
        self._grid = value

    def __getitem__(self, index):
        raise NotImplementedError

    def __setitem__(self, index, value):
        raise NotImplementedError

    def _prepareDataSet(self, value):
        if isinstance(value, DataSet):
            if value.shape != self._shape:
                raise TypeError("Invalid shape.")
            elif value.dtype != self.dtype:
                tmp = DataSet(self._shape, self._dtype)
                tmp.array[...] = value.array
                tmp.attrs.update(value.attrs)
                value = tmp
        else:
            tmp = np.asarray(value, self._dtype)
            if tmp.shape != self._shape:
                raise TypeError("Invalid shape.")
            value = DataSet(self._shape, self._dtype)
            value.array[...] = tmp
        return value

    def show2D(self, interpolation="nearest", cmap=None, vmin=None, vmax=None, conditioner=lambda x: x, q=0.03):
        """Shows series of 2D data in interactive matplotlib figure.
        conditioner(array) is an optional function called to condition the original data and returning the conditioned array.
        If vmin/vmax is None the min(q'th) / max((1-q)th) quantile are used."""
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider
        if len(self.shape) != 2:
            raise ValueError("No 2D data.")
        if cmap is None:
            cmap = plt.cm.gray
        if vmin is None or vmax is None:
            from scipy.stats import scoreatpercentile
            mn = float('+inf')
            mx = float('-inf')
            for i in range(self.indexSize):
                index = np.unravel_index(i, self.indexShape)
                array = self._data.__getitem__(index)
                if array is None:
                    continue
                array = conditioner(array[...]).ravel()
                if vmin is None:
                    mn = min(mn, scoreatpercentile(array, q * 100.0))
                if vmax is None:
                    mx = max(mx, scoreatpercentile(array, (1.0 - q) * 100.0))
            if not np.isfinite(mn):
                mn = 0.0
            if not np.isfinite(mx):
                mx = 0.0
            if vmin is None:
                vmin = mn
            if vmax is None:
                vmax = mx

        dim_scale = np.array(self.attrs.get("dim_scale", [1.0, 1.0]))
        dim_unit = self.attrs.get("dim_unit")
        dim_offset = np.array(self.attrs.get("dim_offset", [0.0, 0.0]))
        tl = np.array((-0.5, -0.5))
        br = np.array((self.shape[0] - 0.5, self.shape[1] - 0.5))
        if dim_scale.size == 1:
            tl *= np.asscalar(dim_scale)
            br *= np.asscalar(dim_scale)
        elif dim_scale.ndim == 1:
            tl *= dim_scale
            br *= dim_scale
        else:
            tl, br = np.dot(dim_scale, tl), np.dot(dim_scale, br)
        tl += dim_offset
        br += dim_offset

        fig = plt.figure()
        fig.add_subplot(111)
        plt.subplots_adjust(bottom=len(self.indexShape) * 0.05 + 0.15)
        if dim_unit is not None:
            plt.xlabel(dim_unit[0])
            plt.ylabel(dim_unit[1])
        array = self._data.__getitem__(np.unravel_index(0, self.indexShape))
        if array is None:
            array = np.zeros(self.shape, dtype=self.dtype)
        array = conditioner(array[...])
        im = plt.imshow(array, extent=[tl[0], br[0], br[1], tl[1]], vmax=vmax, vmin=vmin, interpolation=interpolation,
                        origin="upper", cmap=cmap)

        def updatePlot(xxx):
            index = tuple(int(slider.val) for slider in indexSlider)
            array = self._data.__getitem__(index)
            if array is None:
                array = np.zeros(self.shape, dtype=self.dtype)
            array = conditioner(array[...])
            im.set_array(array)
            fig.canvas.draw()

        indexSlider = []
        for i in range(len(self.indexShape)):
            axes = plt.axes([0.25, (len(self.indexShape) - i) * 0.05, 0.65, 0.03])
            slider = Slider(axes, 'Index%d' % i, 0.0, self.indexShape[i] - 1, valinit=0)
            slider.on_changed(updatePlot)
            indexSlider.append(slider)
        plt.show()

    def saveHDF5(self, fileName, groupName, mode="a", image=None):
        """Save series to a HDF5 file.

        Arguments:
            fileName
                Name of the HDF5 file
            groupName
                Name under which the series is saved.
            mode
                File update mode
            image
                Save elements of DataSeries as image (see saveHDF5)
        """
        import h5py
        from .hdf5 import saveHDF5, setAttrHDF5
        unsaved = []
        with h5py.File(fileName, mode) as fd:
            if groupName in fd:
                del fd[groupName]
            group = fd.create_group(groupName)
            for k, v in self.attrs.items():
                setAttrHDF5(group, unsaved, k, v)
            # Save dataset format
            setAttrHDF5(group, unsaved, '_shape', self.shape)
            setAttrHDF5(group, unsaved, '_indexShape', self.indexShape)
            setAttrHDF5(group, unsaved, '_dtype', self.dtype.str)
            # Save items
            for flatIndex in range(self.indexSize):
                index = np.unravel_index(flatIndex, self.indexShape)
                name = '_'.join('%03d' % i for i in index)
                dataSet = self.__getitem__(index)
                saveHDF5(group, dataSet.array, dataName=name, image=image, warn_attributes=False)
        if len(unsaved) != 0:
            print("The following series attributes weren't saved:")
            for item in unsaved:
                print("\t", item)


class Series(AbstractSeries):
    """Series is implemented as memory arrays."""

    def __init__(self, indexShape, shape, dtype=float):
        """Create empty series of given shape and data type.

        Arguments:
            indexShape : shape of array of datasets
            shape : Shape of datasets
            dtype : Type of datasets
        """
        AbstractSeries.__init__(self, indexShape, shape, dtype)
        self._data = np.zeros(indexShape, dtype=np.object)

    def __getitem__(self, index):
        return self._data.__getitem__(index)

    def __setitem__(self, index, value):
        value = self._prepareDataSet(value)
        self._data.__setitem__(index, value)

    @staticmethod
    def fromFiles(listOfNames, loader, verbose=0):
        """Creates a series from a number of files.

        Arguments
            listOfNames
                List of strings with filenames.
            loader
                Function that takes a string as argument and returns a DataSet as result.
        """
        # First dataset
        first = loader(listOfNames[0])
        if verbose > 0:
            print("Loading...")
            print("\t%d entities, shape=%s, dtype=%s" % (len(listOfNames), first.shape, first.dtype))
            print("\t[%02d] %s" % (0, listOfNames[0]))
        result = Series(len(listOfNames), first.shape, first.dtype)
        result[0] = first
        # Get known attributes
        keys = DataSet.general_attributes + DataSet.imaging_attributes
        keys.remove("dim_offset")  # Allow to vary over series
        for k in keys:
            v = first.attrs.get(k)
            if v is not None:
                result.attrs[k] = v
        # Get others
        for i in range(1, len(listOfNames)):
            if verbose > 0:
                print("\t[%02d] %s" % (i, listOfNames[i]))
            result[i] = loader(listOfNames[i])
        return result

    @staticmethod
    def fromHDF5(fileName, groupName):
        """Load series from a HDF5 file. The series must be stored by
        :meth:`saveHDF5` previously.

        Arguments:
            fileName
                Name of the HDF5 file
            groupName
                Name under which the series was saved.
        """
        import h5py
        from .hdf5 import loadAttrDictHDF5
        with h5py.File(fileName, "r") as fd:
            group = fd[groupName]
            attrs = loadAttrDictHDF5(group.attrs)
            # Get dataset format and create series
            shape = attrs['_shape']
            dtype = attrs['_dtype']
            indexShape = attrs['_indexShape']
            del attrs['_shape']
            del attrs['_dtype']
            del attrs['_indexShape']
            series = Series(indexShape, shape, dtype=dtype)
            series.attrs.update(attrs)
            # Load items
            for flatIndex in range(series.indexSize):
                index = np.unravel_index(flatIndex, series.indexShape)
                name = '_'.join('%03d' % i for i in index)
                savedSet = group[name]
                dataSet = DataSet(savedSet.shape, savedSet.dtype)
                dataSet.array[...] = savedSet[...]
                dataSet.attrs.update(loadAttrDictHDF5(savedSet.attrs))
                series[index] = dataSet
        return series


class LazyLoadingSeries(AbstractSeries):
    """Series only keeps a small number (capacity) of datasets in memory,
    and will load DataSets on demand. However, the individual
    datasets should be considered as read-only."""

    def __init__(self, fileList, loader, shape, dtype=float, capacity=3):
        """Create empty series of given shape and data type.

        Arguments:
            indexShape : shape of array of datasets
            shape : Shape of datasets
            dtype : Type of datasets
        """
        AbstractSeries.__init__(self, len(fileList), shape, dtype)
        self._loader = loader
        self._fileList = fileList
        self._capacity = 3
        self._cache = {}
        self._lru = []

    def __getitem__(self, index):
        import numbers
        if isinstance(index, numbers.Number):
            index = (index,)
        flatIndex = np.ravel_multi_index(index, self.indexShape)
        if flatIndex not in self._cache:
            data = self._load(flatIndex)
            self._insert(flatIndex, data)
        where = self._lru.index(flatIndex)
        del self._lru[where]
        self._lru.append(flatIndex)
        return self._cache[flatIndex]

    def __setitem__(self, index, value):
        raise TypeError("LazyLoaderSeries doesn't allow to set items.")

    def _load(self, flatIndex):
        data = self._loader(self._fileList[flatIndex])
        return data

    def _insert(self, flatIndex, value):
        if flatIndex in self._cache:
            raise RuntimeError("Index already in cache.")
        value = self._prepareDataSet(value)
        if len(self._cache) >= self._capacity:
            oldIndex = self._lru[0]
            del self._lru[0]
            del self._cache[oldIndex]
        self._cache[flatIndex] = value
        self._lru.append(flatIndex)

    @staticmethod
    def fromFiles(listOfNames, loader, verbose=0):
        """Creates a series from a number of files.

        Arguments
            listOfNames
                List of strings with filenames.
            loader
                Function that takes a string as argument and returns a DataSet as result.
        """
        # First dataset
        first = loader(listOfNames[0])
        result = LazyLoadingSeries(listOfNames, loader, first.shape, first.dtype)
        # Get known attributes
        keys = list(DataSet.general_attributes)
        keys.remove("dim_offset")  # Allow to vary over series
        for k in keys:
            v = first.attrs.get(k)
            if v is not None:
                result.attrs[k] = v
        # Insert first and done
        result._insert(0, first)
        return result
