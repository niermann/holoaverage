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
import h5py

from .series import DataSet


def loadAttrDictHDF5(attrs):
    result = {}
    for k, v in attrs.items():
        # Parse attributes with slash into dictionaries
        subDict = result
        while True:
            index = k.find("/")
            if index < 0:
                break
            if index == 0:
                k = k[1:]
                continue
            subDict = subDict.setdefault(k[:index], dict())
            k = k[index + 1:]
        if len(k) > 0:
            subDict[k] = v
    return result


def loadHDF5(fileOrFileName, dataName=None):
    """Load dataset from HDF5 file.

    Arguments:
        fileOrFilename
            Name of file to load or HDF5 file.
        dataName
            Name of dataset to load. If None dataset name should be appended to filename (separated by question mark)

    Returns:
        data as :class:`DataSet`
    """
    if dataName is None:
        items = fileOrFileName.split('?', 1)
        fileOrFileName = items[0]
        dataName = items[1]
    if isinstance(fileOrFileName, h5py.File):
        dataset = fileOrFileName[dataName]
        result = DataSet(dataset.shape, dataset.dtype)
        result.array[...] = dataset[...]
        result.attrs.update(loadAttrDictHDF5(dataset.attrs))
    else:
        with h5py.File(fileOrFileName, "r") as fd:
            dataset = fd[dataName]
            result = DataSet(dataset.shape, dataset.dtype)
            result.array[...] = dataset[...]
            result.attrs.update(loadAttrDictHDF5(dataset.attrs))
    return result


def setAttrHDF5(dataset, unsaved, key, value, prefix=""):
    """Unravel dictionary into slashed pathes"""
    if value is None:
        return
    elif isinstance(value, dict):
        prefix = key + "/"
        for k, v in value.items():
            setAttrHDF5(dataset, unsaved, prefix + k, v)
    else:
        try:
            value = np.asarray(value)
            if value.dtype.kind == 'U':
                value = np.char.encode(value, "utf8")
                dataset.attrs[key] = value
            elif value.dtype.kind == 'O':
                # Save object arrays as numbered keys...
                prefix = key + "/"
                for n, v in enumerate(value):
                    setAttrHDF5(dataset, unsaved, prefix + str(n), v)
            else:
                # Save object lists or array
                dataset.attrs[key] = value
        except TypeError as exc:
            unsaved.append("%s: %s" % (key, str(exc)))


def setH5Image(dataset, vmin, vmax):
    if vmin is None:
        vmin = np.amin(dataset)
    if vmax is None:
        vmax = np.amax(dataset)
    vrange = (vmin, vmax)
    dataset.attrs[b'CLASS'] = np.string_('IMAGE')
    dataset.attrs[b'IMAGE_SUBCLASS'] = np.string_('IMAGE_GRAYSCALE')
    dataset.attrs[b'IMAGE_WHITE_IS_ZERO'] = np.uint8(0)
    dataset.attrs[b'DISPLAY_ORIGIN'] = np.string_('UL')
    dataset.attrs[b'IMAGE_MINMAXRANGE'] = np.array(vrange, dtype=dataset.dtype)


def saveHDF5(groupOrFileName, data, dataName='data', mode="a", image=None, vmin=None, vmax=None, warn_attributes=True):
    """Load dataset from HDF5 file.

    Arguments:
        groupOrFileName
            Name of file to load or HDF5 file/group.
        data
            Data to save
        dataName
            Name of dataset to save.
        mode
            Opening mode of file
        image
            Save as image, by default 2D datasets of integer and float types are saved as image
        vmin, vmax
            Color range if image=True
        warn_attributes
            Warn if there are unsaved attributes
    Returns:
        data as :class:`DataSet`
    """
    if isinstance(data, DataSet):
        attrs = data.attrs
        data = data.array[...]
    else:
        attrs = {}
    if image is None and data.ndim == 2 and data.dtype.kind in ['b', 'i', 'u', 'f']:
        image = True
    if isinstance(groupOrFileName, h5py.Group):
        if dataName in groupOrFileName:
            del groupOrFileName[dataName]
        dataset = groupOrFileName.create_dataset(dataName, data=data)
        unsaved = []
        for k, v in attrs.items():
            setAttrHDF5(dataset, unsaved, k, v)
        if image:
            setH5Image(dataset, vmin, vmax)
    else:
        with h5py.File(groupOrFileName, mode) as fd:
            if dataName in fd:
                del fd[dataName]
            dataset = fd.create_dataset(dataName, data=data)
            unsaved = []
            for k, v in attrs.items():
                setAttrHDF5(dataset, unsaved, k, v)
            if image:
                setH5Image(dataset, vmin, vmax)
    if len(unsaved) != 0 and warn_attributes:
        print("The following attributes weren't saved:")
        for item in unsaved:
            print("\t", item)
