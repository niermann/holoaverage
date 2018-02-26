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

"""
DM3 file loader

For details on the format see:
    Chris Boothroyd's site at http://www.er-c.org/cbb/info/dmformat/
    Greg Jefferis' size at http://rsb.info.nih.gov/ij/plugins/DM3Format.gj.html
"""

import codecs
import struct
import numpy as np

from .series import DataSet


def _get_byteorder():
    import sys
    if sys.byteorder == 'little':
        return '<'
    elif sys.byteorder == 'big':
        return '>'
    raise ValueError("Unknown byte order.")

_host_byteorder = _get_byteorder()

_tagname_errorhandler = 'ignore'
_tagname_decoder = codecs.getdecoder('latin_1')  # just a guess

_latin1_decoder = codecs.getdecoder('latin_1')
_utf16_le_decoder = codecs.getdecoder("'utf_16_le")
_utf16_be_decoder = codecs.getdecoder("'utf_16_be")

class TagBase(object):
    """
    Base class for all DM3 tags.

    Depending on :attr:`tag_type`, the tag is either a :class:`TagItem` (:const:`TAGTYPE_ITEM`)
    or a :class:`TagContainer` (:const:`TAGTYPE_GROUP`). The former contain
    data tags, the latter other tags.
    """
    TAGTYPE_GROUP = 0x14
    TAGTYPE_ITEM = 0x15

    def __init__(self, name, raw, offset, byteorder=_host_byteorder):
        """
        :param name: Name of tag
        :type name: str
        :param raw: Tag in raw form.
        :type raw: bytes
        :param offset: Offset of first byte of tag
        :type offset: int
        :param byteorder: Byteroder of the tag
        :type byteorder: '<' or '>'
        """
        self._raw = raw
        self._offset = offset
        self._byteorder = byteorder
        self._length = 0
        self._name = name
        #print(hex(self._offset), self.tag_type, self._name)

    @property
    def raw(self):
        """Raw bytes representation of tag."""
        return self._raw[self._offset:self._offset + self._length]

    @property
    def name(self):
        """Name of tag."""
        return self._name

    @property
    def byteorder(self):
        """Byte order of the tag. Either "<" for little endian or ">" for big endian."""
        return self._byteorder

    @property
    def tag_type(self):
        """Type of tag."""
        raise NotImplementedError

    def as_object(self):
        """Return tag converted to python object."""
        raise NotImplementedError


class TagItem(TagBase):
    """
    Class representing an entry in a DM3 tag group.

    TagItem instances contain data (opposed to other tags). The type of the
    data is given by :attr:`item_type` and can be further identified by the
    is_xxx group of methods.

    To retrieve the data in a specific object type (number, array, ...) use
    the corresponding as_xxx methods. These will fail if the item type is
    not compatible to the destination type. The :meth:`as_object` method
    will also choose a corresponding python object type for returning the data.
    """

    # Tag entry types
    ITEMTYPE_INT16 = 0x02
    ITEMTYPE_INT32 = 0x03
    ITEMTYPE_UINT16 = 0x04
    ITEMTYPE_UINT32 = 0x05
    ITEMTYPE_FLOAT32 = 0x06
    ITEMTYPE_FLOAT64 = 0x07
    ITEMTYPE_BOOL = 0x08
    ITEMTYPE_INT8 = 0x09
    ITEMTYPE_UINT8 = 0x0a
    ITEMTYPE_INT64 = 0x0b        # ???
    ITEMTYPE_UINT64 = 0x0c       # ???
    ITEMTYPE_RECORD = 0x0f
    ITEMTYPE_STRING = 0x12       # ???
    ITEMTYPE_ARRAY = 0x14

    # Magic value
    ITEM_MAGIC = 0x25252525

    # Number format info: (struct-format, dtype, size in bytes)
    _NUMBER_INFO = {
        ITEMTYPE_INT8:    ('b', 'i1', 1),
        ITEMTYPE_INT16:   ('h', 'i2', 2),
        ITEMTYPE_INT32:   ('i', 'i4', 4),
        ITEMTYPE_INT64:   ('q', 'i8', 8),
        ITEMTYPE_UINT8:   ('B', 'u1', 1),
        ITEMTYPE_UINT16:  ('H', 'u2', 2),
        ITEMTYPE_UINT32:  ('I', 'u4', 4),
        ITEMTYPE_UINT64:  ('Q', 'u8', 8),
        ITEMTYPE_FLOAT32: ('f', 'f4', 4),
        ITEMTYPE_FLOAT64: ('d', 'f8', 8),
        ITEMTYPE_BOOL:    ('?', 'b', 1),
    }

    def __init__(self, raw, offset, byteorder=_host_byteorder):
        """
        :param raw: Tag in raw form.
        :type raw: bytes
        :param offset: Offset of first byte of tag
        :type offset: int
        :param byteorder: Byteroder of the tag
        :type byteorder: '<' or '>'
        """
        name_length = struct.unpack('>H', raw[offset:offset + 2])[0]
        name = _tagname_decoder(raw[offset + 2:offset + 2 + name_length], _tagname_errorhandler)[0]
        super(TagItem, self).__init__(name, raw, offset, byteorder)
        magic = struct.unpack('>I', raw[offset + 2 + name_length:offset + 6 + name_length])[0]
        if magic != TagItem.ITEM_MAGIC:
            raise ValueError("Not a valid tag (magic number missing)")
        self._info_offset = offset + 10 + name_length
        self._info_num = struct.unpack('>I', raw[self._info_offset - 4:self._info_offset])[0]
        item_type = self.item_type
        if TagItem.is_number_type(item_type):
            data_length = TagItem._number_info(item_type)[2]
        elif item_type == TagItem.ITEMTYPE_RECORD:
            data_length = self._get_record_format()[1]
        elif item_type == TagItem.ITEMTYPE_ARRAY:
            param = self._get_array_param()
            data_length = param[0].itemsize * param[1]
        else:
            raise ValueError("Unknown item type.")
        self._length = data_length + self._info_offset - self._offset + 4 * self._info_num

    @property
    def tag_type(self):
        """Type of tag."""
        return TagBase.TAGTYPE_ITEM

    @property
    def item_type(self):
        """Type of this tag."""
        return self.info_field(0)

    def info_field(self, n):
        """Read n-th info field."""
        n = int(n)
        if n < 0 or n >= self._info_num:
            raise IndexError("Invalid info field index.")
        offset = self._info_offset + 4 * n
        return struct.unpack('>I', self._raw[offset:offset + 4])[0]

    @property
    def raw_item(self):
        """Return raw bytes form of item."""
        return self._raw[self._info_offset + self._info_num * 4:self._offset + self._length]

    @staticmethod
    def is_number_type(item_type):
        return TagItem.ITEMTYPE_INT16 <= item_type <= TagItem.ITEMTYPE_UINT64

    def is_number(self):
        """Return whether this is a number type, booleans are considered as numbers also."""
        return TagItem.is_number_type(self.item_type)

    @staticmethod
    def _number_info(item_type):
        if not TagItem.is_number_type(item_type):
            raise ValueError("Expected number item type in record.")
        return TagItem._NUMBER_INFO[item_type]

    def _parse_number(self, item_type):
        """
        Return item as python number.

        :param tag_type: Type
        :type tag_type: ITEMTYPE_...
        :raises TypeError: If not a number (or bool).
        :return: Parsed value
        """
        if not TagItem.is_number_type(self.item_type):
            raise TypeError("Not a number item.")
        return struct.unpack(self._byteorder + TagItem._number_info(self.item_type)[0], self.raw_item)[0]

    def as_number(self):
        """
        Return item as python number.

        :raises TypeError: If not a number (or bool).
        """
        item_type = self.item_type
        value = self._parse_number(item_type)
        if item_type == TagItem.ITEMTYPE_BOOL:
            return int(value)
        return value

    def as_int(self):
        """Like :meth:`as_number` but casts result to :class:`int`."""
        return int(self.as_number())

    def as_float(self):
        """Like :meth:`as_number` but casts result to :class:`float`."""
        return float(self.as_number())

    def is_bool(self):
        """Whether this is a boolean tag."""
        return TagItem.ITEMTYPE_BOOL == self.item_type

    def as_bool(self):
        """
        Return item as boolean.

        :raises TypeError: If not a bool (or number).
        """
        item_type = self.item_type
        value = self._parse_number(item_type)
        if item_type != TagItem.ITEMTYPE_BOOL:
            return bool(value)
        return value

    def is_string(self):
        """
        Return whether this is a string type.

        Valid item types, are the (unsupported) ITEMTYPE_STRING,
        uint8 arrays (interpreted as latin-1 encoding) and uint16 arrays
        (interpreted as utf16 encoding).
        """
        item_type = self.item_type
        if item_type == TagItem.ITEMTYPE_STRING:
            return True
        elif item_type == TagItem.ITEMTYPE_ARRAY:
            base_type = self.info_field(1)
            return base_type == TagItem.ITEMTYPE_UINT8 or base_type == TagItem.ITEMTYPE_UINT16
        else:
            return False

    def as_string(self):
        """
        Return item as string (see :meth:`is_string` for how strings are interpreted).

        :raises TypeError: If not a string type.
        :return: String
        """
        item_type = self.item_type
        if item_type == TagItem.ITEMTYPE_ARRAY:
            base_type = self.info_field(1)
            if base_type == TagItem.ITEMTYPE_UINT8:
                return _latin1_decoder(self.raw_item)[0]
            elif base_type == TagItem.ITEMTYPE_UINT16:
                if self._byteorder == '<':
                    return _utf16_le_decoder(self.raw_item)[0]
                elif self._byteorder == '>':
                    return _utf16_be_decoder(self.raw_item)[0]
        raise TypeError("Invalid tag type.")

    def _get_record_format(self):
        """
        Return struct format string for a record and its size.

        :returns: (struct string, length)
        :rtype: (str, int)
        :raises TypeError: If not a record.
        """
        if not self.is_record():
            raise TypeError("Not a record item.")
        # Ignore record name and field names
        field_num = self.info_field(2)
        if field_num <= 0:
            raise ValueError("Empty record.")
        fmt = self._byteorder
        size = 0
        for n in range(field_num):
            if self.info_field(3 + 2*n) != 0:
                raise ValueError("Unexpected non-zero value in information field.")
            field_type = self.info_field(4 + 2 * n)
            info = TagItem._number_info(field_type)
            fmt += info[0]
            size += info[2]
        return fmt, size

    def is_record(self):
        """Return whether this is a record item."""
        return self.item_type == TagItem.ITEMTYPE_RECORD

    def as_record(self):
        """
        Return item as tuple.

        :raises TypeError: If not a record.
        """
        fmt = self._get_record_format()[0]
        return struct.unpack(fmt, self.raw_item)

    def is_group(self):
        """Return whether this is a group (dict) type."""
        return False

    def is_array(self):
        """Return whether this is an array item."""
        return self.item_type == TagItem.ITEMTYPE_ARRAY

    def _dtype_string(self, item_type):
        if not TagItem.is_number_type(item_type):
            raise ValueError("Expected number item type in record.")
        return self._byteorder + TagItem._NUMBER_INFO[item_type][1]

    def _get_array_param(self):
        """
        Return dtype and length of array.

        :returns: (numpy dtype, length)
        :rtype: (dtype, int)
        :raises TypeError: If not an array.
        """
        if not self.is_array():
            raise TypeError("Not an array item.")
        base_type = self.info_field(1)
        if TagItem.is_number_type(base_type):
            dtype = np.dtype(self._dtype_string(base_type))
            length = self.info_field(2)
        elif base_type == TagItem.ITEMTYPE_RECORD:
            # Ignore record name and field names, however check that they are zero
            if self.info_field(4) != 0:
                raise ValueError("Unexpected non-zero value in information field.")
            field_num = self.info_field(3)
            if field_num <= 0:
                raise ValueError("Empty record.")
            dtype = ''
            for n in range(field_num):
                if self.info_field(4 + 2*n) != 0:
                    raise ValueError("Unexpected non-zero value in information field.")
                dtype += ", " + self._dtype_string(self.info_field(5 + 2 * n))
            dtype = dtype[2:]
            length = self.info_field(4 + 2 * field_num)
        else:
            raise ValueError("Unexpected array base type.")
        return np.dtype(dtype), length

    def get_array_dtype(self):
        """
        Return dtype of array,

        :raises TypeError: If not an array.
        """
        return self._get_array_param()[0]

    def get_array_size(self):
        """
        Return number of items in array,

        :raises TypeError: If not an array.
        """
        return self._get_array_param()[1]

    def as_array(self):
        """
        Return item as numpy array.

        :raises TypeError: If not an array.
        """
        param = self._get_array_param()
        return np.frombuffer(self._raw, dtype=param[0], count=param[1], offset=self._info_offset + self._info_num * 4)

    def as_object(self):
        item_type = self.item_type
        if TagItem.is_number_type(item_type):
            return self.as_number()
        elif item_type == TagItem.ITEMTYPE_RECORD:
            return self.as_record()
        elif item_type == TagItem.ITEMTYPE_ARRAY:
            return self.as_array()
        else:
            raise TypeError("Unknown item type.")
    as_object.__doc__ = TagBase.as_object.__doc__


class TagContainer(TagBase):
    """
    TagContainer instances are tags that contain other tags (opposed to
    actual data).

    TagContainers can be indexed like lists using the [] operator with
    an argument castable to int. The number of tags inside can be retrieved
    using len() and the container can be iterated.

    Named tags can also be retrieved using string arguments to the [] operator.
    Further more tag paths of the form "some_field:5:other" refering to
    successive item retrievals ["some_field"][5]["other"] can also be
    passed directly into the [] operator: ["some_field:5:other"]
    """
    def __init__(self, name, raw, offset, byteorder):
        """
        :param name: Name of tag
        :type name: str
        :param raw: Tag in raw form.
        :type raw: bytes
        :param offset: Offset of first byte of tag
        :type offset: int
        :param length: Length of tag in bytes
        :type length: int
        :param byteorder: Byteroder of the tag
        :type byteorder: '<' or '>'
        """
        super(TagContainer, self).__init__(name, raw, offset, byteorder)
        self._entry_list = []
        self._entry_map = {}

    def _init_container(self, offset):
        # Ignore sorted, open flags
        entry_num = struct.unpack('>I', self._raw[offset + 2:offset + 6])[0]
        data_offset = offset + 6
        for n in range(entry_num):
            tag_type, = struct.unpack('B', self._raw[data_offset:data_offset+1])
            if tag_type == TagItem.TAGTYPE_ITEM:
                entry = TagItem(self._raw, data_offset + 1, self._byteorder)
            elif tag_type == TagItem.TAGTYPE_GROUP:
                entry = TagGroup(self._raw, data_offset + 1, self._byteorder)
            else:
                raise ValueError("Unknown tag type.")
            self._entry_list.append(entry)
            if entry.name:
                self._entry_map[entry.name] = entry
            data_offset += 1 + len(entry.raw)
        self._length = data_offset - self._offset

    @property
    def tag_type(self):
        """Type of tag."""
        return TagBase.TAGTYPE_GROUP

    def as_object(self):
        if len(self._entry_list) > 0 and len(self._entry_list) == len(self._entry_map):
            return dict(self._entry_map)
        return list(self._entry_list)
    as_object.__doc__ = TagBase.as_object.__doc__

    def __len__(self):
        return len(self._entry_list)

    def _get_item(self, key):
        try:
            index = int(key)
        except:
            if not key:
                raise KeyError("Invalid key")
            return self._entry_map[key]
        else:
            return self._entry_list[index]

    def __getitem__(self, key):
        try:
            pos = key.index(":")
        except:
            return self._get_item(key)
        else:
            return self._get_item(key[0:pos])[key[pos + 1:]]

    def _contains_item(self, key):
        try:
            index = int(key)
        except:
            return key in self._entry_map
        else:
            return 0 <= index < len(self._entry_list)

    def __contains__(self, key):
        try:
            pos = key.index(":")
        except:
            return self._contains_item(key)
        else:
            try:
                group = self._get_item(key[0:pos])
            except:
                return False
            else:
                return key[pos + 1:] in group


class TagGroup(TagContainer):
    """A group of tags."""
    def __init__(self, raw, offset, byteorder=_host_byteorder):
        """
        :param raw: Tag in raw form.
        :type raw: bytes
        :param offset: Offset of first byte of tag
        :type offset: int
        :param byteorder: Byteroder of the tag
        :type byteorder: '<' or '>'
        """
        name_length = struct.unpack('>H', raw[offset:offset + 2])[0]
        name = _tagname_decoder(raw[offset + 2:offset + 2 + name_length], _tagname_errorhandler)[0]
        super(TagGroup, self).__init__(name, raw, offset, byteorder)
        self._init_container(offset + 2 + name_length)


class DM3File(TagContainer):
    """
    A DM3 file containing tags.

    Instances can be used like TagContainers. They also implement the context
    manager protocol.

    """
    def __init__(self, filename, mode="r"):
        """
        :param filename: Filename.
        :type filename: str
        :param mode: Opening mode
        :type mode: "r"
        """
        if mode != 'r':
            raise ValueError("Unsupported file mode.")
        with open(filename, "rb") as fd:
            raw = fd.read()
        version, data_length, le_order = struct.unpack('>III', raw[0:12])
        if version != 3:
            raise ValueError("Not a DM3 file.")
        if le_order != 0:
            byteorder = '<'
        else:
            byteorder = '>'
        super(DM3File, self).__init__('', raw, 0, byteorder)
        self._init_container(12)
        self._filename = filename

    def close(self):
        """Close the associated file."""
        pass

    @property
    def filename(self):
        """Name of the file."""
        return self._filename

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def parse_tag(tag, use_string_heuristic=True):
    """
    Parse DM3 tag to full python object.

    This methods recursively converts the tag into python objects. This
    differentiates this method from the :meth:`TagBase.as_object` method
    which, when called on :class:`TagContainer` instances will return
    the container as `dict` or `list` consisting of :class:`TagBase` instances.

    As strings are simply stored as uint16 array within the tag groups, there
    is no way to distinguish them from array data. If `use_string_heuristic` is
    true also uint16 arrays are decoded as unicode string, otherwise they are
    returned as numpy arrays.

    :param tag: Tag to parse to python objects
    :type tag: TagBase
    :param use_string_heuristic: Convert 2 byte arrays to strings
    :type use_string_heuristic: bool
    :return: list or dict
    """
    if tag.tag_type == TagBase.TAGTYPE_ITEM:
        if use_string_heuristic and tag.item_type == TagItem.ITEMTYPE_ARRAY:
            try:
                if tag.get_array_dtype() == np.dtype('<u2'):
                    return tag.as_string()
                elif tag.get_array_dtype() == np.dtype('>u2'):
                    return tag.as_string()
                # Fall through in case of other data types
            except:
                # Fall through in case of error
                pass
        return tag.as_object()
    elif tag.tag_type == TagItem.TAGTYPE_GROUP:
        container = tag.as_object()
        if isinstance(container, dict):
            result = {}
            for key, value in container.items():
                result[key] = parse_tag(value, use_string_heuristic)
        else:
            result = []
            for value in container:
                result.append(parse_tag(value, use_string_heuristic))
        return result
    else:
        return tag.as_object()


def parse_known_metadata_tags(tags):
    """
    Parse DM3 ImageTags to DataSet metadata for known keys.

    :param tags: image tags of DM3 image
    :type tags: TagGroup
    :return: Dict with metadata
    """
    result = {}

    if 'Microscope Info:Voltage' in tags:
        result["voltage(kV)"] = tags['Microscope Info:Voltage'].as_float() / 1000.0
    if 'Acquisition:Parameters:High Level:Exposure (s)' in tags:
        result["exposure(s)"] = tags['Acquisition:Parameters:High Level:Exposure (s)'].as_float()
    if 'Acquisition:Parameters:High Level:Binning' in tags:
        result["binning"] = tags['Acquisition:Parameters:High Level:Binning'].as_object()
    if 'Microscope Info:Microscope' in tags:
        result["microscope"] = tags['Microscope Info:Microscope'].as_string()
    elif 'Tecnai:Microscope Info' in tags:     # FEI microscopes only
        result["microscope"] = tags['Tecnai:Microscope Info'].as_string().split(u'\u2028')[0]
    if 'Acquisition:Device:Name' in tags:
        result["detector"] = tags['Acquisition:Device:Name'].as_string()

    # Timestamp
    if ('DataBar:Acquisition Date' in tags) and ('DataBar:Acquisition Time' in tags):
        date_string = tags['DataBar:Acquisition Date'].as_string()
        time_string = tags['DataBar:Acquisition Time'].as_string().split(' ')
        import datetime
        date = datetime.datetime.strptime(date_string, '%m/%d/%Y').date()
        time = datetime.datetime.strptime(time_string[0], '%I:%M:%S').time()
        timestamp = datetime.datetime.combine(date, time)
        if time_string[1].upper() == 'PM':
            timestamp += datetime.timedelta(hours=12)
        result['timestamp'] = timestamp.isoformat('T')

    return result


_DM3_IMAGE_DATATYPE = {
    1:  'i2',
    2:  'f4',
    3:  'c8',
    6:  'u8',
    7:  'i4',
    9:  'i1',
    10: 'u2',
    11: 'u4',
    12: 'f8',
    13: 'c8',
    39: 'i8',
    40: 'u8'

    # Further (currently unsupported) types
    # 4: packed complex (result of r2h ffts)
    # 8: XRGB quads
    # 14: binary
}


def load_dm3(filename, index=1, include_tags=True):
    """
    Load :class:`DataSet` from DM3 file.

    DM3 image files can contain several images. In all observed files an index
    of 0 corresponds to the thumbnail and an index of 1 corresponds to the data
    itself.

    If `include_tags` is true, the DM3 ImageTags will be parsed into Python
    objects and included as 'dm3-tags' entry in the DataSet's metadata.

    :param filename: Name of the file to load_cel.
    :type filename: str
    :param index: Index of data in file (usually 1)
    :type index: int
    :param include_tags: Whether image tags should be parsed
    :type include_tags: bool
    :rtype: DataSet
    """
    with DM3File(filename, 'r') as tags:
        image = tags['ImageList'][index]
        image_data = image['ImageData']

        # Get data shape
        shape = ()
        for entry in image_data['Dimensions']:
            dim = entry.as_int()
            if dim <= 0:
                raise ValueError("Unexpected non-positive value for image extent.")
            shape = (dim,) + shape

        # Get data type
        image_type = image_data['DataType'].as_int()
        try:
            dtype = _DM3_IMAGE_DATATYPE[image_type]
        except KeyError:
            raise ValueError("Unsupported DM3 image data type.")

        # Read image
        image_array = image_data['Data']
        dtype = np.dtype(image_array.byteorder + dtype)
        dataset = DataSet(shape, dtype)
        dataset.array[...] = np.frombuffer(image_array.raw_item, dtype=dtype).reshape(shape)

        # Get dimension calibrations
        dim_scale = []
        dim_offset = []
        dim_unit = []
        for n, cal in enumerate(image_data['Calibrations:Dimension']):
            scale = cal['Scale'].as_float()
            dim_scale.append(scale)
            dim_offset.append(scale * cal['Origin'].as_float())
            dim_unit.append(cal['Units'].as_string())
        dataset.attrs["dim_unit"] = dim_unit
        dataset.attrs["dim_scale"] = np.array(dim_scale, dtype=float)
        dataset.attrs["dim_offset"] = np.array(dim_offset, dtype=float)

        # Get brightness calibrations
        scale = image_data['Calibrations:Brightness:Scale'].as_float()
        dataset.attrs['unit'] = image_data['Calibrations:Brightness:Units'].as_string()
        dataset.attrs['scale'] = scale
        dataset.attrs['offset'] = scale * image_data['Calibrations:Brightness:Origin'].as_float()

        # Handle image tags
        image_tags = image['ImageTags']
        dataset.attrs.update(parse_known_metadata_tags(image_tags))
        if include_tags:
            dataset.attrs['dm3_tags'] = parse_tag(image_tags)
            dataset.attrs['dm3_unique_id'] = parse_tag(image['UniqueID'])

        # More metadata
        dataset.attrs['source_file'] = filename
        dataset.attrs['source_name'] = image['Name'].as_string()

        return dataset