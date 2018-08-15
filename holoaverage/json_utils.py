# This file is part of holoaverage.
# Copyright (c) 2018 Tore Niermann
#
# holoaverage is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# holoaverage is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with holoaverage.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import re
import uuid
import json


__all__ = ("encode_json", "decode_json")


def remove_comments(json_like):
    """
    Removes C-style comments from *json_like* and returns the result.  Example::

        >>> test_json = '''\
        {
            "foo": "bar", // This is a single-line comment
            "baz": "blah" /* Multi-line
            Comment */
        }'''
        >>> remove_comments('{"foo":"bar","baz":"blah",}')
        '{\n    "foo":"bar",\n    "baz":"blah"\n}'

    Original author:
        Dan McDougall <daniel.mcdougall@liftoffsoftware.com>

    License:
        "Unlicense"

    See:
        https://gist.github.com/liftoff/ee7b81659673eca23cd9fc0d8b8e68b7
    """
    comments_re = re.compile(r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
                             re.DOTALL | re.MULTILINE)
    def replacer(match):
        s = match.group(0)
        if s[0] == '/': return ""
        return s
    return comments_re.sub(replacer, json_like)


def remove_trailing_commas(json_like):
    """
    Removes trailing commas from *json_like* and returns the result.  Example::

        >>> remove_trailing_commas('{"foo":"bar","baz":["blah",],}')
        '{"foo":"bar","baz":["blah"]}'

    Original author:
        Dan McDougall <daniel.mcdougall@liftoffsoftware.com>

    License:
        "Unlicense"

    See:
        https://gist.github.com/liftoff/ee7b81659673eca23cd9fc0d8b8e68b7
    """
    trailing_object_commas_re = re.compile(r'(,)\s*}(?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    trailing_array_commas_re = re.compile(r'(,)\s*\](?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    # Fix objects {} first
    objects_fixed = trailing_object_commas_re.sub("}", json_like)
    # Now fix arrays/lists [] and return the result
    return trailing_array_commas_re.sub("]", objects_fixed)


class SmartJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, uuid.UUID):
            return str(obj)
        return json.JSONEncoder.default(self, obj)


def encode_json(value, return_bytes=False):
    """
    Encode python object to JSON string.

    `value` should be a dict, list, or tuple only containing JSON compatible primitives.
    Special handling is implemented for numpy `ndarray` objects and UUIDs.

    :param value: Value to encode.
    :param return_bytes: Return string as bytes object (otherwise it is a str object).
    :type return_bytes: bool
    """
    result = json.dumps(value, cls=SmartJSONEncoder)
    if return_bytes:
        result = result.encode("utf-8")
    return result


def decode_json(value, strict=False):
    """
    Decode JSON string to python object.

    :param value: Input string
    :param strict: If 'False' comments and trailing commas are allowed
    :return: Equivalent python object.
    """
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if not strict:
        value = remove_comments(value)
        value = remove_trailing_commas(value)
    return json.loads(value)
