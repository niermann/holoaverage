"""
An example of how to remove comments and trailing commas from JSON before
parsing.  You only need the two functions below, `remove_comments()` and
`remove_trailing_commas()` to accomplish this.

Original author:
    Dan McDougall <daniel.mcdougall@liftoffsoftware.com>

License:
    "Unlicense"

See:
    https://gist.github.com/liftoff/ee7b81659673eca23cd9fc0d8b8e68b7
"""

import re, fileinput

__all__ = ("remove_comments", "remove_trailing_commas")

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
    """
    trailing_object_commas_re = re.compile(r'(,)\s*}(?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    trailing_array_commas_re = re.compile(r'(,)\s*\](?=([^"\\]*(\\.|"([^"\\]*\\.)*[^"\\]*"))*[^"]*$)')
    # Fix objects {} first
    objects_fixed = trailing_object_commas_re.sub("}", json_like)
    # Now fix arrays/lists [] and return the result
    return trailing_array_commas_re.sub("]", objects_fixed)
