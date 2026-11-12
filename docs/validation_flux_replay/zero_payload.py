#!/usr/bin/env python3
"""Copy a gfortran sequential unformatted file with every record's payload set to zero.

The 4-byte record markers are kept, so the copy has exactly the layout of the original: the Fortran
reads it back as the same records, holding zeros (the zero-filled flux coefficient files of the
flux replay, docs/validation_flux_replay/README.md). Subrecords (records over 2 GiB, negative
markers) are refused.

Usage: python3 -I -B zero_payload.py IN OUT   (OUT is overwritten)
Output: OUT, and one line "OUT: N records zeroed" on stdout.
Environment: none read. Run time: well under a second for the replay's 1.3 MB files.
"""
import struct
import sys


def zero_payload(data: bytearray) -> int:
    """Zeroes the payload of every record of data in place; returns the record count."""
    position = 0
    records = 0
    while position < len(data):
        if position + 4 > len(data):
            raise ValueError(f"truncated record marker at byte {position}")
        (length,) = struct.unpack_from("<i", data, position)
        if length < 0:
            raise ValueError(f"record {records + 1} has subrecords (negative marker), which are not supported")
        end = position + 4 + length
        if end + 4 > len(data):
            raise ValueError(f"record {records + 1} runs past the end of the file")
        (tail,) = struct.unpack_from("<i", data, end)
        if tail != length:
            raise ValueError(f"record {records + 1} has the markers {length} and {tail}")
        data[position + 4:end] = bytes(length)
        position = end + 4
        records += 1
    return records


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: zero_payload.py IN OUT", file=sys.stderr)
        return 2
    try:
        with open(sys.argv[1], "rb") as stream:
            data = bytearray(stream.read())
        records = zero_payload(data)
    except (OSError, ValueError) as error:
        print(f"error: {sys.argv[1]}: {error}; give a gfortran sequential unformatted file", file=sys.stderr)
        return 1
    with open(sys.argv[2], "wb") as stream:
        stream.write(data)
    print(f"{sys.argv[2]}: {records} records zeroed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
