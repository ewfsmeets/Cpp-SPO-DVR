"""Readers of the files spodvr writes (layer 1): one function per file kind, returning a dataclass of numpy
arrays. Readers never print or draw; they keep NaN data, and the readers of time series note the first
non-finite record (`first_non_finite`), so that a check, not the reader, decides what fails. Each module
that reads a file kind names the C++
writer it follows in a `Shadows:` line."""
