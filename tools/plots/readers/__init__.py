"""Readers of the files spodvr writes (layer 1): one function per file kind, returning a dataclass of numpy
arrays. Readers never print or draw; they keep NaN data and record the first non-finite record
(`first_non_finite`), so that a check, not the reader, decides what fails. Each module names the C++
writer it follows in a `Shadows:` line."""
