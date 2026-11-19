"""The one error the tool raises for a run it cannot read (or an output it may not write)."""


class ReadError(Exception):
    """A file that cannot be read, or written without replacing a file the tool did not write; the message
    says "<file>: <what>; <fix>"."""
