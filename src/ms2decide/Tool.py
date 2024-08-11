from enum import Enum


class Tool(Enum):
    """
    Enumeration representing different tools for analysis.

    This enum defines constants for various tools used in chemical analysis.

    Attributes:
        GNPS (int): Constant representing the GNPS tool (value: 1).
        SIRIUS (int): Constant representing the SIRIUS tool (value: 2).
        ISDB (int): Constant representing the ISDB tool (value: 3).
    """
    GNPS = 1
    SIRIUS = 2
    ISDB = 3
