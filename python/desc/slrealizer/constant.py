"""
All the constants should be accessed through the methods here, so that they would not be confused.
"""

from astropy.coordinates import SkyCoord

def return_zeropoint():
    """
    returns zeropoint for the SDSS magnitude system
    """
    return 22.5

def return_obs_RA_DEC():
    """
    Return OM10's lensed system's position
    """
    return SkyCoord('03h 32m 30s', '10d 00m 24s')
