"""
DecaySpecSim - Nuclear decay spectrum simulator based on Decay2012 database
"""

__version__ = "1.0.0"

from .database import Database
from .nuclide import Nuclide, DecayMode, GammaDecay, BetaDecay, AlphaDecay, ElectronDecay, XRayDecay
from .spectrum import SpectrumGenerator
from .detector import Detector
from .decay_chain import DecayChainSimulator

__all__ = [
    'Database',
    'Nuclide',
    'DecayMode',
    'GammaDecay',
    'BetaDecay',
    'AlphaDecay',
    'ElectronDecay',
    'XRayDecay',
    'SpectrumGenerator',
    'Detector',
    'DecayChainSimulator',
]

