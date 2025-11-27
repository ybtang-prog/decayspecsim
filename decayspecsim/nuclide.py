"""
Nuclide and decay mode classes
"""

from typing import List, Dict, Optional
import numpy as np


class DecayMode:
    """Base class for decay modes"""
    
    def __init__(self, data: Dict):
        self.data = data
        self.lines = data.get('lines', {})
        self.energies = np.array(self.lines.get('energies', []))
        self.intensities = np.array(self.lines.get('intensities', []))
        self.norms = np.array(self.lines.get('norms', []))
        self.number = data.get('number', 0)
        self.mean_energy = data.get('mean_energy', 0.0)
    
    def get_effective_intensities(self) -> np.ndarray:
        """Calculate effective intensities (intensity * norm)"""
        return self.intensities * self.norms


class GammaDecay(DecayMode):
    """Gamma ray decay mode"""
    pass


class BetaDecay(DecayMode):
    """Beta decay mode - continuous spectrum"""
    pass


class AlphaDecay(DecayMode):
    """Alpha particle decay mode"""
    pass


class ElectronDecay(DecayMode):
    """Electron capture/internal conversion decay mode"""
    pass


class XRayDecay(DecayMode):
    """X-ray emission decay mode"""
    pass


class NeutronDecay(DecayMode):
    """Neutron emission decay mode"""
    pass


class Nuclide:
    """Represents a radioactive nuclide with its decay modes"""
    
    def __init__(self, name: str, data: Dict):
        self.name = name
        self.zai = data.get('zai')
        self.halflife = data.get('halflife', 0.0)  # in seconds
        self.decay_modes = {}
        
        # Map decay mode types
        mode_map = {
            'gamma': GammaDecay,
            'beta': BetaDecay,
            'alpha': AlphaDecay,
            'electron': ElectronDecay,
            'x-ray': XRayDecay,
            'neutron': NeutronDecay,
        }
        
        for mode_type, mode_class in mode_map.items():
            if mode_type in data:
                self.decay_modes[mode_type] = mode_class(data[mode_type])
    
    def get_decay_mode(self, mode_type: str) -> Optional[DecayMode]:
        """Get decay mode by type"""
        return self.decay_modes.get(mode_type)
    
    def has_decay_mode(self, mode_type: str) -> bool:
        """Check if nuclide has specific decay mode"""
        return mode_type in self.decay_modes
    
    def list_decay_modes(self) -> List[str]:
        """List all available decay modes"""
        return list(self.decay_modes.keys())

