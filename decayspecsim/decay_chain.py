"""
Decay chain and time evolution simulation
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from .nuclide import Nuclide
from .spectrum import SpectrumGenerator


class DecayChainSimulator:
    """Simulate decay chains and time evolution"""
    
    def __init__(self, database):
        self.database = database
        self.decay_constants = {}  # Cache for lambda values
    
    def get_decay_constant(self, nuclide: Nuclide) -> float:
        """Calculate decay constant from half-life"""
        if nuclide.halflife <= 0:
            return 0.0
        return np.log(2) / nuclide.halflife
    
    def simulate_single_decay(self, nuclide: Nuclide, 
                             initial_activity: float,
                             time: float) -> float:
        """
        Calculate activity at time t for single nuclide
        
        Args:
            nuclide: Nuclide object
            initial_activity: Initial activity (Bq)
            time: Time in seconds
        
        Returns:
            Activity at time t
        """
        lam = self.get_decay_constant(nuclide)
        return initial_activity * np.exp(-lam * time)
    
    def find_daughter(self, parent_name: str) -> Optional[str]:
        """
        Find daughter nuclide (simplified - would need full decay chain data)
        This is a placeholder - in real implementation would use database
        """
        # Common decay chains
        chains = {
            'Sr90': 'Y90',
            'Y90': None,  # stable
            'Cs137': 'Ba137m',
            'Ba137m': 'Ba137',
        }
        
        name_clean = parent_name.replace('-', '').replace('_', '')
        return chains.get(name_clean)
    
    def simulate_chain(self, parent_nuclide: Nuclide,
                      initial_activity: float,
                      time_points: np.ndarray,
                      include_daughters: bool = True) -> Dict[str, np.ndarray]:
        """
        Simulate decay chain evolution
        
        Args:
            parent_nuclide: Parent nuclide
            initial_activity: Initial activity of parent (Bq)
            time_points: Array of time points in seconds
            include_daughters: Whether to include daughter products
        
        Returns:
            Dictionary mapping nuclide names to activity arrays
        """
        results = {}
        
        # Parent decay
        lam_parent = self.get_decay_constant(parent_nuclide)
        parent_activity = initial_activity * np.exp(-lam_parent * time_points)
        results[parent_nuclide.name] = parent_activity
        
        # Daughter decay (simplified - assumes single daughter)
        if include_daughters:
            daughter_name = self.find_daughter(parent_nuclide.name)
            if daughter_name:
                daughter_data = self.database.get_nuclide(daughter_name)
                if daughter_data:
                    daughter = Nuclide(daughter_name, daughter_data)
                    lam_daughter = self.get_decay_constant(daughter)
                    
                    # Bateman equation for parent->daughter
                    if lam_daughter > 0:
                        daughter_activity = (initial_activity * lam_parent / 
                                           (lam_daughter - lam_parent) *
                                           (np.exp(-lam_parent * time_points) - 
                                            np.exp(-lam_daughter * time_points)))
                        results[daughter_name] = daughter_activity
        
        return results
    
    def generate_time_spectrum(self, nuclide: Nuclide,
                              energy_bins: np.ndarray,
                              time: float,
                              initial_activity: float = 1.0,
                              include_daughters: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate spectrum at specific time point
        
        Args:
            nuclide: Nuclide object
            energy_bins: Energy bin edges
            time: Time in seconds
            initial_activity: Initial activity
            include_daughters: Include daughter products
        
        Returns:
            (energies, counts) at time t
        """
        spec_gen = SpectrumGenerator(nuclide)
        bin_centers = (energy_bins[:-1] + energy_bins[1:]) / 2.0
        
        # Get activity at time t
        activity = self.simulate_single_decay(nuclide, initial_activity, time)
        
        # Generate spectrum
        _, counts = spec_gen.generate_total_spectrum(energy_bins)
        
        # Scale by activity
        counts = counts * activity
        
        # Add daughter contributions
        if include_daughters:
            daughter_name = self.find_daughter(nuclide.name)
            if daughter_name:
                daughter_data = self.database.get_nuclide(daughter_name)
                if daughter_data:
                    daughter = Nuclide(daughter_name, daughter_data)
                    daughter_activity = self.simulate_chain(
                        nuclide, initial_activity, np.array([time]), 
                        include_daughters=True
                    ).get(daughter_name, np.array([0.0]))[0]
                    
                    daughter_spec_gen = SpectrumGenerator(daughter)
                    _, daughter_counts = daughter_spec_gen.generate_total_spectrum(energy_bins)
                    counts = counts + daughter_counts * daughter_activity
        
        return bin_centers, counts

