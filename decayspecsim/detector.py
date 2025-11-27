"""
Detector response simulation
"""

import numpy as np
from typing import Tuple, Optional, Callable
from scipy.ndimage import gaussian_filter1d


class Detector:
    """Simulates detector response to radiation"""
    
    def __init__(self, name: str = "generic",
                 resolution_params: Optional[Tuple[float, float, float]] = None,
                 efficiency_curve: Optional[np.ndarray] = None,
                 efficiency_energies: Optional[np.ndarray] = None):
        """
        Initialize detector
        
        Args:
            name: Detector name (e.g., 'hpge', 'nai')
            resolution_params: (a, b, c) for FWHM = sqrt(a + b*E + c*E^2) in keV
            efficiency_curve: Efficiency values
            efficiency_energies: Energy values for efficiency curve (keV)
        """
        self.name = name
        
        # Default resolution parameters
        if resolution_params is None:
            if name.lower() == 'hpge':
                # HPGe: good resolution
                resolution_params = (0.1, 0.0005, 0.0)
            elif name.lower() == 'nai':
                # NaI: poor resolution
                resolution_params = (0.0, 0.006, 0.0)
            else:
                resolution_params = (0.0, 0.003, 0.0)
        
        self.resolution_a, self.resolution_b, self.resolution_c = resolution_params
        
        # Efficiency
        self.efficiency_curve = efficiency_curve
        self.efficiency_energies = efficiency_energies
    
    def get_fwhm(self, energy: np.ndarray) -> np.ndarray:
        """Calculate FWHM at given energies"""
        return np.sqrt(self.resolution_a + self.resolution_b * energy + 
                      self.resolution_c * energy ** 2)
    
    def get_sigma(self, energy: np.ndarray) -> np.ndarray:
        """Convert FWHM to sigma for Gaussian"""
        return self.get_fwhm(energy) / 2.355
    
    def get_efficiency(self, energy: np.ndarray) -> np.ndarray:
        """
        Get efficiency at given energies
        
        Args:
            energy: Energy values in keV
        
        Returns:
            Efficiency values (1.0 if no curve provided)
        """
        if self.efficiency_curve is None:
            return np.ones_like(energy)
        
        # Interpolate efficiency curve
        return np.interp(energy, self.efficiency_energies, self.efficiency_curve)
    
    def apply_resolution(self, energies: np.ndarray, counts: np.ndarray) -> np.ndarray:
        """
        Apply energy resolution (Gaussian smearing)
        
        Uses convolution with variable-width Gaussian kernel
        
        Args:
            energies: Energy bin centers
            counts: Counts in each bin
        
        Returns:
            Smeared spectrum
        """
        if len(energies) == 0:
            return counts
        
        # Use simplified approach: average sigma for convolution
        # More accurate would be variable-width convolution, but slower
        sigmas = self.get_sigma(energies)
        avg_sigma = np.mean(sigmas[sigmas > 0])
        
        if avg_sigma <= 0:
            return counts
        
        bin_width = energies[1] - energies[0] if len(energies) > 1 else 1.0
        
        # Create Gaussian kernel
        kernel_size = int(6 * avg_sigma / bin_width)
        if kernel_size < 3:
            kernel_size = 3
        if kernel_size > len(energies):
            kernel_size = len(energies)
        
        kernel = np.arange(-kernel_size, kernel_size + 1) * bin_width
        kernel = np.exp(-0.5 * (kernel / avg_sigma) ** 2)
        kernel = kernel / np.sum(kernel)
        
        # Convolve
        result = np.convolve(counts, kernel, mode='same')
        
        return result
    
    def apply_efficiency(self, energies: np.ndarray, counts: np.ndarray) -> np.ndarray:
        """Apply detection efficiency"""
        efficiency = self.get_efficiency(energies)
        return counts * efficiency
    
    def add_compton_continuum(self, energies: np.ndarray, counts: np.ndarray,
                             gamma_energies: np.ndarray,
                             gamma_intensities: np.ndarray) -> np.ndarray:
        """
        Add simplified Compton continuum for gamma rays
        
        Args:
            energies: Energy bin centers
            counts: Current counts
            gamma_energies: Gamma ray energies in keV
            gamma_intensities: Gamma ray intensities
        
        Returns:
            Spectrum with Compton continuum added
        """
        result = counts.copy()
        bin_width = energies[1] - energies[0] if len(energies) > 1 else 1.0
        
        for E_gamma, intensity in zip(gamma_energies, gamma_intensities):
            if E_gamma <= 0 or intensity <= 0:
                continue
            
            # Compton edge energy
            m_e = 511.0  # electron rest mass in keV
            E_compton = E_gamma / (1 + m_e / (2 * E_gamma))
            
            # Add flat continuum from 0 to Compton edge
            # Area proportional to photopeak area
            compton_area = intensity * 0.3  # rough estimate
            compton_counts_per_bin = compton_area / (E_compton / bin_width) if E_compton > 0 else 0
            
            for i, E in enumerate(energies):
                if 0 < E < E_compton:
                    result[i] += compton_counts_per_bin
        
        return result
    
    def add_escape_peaks(self, energies: np.ndarray, counts: np.ndarray,
                        gamma_energies: np.ndarray,
                        gamma_intensities: np.ndarray) -> np.ndarray:
        """
        Add single and double escape peaks for gamma rays (HPGe only)
        
        Args:
            energies: Energy bin centers
            counts: Current counts
            gamma_energies: Gamma ray energies in keV
            gamma_intensities: Gamma ray intensities
        
        Returns:
            Spectrum with escape peaks added
        """
        if self.name.lower() != 'hpge':
            return counts
        
        result = counts.copy()
        bin_width = energies[1] - energies[0] if len(energies) > 1 else 1.0
        
        for E_gamma, intensity in zip(gamma_energies, gamma_intensities):
            if E_gamma <= 511.0:
                continue
            
            # Single escape peak
            E_single = E_gamma - 511.0
            # Double escape peak
            E_double = E_gamma - 1022.0
            
            # Escape peak intensities (rough estimates)
            single_intensity = intensity * 0.01
            double_intensity = intensity * 0.001 if E_double > 0 else 0
            
            # Add as narrow peaks
            for E_peak, peak_intensity in [(E_single, single_intensity), 
                                          (E_double, double_intensity)]:
                if E_peak > 0:
                    idx = np.argmin(np.abs(energies - E_peak))
                    if 0 <= idx < len(result):
                        sigma = self.get_sigma(E_peak)
                        # Add Gaussian peak
                        for i in range(max(0, idx - 5), min(len(energies), idx + 6)):
                            dE = energies[i] - E_peak
                            gaussian = np.exp(-0.5 * (dE / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
                            result[i] += peak_intensity * gaussian * bin_width
        
        return result
    
    def simulate(self, energies: np.ndarray, counts: np.ndarray,
                gamma_energies: Optional[np.ndarray] = None,
                gamma_intensities: Optional[np.ndarray] = None,
                include_compton: bool = True,
                include_escape: bool = True) -> np.ndarray:
        """
        Full detector simulation
        
        Args:
            energies: Energy bin centers
            counts: Theoretical counts
            gamma_energies: Gamma ray energies (for Compton/escape peaks)
            gamma_intensities: Gamma ray intensities
            include_compton: Add Compton continuum
            include_escape: Add escape peaks (HPGe only)
        
        Returns:
            Simulated detector spectrum
        """
        result = counts.copy()
        
        # Apply efficiency
        result = self.apply_efficiency(energies, result)
        
        # Add Compton continuum
        if include_compton and gamma_energies is not None:
            result = self.add_compton_continuum(energies, result, 
                                               gamma_energies, gamma_intensities)
        
        # Add escape peaks
        if include_escape and gamma_energies is not None:
            result = self.add_escape_peaks(energies, result,
                                         gamma_energies, gamma_intensities)
        
        # Apply resolution (do this last)
        result = self.apply_resolution(energies, result)
        
        return result

