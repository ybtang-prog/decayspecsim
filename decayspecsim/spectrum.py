"""
Spectrum generation from nuclide decay data
"""

import numpy as np
from typing import Tuple, Optional
from scipy import special
from .nuclide import Nuclide


class SpectrumGenerator:
    """Generate theoretical spectra from nuclide decay data"""
    
    def __init__(self, nuclide: Nuclide):
        self.nuclide = nuclide
    
    def generate_discrete_spectrum(self, decay_type: str, 
                                   energy_bins: np.ndarray,
                                   normalize: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate discrete spectrum (gamma, alpha, x-ray, electron)
        
        Args:
            decay_type: 'gamma', 'alpha', 'x-ray', or 'electron'
            energy_bins: Energy bin edges in keV
            normalize: If True, normalize to total intensity
        
        Returns:
            (energies, counts) arrays
        """
        decay_mode = self.nuclide.get_decay_mode(decay_type)
        if decay_mode is None:
            return np.array([]), np.array([])
        
        energies = decay_mode.energies / 1000.0  # Convert to keV
        intensities = decay_mode.get_effective_intensities()
        
        # Bin the spectrum
        counts = np.zeros(len(energy_bins) - 1)
        bin_centers = (energy_bins[:-1] + energy_bins[1:]) / 2.0
        
        for e, intensity in zip(energies, intensities):
            if intensity > 0 and e > 0:
                idx = np.searchsorted(energy_bins, e) - 1
                if 0 <= idx < len(counts):
                    counts[idx] += intensity
        
        if normalize and np.sum(counts) > 0:
            counts = counts / np.sum(counts) * np.sum(intensities)
        
        return bin_centers, counts
    
    def generate_beta_spectrum(self, energy_bins: np.ndarray,
                              normalize: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate continuous beta spectrum using Fermi theory
        
        Args:
            energy_bins: Energy bin edges in keV
            normalize: If True, normalize to total intensity
        
        Returns:
            (energies, counts) arrays
        """
        decay_mode = self.nuclide.get_decay_mode('beta')
        if decay_mode is None:
            return np.array([]), np.array([])
        
        bin_centers = (energy_bins[:-1] + energy_bins[1:]) / 2.0
        bin_width = energy_bins[1] - energy_bins[0]
        counts = np.zeros(len(bin_centers))
        
        # Process each beta branch
        endpoint_energies = decay_mode.energies / 1000.0  # keV
        intensities = decay_mode.get_effective_intensities()
        
        for Q, intensity in zip(endpoint_energies, intensities):
            if Q <= 0 or intensity <= 0:
                continue
            
            # Calculate beta spectrum for this branch
            # Use vectorized calculation for better performance
            valid_mask = (bin_centers > 0) & (bin_centers < Q)
            if not np.any(valid_mask):
                continue
            
            E_valid = bin_centers[valid_mask]
            
            # Fermi theory: N(E) ~ F(Z,E) * p * E_total * (Q - E)^2
            # where p is momentum, E_total is total energy
            m_e = 0.511  # electron rest mass in MeV
            E_MeV = E_valid / 1000.0  # Convert to MeV
            Q_MeV = Q / 1000.0
            
            # Momentum in MeV/c
            p = np.sqrt(E_MeV * (E_MeV + 2 * m_e))
            E_total = E_MeV + m_e
            
            # Simplified Fermi function approximation
            # F(Z,E) = 2 * pi * eta / (1 - exp(-2*pi*eta))
            # where eta = alpha * Z * E_total / p
            Z = self.nuclide.zai // 10000 if self.nuclide.zai else 1
            alpha = 1.0 / 137.036  # fine structure constant
            eta = alpha * Z * E_total / np.maximum(p, 1e-10)
            
            # Avoid division by zero
            fermi_func = np.where(np.abs(eta) < 1e-10, 
                                 1.0,
                                 2 * np.pi * eta / (1 - np.exp(-2 * np.pi * eta)))
            
            # Beta spectrum shape
            spectrum_values = fermi_func * p * E_total * (Q_MeV - E_MeV) ** 2
            counts[valid_mask] += intensity * spectrum_values * bin_width
        
        # Normalize
        if normalize and np.sum(counts) > 0:
            # Normalize so integral matches total intensity
            total_intensity = np.sum(intensities)
            integral = np.sum(counts) * bin_width
            if integral > 0:
                counts = counts / integral * total_intensity
        
        return bin_centers, counts
    
    def generate_total_spectrum(self, energy_bins: np.ndarray,
                               decay_types: Optional[list] = None,
                               include_beta: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate total spectrum combining all decay modes
        
        Args:
            energy_bins: Energy bin edges in keV
            decay_types: List of discrete decay types to include. If None, includes all
            include_beta: Whether to include beta continuous spectrum
        
        Returns:
            (energies, counts) arrays
        """
        bin_centers = (energy_bins[:-1] + energy_bins[1:]) / 2.0
        total_counts = np.zeros(len(bin_centers))
        
        # Add discrete spectra
        if decay_types is None:
            decay_types = ['gamma', 'alpha', 'x-ray', 'electron']
        
        for decay_type in decay_types:
            if self.nuclide.has_decay_mode(decay_type):
                _, counts = self.generate_discrete_spectrum(decay_type, energy_bins, normalize=False)
                if len(counts) == len(total_counts):
                    total_counts += counts
        
        # Add beta spectrum
        if include_beta and self.nuclide.has_decay_mode('beta'):
            _, beta_counts = self.generate_beta_spectrum(energy_bins, normalize=False)
            if len(beta_counts) == len(total_counts):
                total_counts += beta_counts
        
        return bin_centers, total_counts

