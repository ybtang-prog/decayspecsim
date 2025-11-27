"""
Command-line interface for DecaySpecSim
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional

from .database import Database
from .nuclide import Nuclide
from .spectrum import SpectrumGenerator
from .detector import Detector
from .decay_chain import DecayChainSimulator


def parse_resolution(s: str) -> tuple:
    """Parse resolution parameters from string"""
    try:
        parts = [float(x.strip()) for x in s.split(',')]
        if len(parts) == 3:
            return tuple(parts)
        elif len(parts) == 1:
            return (0.0, parts[0], 0.0)
        else:
            raise ValueError("Resolution must have 1 or 3 values")
    except:
        raise ValueError(f"Invalid resolution format: {s}")


def parse_time_points(s: str) -> np.ndarray:
    """Parse time points from string (e.g., '0, 7d, 30d')"""
    import re
    times = []
    for part in s.split(','):
        part = part.strip()
        # Parse with units
        match = re.match(r'([\d.]+)\s*([dhms]?)', part.lower())
        if match:
            value = float(match.group(1))
            unit = match.group(2)
            # Convert to seconds
            multipliers = {'d': 86400, 'h': 3600, 'm': 60, 's': 1, '': 1}
            times.append(value * multipliers.get(unit, 1))
        else:
            times.append(float(part))
    return np.array(times)


def print_nuclide_info(nuclide: Nuclide):
    """Print information about nuclide"""
    print(f"\nNuclide: {nuclide.name} (ZAI: {nuclide.zai})")
    print(f"Half-life: {nuclide.halflife:.6e} s ({nuclide.halflife / 86400:.2f} days)")
    print(f"\nDecay Modes:")
    
    for mode_type in nuclide.list_decay_modes():
        mode = nuclide.get_decay_mode(mode_type)
        print(f"  - {mode_type.capitalize()}: {mode.number} lines")
        if mode.number > 0:
            energies = mode.energies / 1000.0  # keV
            intensities = mode.intensities
            norms = mode.norms
            for i in range(min(5, len(energies))):  # Show first 5
                eff_int = intensities[i] * norms[i]
                if eff_int > 0.01:
                    print(f"    - E={energies[i]:.2f} keV, I={eff_int:.2f}%")
                else:
                    print(f"    - E={energies[i]:.2f} keV, I={eff_int:.6f}%")
            if len(energies) > 5:
                print(f"    ... and {len(energies) - 5} more")


def save_spectrum(energies: np.ndarray, counts: np.ndarray, filename: str):
    """Save spectrum to CSV file"""
    data = np.column_stack([energies, counts])
    np.savetxt(filename, data, delimiter=',', header='Energy (keV),Counts', 
               comments='', fmt='%.6e')


def plot_spectrum(energies: np.ndarray, counts: np.ndarray, 
                 filename: Optional[str] = None,
                 title: str = "Spectrum",
                 xlabel: str = "Energy (keV)",
                 ylabel: str = "Counts"):
    """Plot and optionally save spectrum"""
    plt.figure(figsize=(10, 6))
    plt.plot(energies, counts, 'b-', linewidth=1.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {filename}")
    else:
        plt.show()
    
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='DecaySpecSim - Nuclear decay spectrum simulator',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--nuclide', '-n', type=str, required=True,
                       help='Nuclide name (e.g., Co-60, Cs-137, H-3)')
    parser.add_argument('--info', action='store_true',
                       help='Print nuclide information')
    parser.add_argument('--type', '-t', type=str, default='all',
                       choices=['all', 'gamma', 'beta', 'alpha', 'x-ray', 'electron'],
                       help='Type of spectrum to generate')
    parser.add_argument('--output', '-o', type=str,
                       help='Output CSV file for spectrum data')
    parser.add_argument('--plot', '-p', type=str,
                       help='Output plot file (PNG/PDF)')
    parser.add_argument('--detector', '-d', type=str,
                       choices=['hpge', 'nai', 'generic'],
                       help='Detector type for simulation')
    parser.add_argument('--resolution', type=str,
                       help='Resolution parameters: "a,b,c" for FWHM=sqrt(a+b*E+c*E^2)')
    parser.add_argument('--energy_range', type=str, default='0,3000',
                       help='Energy range: "min,max" in keV')
    parser.add_argument('--bins', type=int, default=1000,
                       help='Number of energy bins')
    parser.add_argument('--decay_chain', action='store_true',
                       help='Simulate decay chain')
    parser.add_argument('--time_points', type=str,
                       help='Time points for evolution: "0, 7d, 30d"')
    parser.add_argument('--activity', type=float, default=1.0,
                       help='Initial activity in Bq')
    
    args = parser.parse_args()
    
    # Load database
    try:
        db = Database()
    except FileNotFoundError:
        print("Error: Database file not found. Please ensure lines_decay_2012.min.json is in the current directory.")
        sys.exit(1)
    
    # Get nuclide
    nuclide_data = db.get_nuclide(args.nuclide)
    if nuclide_data is None:
        print(f"Error: Nuclide '{args.nuclide}' not found in database")
        sys.exit(1)
    
    nuclide = Nuclide(args.nuclide, nuclide_data)
    
    # Print info if requested
    if args.info:
        print_nuclide_info(nuclide)
        if not args.output and not args.plot:
            return
    
    # Parse energy range
    try:
        emin, emax = map(float, args.energy_range.split(','))
    except:
        emin, emax = 0.0, 3000.0
    
    energy_bins = np.linspace(emin, emax, args.bins + 1)
    
    # Generate spectrum
    spec_gen = SpectrumGenerator(nuclide)
    
    if args.decay_chain and args.time_points:
        # Time evolution
        chain_sim = DecayChainSimulator(db)
        time_points = parse_time_points(args.time_points)
        
        # Generate spectra for each time point
        all_energies = []
        all_counts = []
        labels = []
        
        for t in time_points:
            energies, counts = chain_sim.generate_time_spectrum(
                nuclide, energy_bins, t, args.activity
            )
            all_energies.append(energies)
            all_counts.append(counts)
            labels.append(f"t = {t/86400:.1f} days" if t > 86400 else f"t = {t:.0f} s")
        
        # Plot multiple curves
        if args.plot:
            plt.figure(figsize=(10, 6))
            for energies, counts, label in zip(all_energies, all_counts, labels):
                plt.plot(energies, counts, label=label, linewidth=1.5)
            plt.xlabel('Energy (keV)')
            plt.ylabel('Counts')
            plt.title(f'{nuclide.name} Decay Chain Evolution')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            plt.savefig(args.plot, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {args.plot}")
            plt.close()
        
        # Save last time point
        if args.output:
            save_spectrum(all_energies[-1], all_counts[-1], args.output)
    else:
        # Single spectrum
        if args.type == 'all':
            energies, counts = spec_gen.generate_total_spectrum(energy_bins)
        elif args.type == 'beta':
            energies, counts = spec_gen.generate_beta_spectrum(energy_bins)
        else:
            energies, counts = spec_gen.generate_discrete_spectrum(args.type, energy_bins)
        
        # Apply detector simulation if requested
        if args.detector:
            detector_params = {}
            if args.resolution:
                detector_params['resolution_params'] = parse_resolution(args.resolution)
            
            detector = Detector(name=args.detector, **detector_params)
            
            # Get gamma data for Compton/escape peaks
            gamma_energies = None
            gamma_intensities = None
            if nuclide.has_decay_mode('gamma'):
                gamma_mode = nuclide.get_decay_mode('gamma')
                gamma_energies = gamma_mode.energies / 1000.0
                gamma_intensities = gamma_mode.get_effective_intensities()
            
            counts = detector.simulate(energies, counts, gamma_energies, gamma_intensities)
        
        # Save output
        if args.output:
            save_spectrum(energies, counts, args.output)
            print(f"Spectrum saved to {args.output}")
        
        # Plot
        if args.plot:
            title = f"{nuclide.name} {args.type.capitalize()} Spectrum"
            if args.detector:
                title += f" ({args.detector.upper()} Detector)"
            plot_spectrum(energies, counts, args.plot, title=title)


if __name__ == '__main__':
    main()

