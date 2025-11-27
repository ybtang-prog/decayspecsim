# DecaySpecSim

Nuclear decay spectrum simulator based on Decay2012 database.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Print nuclide information

```bash
python decayspecsim.py --nuclide Co-60 --info
```

### Generate theoretical spectrum

```bash
python decayspecsim.py --nuclide Co-60 --type gamma --output co60_theoretical.csv
```

### Simulate detector response

```bash
# HPGe detector
python decayspecsim.py --nuclide Co-60 --type gamma --detector hpge --plot co60_hpge.png

# NaI detector
python decayspecsim.py --nuclide Co-60 --type gamma --detector nai --plot co60_nai.png
```

### Beta spectrum

```bash
python decayspecsim.py --nuclide H-3 --type beta --plot h3_beta.png
```

### Decay chain evolution

```bash
python decayspecsim.py --nuclide Sr-90 --decay_chain --time_points "0, 7d, 30d, 365d" --type beta --plot sr90_evolution.png
```

## Command-line Options

- `--nuclide, -n`: Nuclide name (required)
- `--info`: Print nuclide information
- `--type, -t`: Spectrum type (all, gamma, beta, alpha, x-ray, electron)
- `--output, -o`: Output CSV file
- `--plot, -p`: Output plot file
- `--detector, -d`: Detector type (hpge, nai, generic)
- `--resolution`: Resolution parameters "a,b,c" for FWHM=sqrt(a+b*E+c*E^2)
- `--energy_range`: Energy range "min,max" in keV
- `--bins`: Number of energy bins
- `--decay_chain`: Enable decay chain simulation
- `--time_points`: Time points for evolution (e.g., "0, 7d, 30d")
- `--activity`: Initial activity in Bq

## Examples

See the test cases in the documentation for detailed examples.

