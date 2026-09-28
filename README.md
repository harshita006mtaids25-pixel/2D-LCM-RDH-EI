
# 2D Lozi–Chebyshev Chaotic Map Based RDH-EI

MATLAB implementation of a **Reversible Data Hiding in Encrypted Images (RDH-EI)** scheme using the proposed **2D Lozi–Chebyshev Chaotic Map (2D-LCM)** and **Daubechies 5/3 Integer Wavelet Transform (IWT)**.

## Overview

The proposed framework combines:

* 2D Lozi–Chebyshev chaotic map for encryption
* Daubechies 5/3 Integer Wavelet Transform
* Histogram-shifting based reversible data embedding
* Chaotic permutation and XOR-based encryption
* Plaintext-dependent diffusion
* Exact secret-data extraction
* Bit-exact cover-image recovery

The complete implementation also performs chaotic-map characterization and security/performance analysis.

## Main Pipeline

```text
Cover Image
     ↓
FOB Preprocessing
     ↓
Daubechies 5/3 IWT
     ↓
Histogram-Shifting Embedding
     ↓
Inverse IWT
     ↓
Marked Image
     ↓
2D-LCM Chaotic Encryption
     ↓
Encrypted Stego Image
```

The receiver performs the inverse operations to recover both the secret data and the original cover image.

## Chaotic Map

The proposed map is:

$$
x_{i+1}=1-a|\cos(w\cos^{-1}(x_i))|+y_i
$$

$$
y_{i+1}=bx_i
$$

The implementation includes Lyapunov-exponent verification, bifurcation analysis, phase-space analysis, ergodicity analysis, 0–1 chaos testing, and comparison with existing chaotic maps.

## Security Analysis

The implementation evaluates:

* Shannon entropy
* Horizontal/vertical/diagonal correlation
* NPCR
* UACI
* Key sensitivity
* Chi-square uniformity
* 0–1 chaos test
* Lyapunov exponents

## RDH Performance

The implementation evaluates:

* Embedding capacity (bpp)
* PSNR
* Capacity–PSNR relationship
* Exact secret recovery
* Bit-exact cover-image recovery

## Requirements

* MATLAB R2016b or later
* No additional toolbox required for the core implementation

## Usage

Set the cover and secret-image paths in the configuration section:

```matlab
cfg.imagePath = 'data/cover/Jetplane.tiff';
cfg.secretImagePath = 'data/secret/secret.tiff';
```

Then run:

```matlab
sidflozi_chebyshev_RDH_LCM_only
```

All generated figures and numerical results are saved in the configured `results` directory.

## Repository Structure

```text
├── README.md
├── src/
│   └── sidflozi_chebyshev_RDH_LCM_only.m
├── data/
│   ├── cover/
│   └── secret/
├── results/
└── docs/
```

## Reproducibility

The implementation uses fixed random initialization and reports the chaotic parameters, embedding rate, security metrics, PSNR, capacity, and reversibility results.

For reproducible comparison, always report the **embedding rate (bpp)** together with PSNR.

## Citation

If you use this implementation in your research, please cite the associated research paper.

## License

See `LICENSE` for the terms of use.
