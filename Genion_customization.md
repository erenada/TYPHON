# Using a Custom Genion Build with Typhon

## Context
Typhon requires a version of Genion that outputs read IDs for all potential chimeric RNAs, including those that do not pass Genion's default filtering steps. This is achieved by running Genion in a 'debug' mode, which required modifying the Genion source code (see `annotate.cpp` and `annotate_readid_report.txt` for details).

## Why Not Upstream?
An upstream contribution to the official Genion repository may not be suitable in this case, as the modification fundamentally changes Genion's output to report all potential chimeric RNAs, not just those passing its filters. This is more of a debug or research mode rather than a general feature.

## Options for Integrating Custom Genion with Typhon

### 1. Maintain a Forked Version
- Fork the Genion repository and apply your modifications.
- Build and host your own conda package for your forked Genion.
- Update Typhon's `environment.yml` to install your custom Genion from your conda channel.
- **Recommended for reproducibility and ease of use.**

### 2. Build Genion from Source as Part of Typhon Setup
- Include a script or module in Typhon that clones your fork, applies the patch, and builds Genion from source.
- Less user-friendly, but ensures the correct version is used.

### 3. Provide a Pre-built Binary
- Build your custom Genion binary for common platforms.
- Host it (e.g., in your GitHub releases).
- Instruct users to download and use this binary, or automate this in Typhon's setup.

## Recommended Approach
- **Fork and build a custom conda package** for your modified Genion.
- Host it on Anaconda Cloud or a private channel.
- Update Typhon's documentation and `environment.yml` to use your custom package.

## Documentation and Reproducibility
- Keep `annotate_readid_report.txt` up to date with a summary of your changes.
- Add a section in Typhon's documentation about the custom Genion, why it's needed, and how it's built/installed.
- Add a check in Typhon to verify the correct Genion version is installed.

## Summary Table
| Approach                | User Effort | Reproducibility | Recommended? |
|-------------------------|-------------|-----------------|--------------|
| Upstream contribution   | Low         | High            | Not suitable |
| Fork + conda package    | Low         | High            | Good         |
| Build from source       | Medium      | Medium          | Acceptable   |
| Pre-built binary        | Low         | Medium          | Acceptable   |

If you need help with any of these steps (e.g., writing a conda recipe, setup script, or documentation), just ask! 