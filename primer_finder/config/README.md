# Configuration System

This directory contains the configuration system for the primer-finder project. The configuration system provides a way to configure the application using a YAML file and environment variables.

## Overview

The configuration system consists of the following components:

- `default_config.yaml`: The default configuration file.
- `config_loader.py`: A module that loads and validates the configuration.
- `example_config_usage.py`: An example script demonstrating how to use the configuration system.

## Configuration File

The configuration is stored in a YAML file. The default configuration file is `default_config.yaml`. You can create your own configuration file and specify its path when initializing the configuration loader.

The configuration file is organized into the following sections:

- `paths`: File paths for various resources.
- `logging`: Logger settings.
- `features`: Feature toggles.
- `algorithm`: Algorithm parameters for primer finding and ORF matching.
- `parallelization`: Settings for parallel processing.

## Using the Configuration System

### Basic Usage

```python
from primer_finder.config import get_config_loader

# Get the configuration loader
config_loader = get_config_loader()

# Get the complete configuration
config = config_loader.get_config()

# Get a specific configuration value
input_file = config_loader.get("paths", "input_file")
```

### Using a Custom Configuration File

```python
from primer_finder.config import get_config_loader

# Get the configuration loader with a custom configuration file
config_loader = get_config_loader("path/to/custom_config.yaml")

# Get the complete configuration
config = config_loader.get_config()
```

### Using Environment Variables

You can override configuration values using environment variables. Environment variables should be prefixed with `PRIMER_FINDER_` and use double underscores to separate nested keys. For example, to override `paths.input_file`, use `PRIMER_FINDER_PATHS__INPUT_FILE`.

```bash
# Set environment variables
export PRIMER_FINDER_PATHS__INPUT_FILE="./data/custom_input.fna"
export PRIMER_FINDER_ALGORITHM__E_VALUE=1000
```

## Configuration Validation

The configuration is validated when loaded. If a required value is missing, a `ConfigurationError` will be raised.

## Example

See `example_config_usage.py` for a complete example of how to use the configuration system.
