"""Configuration package for primer-finder."""

from primer_finder.config.config_loader import get_config_loader, ConfigLoader, ConfigurationError

__all__ = ["get_config_loader", "ConfigLoader", "ConfigurationError"]

# Initialize the configuration loader
get_config_loader()