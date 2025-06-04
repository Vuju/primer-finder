"""Configuration loader for primer-finder."""

import os
import yaml
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Union, List

from primer_finder.matching.dtos.primer_data_dto import primer_info_from_string, PrimerDataDTO

# Type alias for configuration
Config = Dict[str, Any]

class ConfigurationError(Exception):
    """Exception raised for configuration errors."""
    pass

class ConfigLoader:
    """
    Configuration loader for primer-finder.
    
    This class loads configuration from a YAML file, overrides values with
    environment variables, and validates the configuration.
    """
    
    # Default configuration file path
    DEFAULT_CONFIG_PATH = Path(__file__).parent / "default_config.yaml"
    
    # Environment variable prefix for overriding configuration
    ENV_PREFIX = "PRIMER_FINDER_"
    
    def __init__(self, config_path: Optional[Union[str, Path]] = None):
        """
        Initialize the configuration loader.
        
        Args:
            config_path: Path to the configuration file. If None, the default
                configuration file is used.
        """
        self.config_path = Path(config_path) if config_path else self.DEFAULT_CONFIG_PATH
        self.config = self._load_config()
        self._override_from_env()
        self._validate_config()
    
    def _load_config(self) -> Config:
        """
        Load configuration from a YAML file.
        
        Returns:
            The loaded configuration.
            
        Raises:
            ConfigurationError: If the configuration file cannot be loaded.
        """
        try:
            with open(self.config_path, "r") as f:
                return yaml.safe_load(f)
        except Exception as e:
            raise ConfigurationError(f"Failed to load configuration from {self.config_path}: {e}")
    
    def _override_from_env(self) -> None:
        """
        Override configuration values with environment variables.
        
        Environment variables should be prefixed with PRIMER_FINDER_ and use
        double underscores to separate nested keys. For example, to override
        paths.input_file, use PRIMER_FINDER_PATHS__INPUT_FILE.
        """
        for env_var, value in os.environ.items():
            if not env_var.startswith(self.ENV_PREFIX):
                continue
            
            # Remove prefix and split by double underscore to get nested keys
            key_path = env_var[len(self.ENV_PREFIX):].lower().split("__")
            
            # Convert value to appropriate type
            if value.lower() == "true":
                value = True
            elif value.lower() == "false":
                value = False
            elif value.lower() == "null" or value.lower() == "none":
                value = None
            else:
                try:
                    # Try to convert to int or float
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    # Keep as string if conversion fails
                    pass
            
            # Update config with the environment variable value
            self._set_nested_value(self.config, key_path, value)
    
    def _set_nested_value(self, config: Dict[str, Any], key_path: List[str], value: Any) -> None:
        """
        Set a nested value in the configuration.
        
        Args:
            config: The configuration dictionary.
            key_path: List of keys to navigate to the target value.
            value: The value to set.
            
        Raises:
            ConfigurationError: If the key path is invalid.
        """
        if not key_path:
            return
        
        if len(key_path) == 1:
            config[key_path[0]] = value
            return
        
        if key_path[0] not in config:
            config[key_path[0]] = {}
        
        if not isinstance(config[key_path[0]], dict):
            raise ConfigurationError(f"Cannot set nested value for non-dict key: {key_path[0]}")
        
        self._set_nested_value(config[key_path[0]], key_path[1:], value)
    
    def _validate_config(self) -> None:
        """
        Validate the configuration.
        
        Raises:
            ConfigurationError: If the configuration is invalid.
        """
        required_sections = ["paths", "database", "logging", "features", "algorithm", "parallelization"]
        for section in required_sections:
            if section not in self.config:
                raise ConfigurationError(f"Missing required configuration section: {section}")
        
        # Validate paths
        required_paths = ["primer_information", "input_file", "output_file", "log_file"]
        for path_key in required_paths:
            if path_key not in self.config["paths"]:
                raise ConfigurationError(f"Missing required path: {path_key}")

        # Validate database
        required_paths = ["input_table_name", "id_column_name", "sequence_column_name"]
        for path_key in required_paths:
            if path_key not in self.config["database"]:
                raise ConfigurationError(f"Missing required database path: {path_key}")
        
        # Validate logging
        if "level" not in self.config["logging"]:
            raise ConfigurationError("Missing required logging level")
        
        # Validate features
        required_features = ["enable_primer_finder", "enable_orf_finder"]
        for feature in required_features:
            if feature not in self.config["features"]:
                raise ConfigurationError(f"Missing required feature toggle: {feature}")
        
        # Validate algorithm parameters
        required_algorithm_params = [
            "search_area", "smith_waterman_score_cutoff", "gap_penalty",
            "triplet_gap_penalty", "end_of_read_bonus", "orf_matching_lower_threshold",
            "orf_matching_upper_threshold", "protein_translation_table", "e_value"
        ]
        for param in required_algorithm_params:
            if param not in self.config["algorithm"]:
                raise ConfigurationError(f"Missing required algorithm parameter: {param}")
        
        # Validate parallelization settings
        required_parallelization_params = ["chunk_size"]
        for param in required_parallelization_params:
            if param not in self.config["parallelization"]:
                raise ConfigurationError(f"Missing required parallelization parameter: {param}")
    
    def get_config(self) -> Config:
        """
        Get the complete configuration.
        
        Returns:
            The complete configuration.
        """
        return self.config
    
    def get(self, section: str, key: str, default: Any = None) -> Any:
        """
        Get a configuration value.
        
        Args:
            section: The configuration section.
            key: The configuration key.
            default: The default value to return if the key is not found.
            
        Returns:
            The configuration value, or the default value if not found.
        """
        if section not in self.config:
            return default
        
        return self.config[section].get(key, default)


# Create a singleton instance of the configuration loader
_config_loader = None

def get_config_loader(config_path: Optional[Union[str, Path]] = None) -> ConfigLoader:
    """
    Get the configuration loader instance.
    
    Args:
        config_path: Path to the configuration file. If None, the default
            configuration file is used.
            
    Returns:
        The configuration loader instance.
    """
    global _config_loader
    if _config_loader is None or config_path is not None:
        _config_loader = ConfigLoader(config_path)
    return _config_loader


def get_primer_information(primer_information_file, primer_data):
    with open(primer_information_file, "r") as primer_info_file:
        for line in primer_info_file.readlines():
            entry = primer_info_from_string(line)
            primer_data.append(entry)