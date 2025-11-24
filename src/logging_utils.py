"""
Logging utilities for the promoter bashing analysis pipeline.
"""

import logging
import sys
from pathlib import Path
from datetime import datetime

def setup_logger(name: str, log_file: str = None, level: int = logging.INFO):
    """
    Set up a logger with consistent formatting.
    
    Parameters:
    -----------
    name : str
        Name of the logger (usually __name__)
    log_file : str, optional
        Path to log file. If None, only logs to console
    level : int
        Logging level (default: INFO)
    
    Returns:
    --------
    logging.Logger
    """
    
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (optional)
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def log_function_call(func):
    """
    Decorator to log function calls with parameters and execution time.
    """
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(func.__module__)
        logger.info(f"Calling {func.__name__} with args={args}, kwargs={kwargs}")
        
        start_time = datetime.now()
        try:
            result = func(*args, **kwargs)
            duration = datetime.now() - start_time
            logger.info(f"{func.__name__} completed successfully in {duration}")
            return result
        except Exception as e:
            duration = datetime.now() - start_time
            logger.error(f"{func.__name__} failed after {duration}: {e}")
            raise
    
    return wrapper
