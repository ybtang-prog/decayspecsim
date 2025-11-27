"""
Database module for loading and querying Decay2012 data
"""

import json
import os
from typing import Dict, Optional, Union


class Database:
    """Loads and provides access to Decay2012 database"""
    
    def __init__(self, db_path: Optional[str] = None):
        """
        Initialize database from JSON file
        
        Args:
            db_path: Path to Decay2012Database.json. If None, looks for 
                    lines_decay_2012.min.json in current directory
        """
        if db_path is None:
            # Try common names
            for name in ['lines_decay_2012.min.json', 'Decay2012Database.json']:
                if os.path.exists(name):
                    db_path = name
                    break
            if db_path is None:
                raise FileNotFoundError("Database file not found")
        
        with open(db_path, 'r', encoding='utf-8') as f:
            self._data = json.load(f)
    
    def get_nuclide(self, identifier: Union[str, int]) -> Optional[Dict]:
        """
        Get nuclide data by name or ZAI
        
        Args:
            identifier: Nuclide name (e.g., 'Co-60', 'Co60', 'H-3') or ZAI value
        
        Returns:
            Dictionary with nuclide data or None if not found
        """
        if isinstance(identifier, int):
            # Search by ZAI
            for name, data in self._data.items():
                if data.get('zai') == identifier:
                    return data
            return None
        
        # Search by name - try different formats
        name = identifier.replace('-', '').replace('_', '')
        
        # Try exact match first
        if name in self._data:
            return self._data[name]
        
        # Try case-insensitive
        name_lower = name.lower()
        for key in self._data.keys():
            if key.lower() == name_lower:
                return self._data[key]
        
        # Try with 'm' suffix for metastable states
        if name + 'm' in self._data:
            return self._data[name + 'm']
        
        return None
    
    def list_nuclides(self) -> list:
        """Return list of all nuclide names in database"""
        return list(self._data.keys())
    
    def has_nuclide(self, identifier: Union[str, int]) -> bool:
        """Check if nuclide exists in database"""
        return self.get_nuclide(identifier) is not None

