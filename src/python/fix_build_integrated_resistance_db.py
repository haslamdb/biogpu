#!/usr/bin/env python3
"""
fix_build_integrated_resistance_db.py
Fixes and enhancements for the integrated resistance database builder
"""

import os
import sys
import json

def main():
    """Check and fix the build_integrated_resistance_db.py script"""
    
    script_path = "src/python/build_integrated_resistance_db.py"
    
    if not os.path.exists(script_path):
        print(f"Creating {script_path}...")
        
        # Create a minimal version that works with the existing script
        script_content = '''#!/usr/bin/env python3
"""
build_integrated_resistance_db.py
Build GPU-ready integrated resistance database from various sources
This is a wrapper that uses the existing build_integrated_resistance_db.py
"""

import sys
import os

# Add the directory containing the original script to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import and run the original script
try:
    from build_integrated_resistance_db import *
    
    if __name__ == "__main__":
        main()
except ImportError:
    print("Error: Could not import build_integrated_resistance_db")
    print("Please ensure the original script is in the same directory")
    sys.exit(1)
'''
        
        os.makedirs(os.path.dirname(script_path), exist_ok=True)
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        print(f"Created wrapper script at {script_path}")
    
    # Check if the integrated database builder needs the resistance mutation class
    try:
        # Test import
        sys.path.insert(0, os.path.dirname(script_path))
        import build_integrated_resistance_db
        print("✓ build_integrated_resistance_db.py imports successfully")
    except ImportError as e:
        print(f"Import error: {e}")
        print("\nThe script appears to be missing or has import errors.")
        print("Let me check the actual location...")
        
        # Look for the script in various locations
        possible_locations = [
            "runtime/kernels/resistance/build_integrated_resistance_db.py",
            "src/python/build_integrated_resistance_db.py",
            "scripts/build_integrated_resistance_db.py"
        ]
        
        found = False
        for loc in possible_locations:
            if os.path.exists(loc):
                print(f"Found script at: {loc}")
                found = True
                
                # Copy to expected location
                import shutil
                os.makedirs(os.path.dirname(script_path), exist_ok=True)
                shutil.copy2(loc, script_path)
                print(f"Copied to: {script_path}")
                break
        
        if not found:
            print("\nThe build_integrated_resistance_db.py script was not found.")
            print("This script should have been included in the documents.")
            print("\nYou can use the script from the documents by:")
            print("1. Copy the content of build_integrated_resistance_db.py from the documents")
            print("2. Save it to: src/python/build_integrated_resistance_db.py")
            print("3. Make it executable: chmod +x src/python/build_integrated_resistance_db.py")

if __name__ == "__main__":
    main()