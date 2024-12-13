import subprocess
import os

def test_main_execution():
    # Run the script as a subprocess
    result = subprocess.run(
        ["python3", "main.py"],
        capture_output=True,
        text=True
    )

    # Assert that it executed without errors
    assert result.returncode == 0, f"Error: {result.stderr}"
    print(result.stdout)
