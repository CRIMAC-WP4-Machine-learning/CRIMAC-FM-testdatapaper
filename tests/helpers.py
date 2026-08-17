import subprocess


def run(command):
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"Command failed: {' '.join(command)}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    return result
