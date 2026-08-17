import re
import sys
from pathlib import Path
import os

label = sys.argv[1]

version_file = Path("src/mmaseq/__version__.py")

text = version_file.read_text()

m = re.search(r'__version__ = "(\d+)\.(\d+)\.(\d+)"', text)

major, minor, patch = map(int, m.groups())

print(
    "Before",
    f"label: {label}",
    f"major: {major}",
    f"minor: {minor}",
    f"patch: {patch}",
    sep = "\n"
)

if label == "major":
    major += 1
    minor = 0
    patch = 0

elif label == "minor":
    minor += 1
    # Determine every tens of minor version increments (10, 20, 30 etc.)
    if (minor % 10) == 0:
        # Set patch version to 1 directly
        patch = 1
    patch = 0

elif label == "patch":
    patch += 1


new_version = f'{major}.{minor}.{patch}'

print(
    "After",
    "new_version: {new_version}",
    f"major: {major}",
    f"minor: {minor}",
    f"patch: {patch}",
    sep = "\n"
)

version_file.write_text(new_version)

# Write new version to Github actions
if "GITHUB_OUTPUT" in os.environ:
    with open(os.environ["GITHUB_OUTPUT"], "a") as f:
        f.write(f"new_version={new_version}\n")
