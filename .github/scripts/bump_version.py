import re
import sys
from pathlib import Path

label = sys.argv[1]

version_file = Path("src/mmaseq/__version__.py")

text = version_file.read_text()

m = re.search(r'__version__ = "(\d+)\.(\d+)\.(\d+)"', text)

major, minor, patch = map(int, m.groups())

if label == "major":
    major += 1
    minor = 0
    patch = 0
elif label == "feature":
    minor += 1
    patch = 0
    # Determine every tens of minor version increments (10, 20, 30 etc.)
    if (minor % 10) == 0:
        # Set patch version to 1 directly
        patch = 1
    build = 0

elif label == "patch":
    patch += 1

new_version = f'{major}.{minor}.{patch}'

text = re.sub(
    r'__version__ = ".*"',
    f'__version__ = "{new_version}"',
    text,
)

version_file.write_text(text)

print(new_version)
