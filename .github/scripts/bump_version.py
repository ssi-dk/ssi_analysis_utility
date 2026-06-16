import re
import sys
from pathlib import Path

label = sys.argv[1]

version_file = Path("src/mmaseq/_version.py")

text = version_file.read_text()

m = re.search(r'__version__ = "(\d+)\.(\d+)\.(\d+)\+(\d+)"', text)

major, minor, patch, build = map(int, m.groups())

if label == "feature":
    minor += 1
    patch = 0
    build = 0

elif label == "patch":
    patch += 1
    build = 0

elif label == "build":
    build += 1

new_version = f'{major}.{minor}.{patch}+{build}'

text = re.sub(
    r'__version__ = ".*"',
    f'__version__ = "{new_version}"',
    text,
)

version_file.write_text(text)

print(new_version)
