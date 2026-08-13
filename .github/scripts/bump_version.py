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

elif label == "patch":
    patch += 1


new_version = f'{major}.{minor}.{patch}'

text = re.sub(
    r'__version__ = ".*"',
    f'__version__ = "{new_version}"',
    text,
)

version_file.write_text(text)

# Write new version to Github actions
if "GITHUB_OUTPUT" in os.environ:
    with open(os.environ["GITHUB_OUTPUT"], "a") as f:
        f.write(f"new_version={new_version}\n")
