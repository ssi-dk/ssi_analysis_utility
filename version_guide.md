# How to
git tag <tagname (e.g. v1.2.3)> <commit>
git push origin --tags 

# Hierarki
Whenever a PR is merged into Dev or Main, the version must be bumped (by 1 e.g. 9 -> 10).

## Description
* **Major**: Versions that breaks backwards compatability
* **Minor**: Whenever something is merged into main without effecting backwards compatability
* **Patches**: Whenever any PR is merged into dev

Major: +.
Minor: .+
Patches: ..+

Major versions resets Minor and Patches versions
Minor versions resets Patches versions

## Examples
Branch `xyz` is merged into dev, version starts 1.4.8 -> After merge, version must be bumped to 1.4.9

Branch dev is merged into main, versions starts 1.4.9 -> After merge, version must be bumped to 1.5.0

Branch dev is merged into main and samplesheet format has now changed, which breaks backwards compatability, version starts 1.5.6 -> After merge, version must be bumped to 2.0.0


