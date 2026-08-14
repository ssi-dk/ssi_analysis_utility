## Description
<!-- Briefly describe the changes in this PR -->

## Type of Change
Add a label to this PR by commenting one of the following:
- `/build` - Minor documentation or general repository. These changes will not increment the software version upon successful merge
- `/patch` - Bugfixes and typos
- `/minor` - New features and software enhancements
- `/major` - Backwards incompatible changes (collaborators only)

**Note:** The label check waits 20 seconds before reading labels. Comment your label immediately after creating the PR to ensure the check passes.

## Checklist (Does NOT apply for changes to Github actions and documentations)
Before creating a Pull request for your changes, you can save time by completing the following steps in advance, on you local repository.
### Setup (Required for all PRs)
- [ ] I have installed the local branch (e.g. using `pip install .` in an appropriate virtual environment)

### For Adding or Removing Modules
- [ ] I have added my new module to / removed the old module from - the appropriate `.smk` file:
  - `PR_analysis.smk` for paired-end reads
  - `SR_analysis.smk` for single-end reads
  - `analysis.smk` for assembly-based analysis
- [ ] I have checked that the rule Output follows correct file naming conventions
- [ ] I have created/updated the appropriate species config in `mmaseq/src/mmaseq/config/species_configs/`
- [ ] I have added/removed the module alongside its configuration to both `test.yaml` and `all.yaml`
- [ ] I have included a sample and added its ftp link to `mmaseq/src/mmaseq/data/reads/reads.urls`
- [ ] I have updated the samplesheet at `mmaseq/src/mmaseq/data/samplesheet.tsv` (if needed)
- [ ] I have run a minimal test using `mmadeploy --test` and confirmed that my module completed succesfully

### For other pipeline Changes
- [ ] I have that all my changes work using `mmadeploy --update`

## Related Issues
<!-- Link any related issues: Closes #123 -->
