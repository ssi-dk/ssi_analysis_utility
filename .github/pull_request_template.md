## Description
<!-- Briefly describe the changes in this PR -->

## Type of Change
Add a label to this PR by commenting one of the following:
- `/build` - Minor documentation or general repository. These changes will not increment the software version upon successful merge
- `/patch` - Bugfixes and typos
- `/minor` - New features and software enhancements
- `/major` - Backwards incompatible changes (collaborators only)

**Note:** The label check waits 20 seconds before reading labels. Comment your label immediately after creating the PR to ensure the check passes.

## Checklist

### Setup (Required for all PRs)
- [ ] Installed local version: `pip install .` in a virtual environment

### For Adding or Removing Modules
- [ ] Added module to appropriate `.smk` file:
  - `PR_analysis.smk` for paired-end reads
  - `SR_analysis.smk` for single-end reads
  - `analysis.smk` for assembly-based analysis
- [ ] Output follows correct file naming conventions
- [ ] Created/updated species config in `mmaseq/src/mmaseq/config/species_configs/`
- [ ] Added module configuration to both `test.yaml` and `all.yaml`
- [ ] Added public test dataset to `mmaseq/src/mmaseq/data/samplesheet.tsv` (if needed)
- [ ] Tested with `mmadeploy --test`

### For Other Changes
- [ ] Verified changes with `mmadeploy --update`

## Related Issues
<!-- Link any related issues: Closes #123 -->
