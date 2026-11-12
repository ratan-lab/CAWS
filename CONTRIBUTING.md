# Contributing to CAWS

Thanks for your interest in improving CAWS. Bug reports, questions, documentation fixes and code contributions are all welcome.

## Reporting bugs and asking questions

Open an issue at https://github.com/ratan-lab/CAWS/issues. For bugs, please include:

- what you ran (the `snakemake` command line) and what you expected to happen
- your `config.json` and samplesheet (file paths can be redacted)
- the log of the failing rule from `logs/<rule>/`, and the error printed by Snakemake
- the Snakemake and conda versions (`snakemake --version`, `conda --version`) and the CAWS commit or release

## Contributing code

1. Fork the repository and create a branch from `main`.
2. Make your change. Keep each pull request focused on one fix or feature.
   - Tools used by rules go in a conda environment under `workflow/envs/`, with pinned versions.
   - New configuration options need a default in `config.json` and an entry in the README configuration table.
   - R code follows tidyverse style with the native `|>` pipe.
3. Run the test dataset locally (see "Testing the Installation" in the README):
   ```bash
   conda run -n caws-setup bash .test/setup.sh
   snakemake -s workflow/Snakefile --configfile .test/run/config.json --use-conda --cores 4
   bash .test/check_outputs.sh
   ```
   If your change adds an output, add a check for it to `.test/check_outputs.sh`.
4. Open a pull request describing what changed and why. GitHub Actions runs the same test on every pull request.

## Use of AI tools

If you used AI tools to write a substantial part of your contribution, say so in the pull request. You are responsible for reviewing and testing all code you submit.
