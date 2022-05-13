`data` directory organization:

- `raw` contains raw files as they came from experiment results or an external source.
These files should generally not be modified once created.
Examples of such data are sequencing data from a core facility, expression data downloaded from a repository,

- `processed` may contain modified versions of the `raw` files, after steps like filtering, sorting, trimming, etc.

- `reference` contains files like genome annotations, mapping indices, or gene lists that are not specific to the project, but may be required for analysis.
