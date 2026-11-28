"""Write the CI papermill parameters for one notebook to a YAML file.

Usage: python .github/scripts/ci_params.py <notebook> <output.yml>

Looks up <notebook> (path relative to the repository root) in
.github/ci-parameters.yml and writes its parameters to <output.yml>, to be
passed to papermill with `-f <output.yml>`. Notebooks that are not listed get
an empty parameter set, so they run with their own defaults.

Fails if a listed notebook does not declare every requested parameter in a
cell tagged `parameters`: papermill would then inject the values at the top of
the notebook, where they would be silently overwritten.
"""

import os
import sys

import papermill
import yaml

CONFIG = os.path.join(os.path.dirname(__file__), "..", "ci-parameters.yml")


def main(notebook, output):
    with open(CONFIG) as f:
        config = yaml.safe_load(f) or {}
    config = {os.path.normpath(path): params for path, params in config.items()}

    params = config.get(os.path.normpath(notebook)) or {}

    if params:
        declared = papermill.inspect_notebook(notebook)
        missing = sorted(set(params) - set(declared))
        if missing:
            sys.exit(
                f"{notebook}: parameters {missing} are set in .github/ci-parameters.yml "
                "but not defined in a cell tagged 'parameters'"
            )
        print(f"{notebook}: running with CI parameters {params}")
    else:
        print(f"{notebook}: running with notebook defaults")

    with open(output, "w") as f:
        yaml.safe_dump(params, f)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
