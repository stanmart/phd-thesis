[![checks](https://github.com/stanmart/phd-thesis/actions/workflows/ci.yml/badge.svg)](https://github.com/stanmart/phd-thesis/actions/workflows/ci.yml)
[![publish](https://github.com/stanmart/phd-thesis/actions/workflows/publish.yml/badge.svg)](https://github.com/stanmart/phd-thesis/actions/workflows/publish.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/licenses/MIT)

# PhD Dissertation

This repository contains code to compile my PhD dissertation, as well as the three papers that it consists of. The main outputs are the following:

 * [`dissertation.pdf`](stanmart.github.io/phd-thesis/dissertation.pdf): The full dissertation, the main topic of which is the bargaining between a few central and many fringe players.
 * [`theory.pdf`](stanmart.github.io/phd-thesis/theory.pdf): A paper with an abstract, theoretical treatment of the problem.
 * [`application.pdf`](stanmart.github.io/phd-thesis/application.pdf): An IO paper about a hybrid platform, which engages in bargaining with the fringe entrants.
 * [`experiment.pdf`](stanmart.github.io/phd-thesis/experiment.pdf): A lab experiment that examines the outcomes of bargaining and coalition formation in a setting with one indispensable and two smaller players.

> [!TIP]
> The code for the papers and the analysis were tracked in different repositories during the majority of their development.
> You can find the original repositories with the full histories at the following links:
> * [stanmart/ValueOfIntermediation](https://github.com/stanmart/ValueOfIntermediation)
> * [stanmart/unstructured-bargaining-analysis](https://github.com/stanmart/unstructured-bargaining-analysis)

## How to compile

The project is set up so that `pixi` installs required dependencies into a local virtual environment. The exception is latex, which has to be install manually, and `latexmk` must be on the search path (required texlive packages are listed in `tl_packages.txt`).

Other than these, simply [install `pixi`](https://pixi.sh/latest/#installation), and then you can use the following commands to compile the outputs:

 * `pixi run dissertation` to compile the dissertation
 * `pixi run papers` to compile three papers
 * `pixi run prepare-to-publish` to compile everything and place it in the `dist` folder

<details>
<summary>Other pixi tasks</summary>
The following commands are available to check the code:

 * `pixi run format` to format the Python code using `ruff`
 * `pixi run lint` to lint the Python code using `ruff`
 * `pixi run typecheck` to typecheck the Python code using `pyright`
 * `pixi run spell` to check the spelling using `codespell`
 * `pixi run check` to run all the checks

The following commands are available to create graphs of snakemake's execution plan:

 * `pixi run dag` to create a directed acyclic graph of the snakemake workflow
 * `pixi run filegraph` to create a file graph of the snakemake workflow
 * `pixi run rulegraph` to create a rule graph of the snakemake workflow

The following commands are used for the CI publish job:

 * `pixi run update-latex-deps` to collect the texlive packages needed for the project and write them to `tl_packages.txt`
 * `pixi run prepare-to-publish` to collect every output file and place them into the `gh-pages` folder

</details>

## Continuous integration

The project uses Github Actions to automatically compile the outputs and upload them to Github pages. First, the following checks must pass:

 * `ruff` for Python linting and formatting
 * `pyright` for Python type checking
 * `codespell` for spell checking

If they do, the `publish` job is triggered, which compiles the outputs and uploads them to the `gh-pages` branch. The results can be found at `stanmart.github.io/phd-thesis/{output}`.
