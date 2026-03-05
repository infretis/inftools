# &infin;tools 🧰
👷‍♂️ Work in progress python package that provide various &infin;RETIS tools for tasks like setting up simulations and analyzing data. 📈

🤝 Contributions of tools and modifications are welcome – feel free to submit your enhancements! 🚀
## Installation
Inftools can be installed with:
```
pip install inftools
```
To install the latest development version:
1. Clone the repository
2. (Optional) Create a virtual environment
3. Install with:
   ```bash
   python -m pip install -e .
   ```
## Usage
```bash
inft -h
```

#### `inft center_periodic`

The command `inft center_periodic` can be run to re-center an .xyz trajectory to a certain atom index. For example, to center the first atom `C` at index 0 in the `co2.xyz` frame in [`examples/co2.xyz`](examples/co2.xyz) with cubic box length of 12.4138, the following command can be ran:

```bash
inft center_periodic -i co2.xyz -o co2_c.xyz -c 12.4138 -idx 0
```

## Development practice

Current CLI capabilities is powered by the [Typer](https://typer.tiangolo.com/) python library.

### Python files

Any function `function_name` in a python file `*.py` with the following structure

```python3
from typing import Annotated

def function_name(
    args1: Annotated[str, typer.Option("-args1", help="var1")],
    args2: Annotated[int, typer.Option("-args2", help="var2")],
    ):
    """
    Function description
    """
```

added to the folders [`inftools/exercises`](inftools/exercises), [`inftools/tistools`](inftools/tistools) and [`inftools/xyz`](inftools/xyz) will automatically be callable via the terminal as

```bash
inft function_name -args1 var1 -args2 var2
```

### Individual functions

Individual functions can also be added to the CLI library by including them in [`inftools/bin.py`](inftools/bin.py), as done for the wham function in the [`inftools/analysis`](inftools/analysis) folder:

```python3
from inftools.analysis.wham import wham

MAPPER["wham"] = wham
```

### Testing

The CLI functions are tested by utilizing the [pytest](https://docs.pytest.org/en/stable/) python library.

The tests included in the [`inftools/test`](inftools/) folder can be ran by running the following command in the terminal:

```bash
pytest
```

We currently do not have a testing policy but please do include tests for the included CLI functions for higher coverage.
