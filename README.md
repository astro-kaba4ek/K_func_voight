# Calculate the core function `K`, `K_0`, `L` and `L_0` for the Voigt absorption profile at different values of parameters `a`, `\beta` and `\tau`.

The program for calculating the core functions `K`, `K_0`, `L` and `L_0` for the Voigt absorption profile depending on the parameters `a`, `\beta` and `\tau`. The program is not able to calculate values for large (depending on `a` and `\beta`) `\tau`, due to the peculiarities of the formulas used and the overflow of variables. Use carefully when `\tau` is greater than `~10`. Accuracy `~1e-5`.

Program written as part of the course "Special Workshop in Theoretical Astrophysics" under the guidance of Dr. Nagirner. Author S.I. Laznevoi.

## Files
- **main.py** — main program
- **utils.py** — module with auxiliary functions (reading arguments, saving results, etc.)
- **analytical_functions.py** — module with "analytical" core functions
- **K_nl_funcs.py** — module with "numerical" core functions
- **calculation_functions.py** — module with auxiliary functions for numerical calculation
- **constants.py** — module with constants used
- **input.txt** — example of an input file
- **kvoigt.pdf** — a file with theoretical information on which this code is based

## The content of the input.txt
1. **marker_a** — if not 0 (k), then an array of k values, if 0 - from a_1 to a_N N evenly spaced points in the interval
2. **a** — `a` task parameter
3. **marker_b** — if not 0 (k), then an array of k values, if 0 - from b_1 to b_N N evenly spaced points in the interval
4. **b** — `\beta` the ratio of absorption coefficients in the continuum and the center of the line. Task parameter
5. **marker_t** — if not 0 (k), then an array of k values, if 0 - from t_1 to t_N N evenly spaced points in the interval
6. **t** — `\tau` the optical distance. Task parameter

## Usage
To run the program, you can use the command 
```bash
python3 ./main.py --input input.txt
``` 
You can use other flags in conjunction with `--input` or instead of it. To learn more, use the 
```bash
python3 ./main.py --help
```
## Additional information
- The following packages are required for the program to work:
	- **scipy**
	- **numpy**
	- **math**
	- **matplotlib**
	- **argparse**
	- **subprocess**
