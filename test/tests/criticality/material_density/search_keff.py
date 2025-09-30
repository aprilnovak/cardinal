import openmc
import argparse

from make_openmc_model import make_pu_cube

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target", type=float, default=1.0)
args = parser.parse_args()

crit_density, guesses, keffs = openmc.search_for_keff(
        make_pu_cube,
        target=args.target,
        bracketed_method="brentq",
        bracket=[1, 30],
        tol=1e-5,
        print_iterations=True)

for guess, keff in zip(guesses, keffs):
    print( guess, keff)

print("final density ", crit_density)

