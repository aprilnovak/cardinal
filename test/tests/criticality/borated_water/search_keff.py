import openmc

from make_openmc_model import make_uo2_bwater_box

crit_boron, guesses, keffs = openmc.search_for_keff(
        make_uo2_bwater_box,
        bracketed_method="brentq",
        bracket=[0, 1000],
        tol=1e-3,
        print_iterations=True)

for guess, keff in zip(guesses, keffs):
    print( guess, keff)

print("final boron concentration ", crit_boron)