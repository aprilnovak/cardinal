import os
import re
import shutil

mdot_values = [1.0, 1.2, 10.0]
fuel_values = [0.5, 2.0, 20.0]

os.system('rm -rf results*')


def modify_common(mdot, fuel):
    with open("common.i") as f:
        text = f.read()

    text = re.sub(r"mdot\s*=\s*.*", f"mdot = {mdot}", text)

    text = re.sub(r"fuel_conductivity\s*=\s*.*", f"fuel_conductivity = {fuel}", text)

    with open("common.i", "w") as f:
        f.write(text)


for mdot in mdot_values:
    for fuel in fuel_values:

        modify_common(mdot, fuel)

        print(f"Running mdot={mdot}, fuel={fuel}")
        for file in os.listdir("."):
            if file.startswith("openmc_out_axial_power"):
                os.remove(file)

        os.system("cardinal-opt -i openmc.i --n-threads=10")

        folder = f"results_m{mdot}_f{fuel}"
        os.makedirs(folder, exist_ok=True)

        for file in os.listdir("."):
            if file.startswith("openmc_out_axial_power"):
                shutil.move(file, os.path.join(folder, file))
