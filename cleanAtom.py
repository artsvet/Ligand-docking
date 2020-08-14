import os

to_walk = os.getcwd() + '\\targets'

for root, dirs, files in os.walk(to_walk):
    for file in files:
        f_in = root + '\\' + file
        f_out = f_in[:-4] + '.clean.pdb'

        with open(f_in, "r") as f:
            lines = f.readlines()

        with open(f_out, "w") as f:
            for line in lines:
                if line.strip("\n")[:4] == "ATOM":
                    f.write(line)
