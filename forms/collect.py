import sys
import os
import shutil
import ast

count_opts = {}

# Init
if len(sys.argv) == 1 or len(sys.argv) > 3:
    print "Usage: name [q_suffix]"
    sys.exit(0)
if len(sys.argv) == 3:
    suffix = sys.argv[2]
    to_change = []
    for f in os.listdir("plots"):
        if f.endswith(".html"):
            filepath = os.path.join("plots", f)
            os.rename(filepath, filepath.replace(".html", "." + suffix + ".html"))
    for f in os.listdir("results"):
        if f.endswith(".dat"):
            filepath = os.path.join("results", f)
            os.rename(filepath, filepath.replace(".dat", "." + suffix + ".dat"))
name = sys.argv[1]
default_dir = "/data/FIREDRAKE_RESULTS/BENCH_SUITE_PLOTS_3"

# Create directory where to move stuff
default_dir = os.path.join(default_dir, name)
if not os.path.exists(default_dir):
    os.makedirs(default_dir)

# Move stuff
for f in os.listdir("plots"):
    filepath = os.path.join("plots", f)
    shutil.move(filepath, default_dir)
for f in os.listdir("results"):
    filepath = os.path.join("results", f)
    shutil.move(filepath, default_dir)
coffee_dir = "/tmp/coffee_dump/"
for d in os.listdir(coffee_dir):
    dirpath = os.path.join(coffee_dir, d)
    for f in os.listdir(dirpath):
        filepath = os.path.join(dirpath, f)
        if f.endswith(".out"):
            with open(filepath) as _f:
                lines = _f.readlines()
            for l in lines:
                # Parse the autotuner output and track used optimizations
                if l.startswith("***"):
                    opts = ast.literal_eval(l[29:-4])
                    if opts in count_opts:
                        count_opts[opts] += 1
                    else:
                        count_opts[opts] = 1
        elif not (f.endswith(".log") or f.endswith(".err")):
            os.remove(filepath)
            continue
        # Concatenate the files if already present
        if os.path.exists(os.path.join(default_dir, f)):
            with open(os.path.join(default_dir, f), "a") as outfile:
                with open(filepath, "r") as infile:
                    outfile.write(infile.read())
            os.remove(filepath)
        else:    
            shutil.move(filepath, default_dir)

# Record used optimizations
filepath = os.path.join(default_dir, "track_optimazations.txt")
with open(filepath, "a") as f:
    for k, v in count_opts.items():
        f.write("%s: %d\n" % (str(k), v))
