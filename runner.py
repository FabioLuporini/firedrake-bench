#!/usr/bin/env python
import subprocess
import os
import sys

try:
    with open("runner_progress", "r") as f:
        progress = int(f.read())
except:
    progress = 0

count = 0
maxcount = 16


def qsub(*args):
    global count
    try:
        if count >= progress:
            print "Starting job %d" % count
            print "submitting %s" % " ".join(("qsub",) + args)
            subprocess.check_call(("../qsub.py",) + args)
        count += 1
        if count >= progress + maxcount:
            print "Maximum jobs achieved"
            raise Exception
    except:
        with open("../runner_progress", "w") as f:
            f.write(str(count))
        sys.exit(0)


os.chdir("poisson")

# strong scaling
for variant in ("firedrake", "dolfin"):
    for degree in (1, 2, 3):
        # intra node
        for np in (1, 3, 6, 12, 24):
            qsub("-t", "../%s.tpl" % variant, "-j", "%s_s80_d%s_" %
                 (variant[:2], degree), "-w", "1:0", "-r", "-n", "1",
                 "--np", str(np), "--", "%s_poisson" % variant,
                 "verbose=True", "size=80", "degree=%s" % degree, "-b",
                 "-s")
        # inter node
        for n in (2, 4, 8, 16, 32, 64):
            qsub("-t", "../%s.tpl" % variant, "-j", "%s_s80_d%s_" %
                 (variant[:2], degree), "-w", "1:0", "-r", "-n", str(n),
                 "--np", "24", "--", "%s_poisson" % variant, "verbose=True",
                 "size=80", "degree=%s" % degree, "-b", "-s")


# Poisson weak scaling
for variant in ("firedrake", "dolfin"):
    for degree in (1, 2, 3):
        # intra node
        for np in (1, 3, 6, 12, 24):
            qsub("-t", "../%s.tpl" % variant, "-j", "%s_w1k_d%s_" %
                 (variant[:2], degree), "-w", "0:30", "-r", "-n", "1",
                 "--np", str(np), "--", "%s_poisson" % variant,
                 "verbose=True", "size=1000", "weak=True", "degree=%s" %
                 degree, "-b", "-s")
        # inter node
        for n in (2, 4, 8, 16, 32, 64):
            qsub("-t", "../%s.tpl" % variant, "-j", "%s_w1k_d%s_" %
                 (variant[:2], degree), "-w", "0:30", "-r", "-n",
                 str(n), "--np", "24", "--", "%s_poisson" % variant,
                 "verbose=True", "size=1000", "weak=True", "degree=%s"
                 % degree, "-b", "-s")


os.chdir("..")

# Wave

os.chdir("wave")

# strong scaling, intra node
for np in (1, 3, 6, 12, 24):
    qsub("-t", "../firedrake.tpl", "-j", "f_s0.072_", "-w", "1:0", "-r",
         "-n", "1", "--np", str(np), "--", "firedrake_wave", "verbose=True",
         "scale=0.072", "-b", "-s")
# inter node
for n in (2, 4, 8, 16, 32):
    qsub("-t", "../firedrake.tpl", "-j", "f_s0.072_", "-w", "1:0", "-r",
         "-n", str(n), "--np", "24", "--", "firedrake_wave", "verbose=True",
         "scale=0.072", "-b", "-s")


# weak scaling, intra node
for np in (1, 3, 6, 12, 24):
    qsub("-t", "../firedrake.tpl", "-j", "f_w0.5_", "-w", "1:0", "-r", "-n",
         "1", "--np", str(np), "--", "firedrake_wave", "verbose=True",
         "scale=0.5", "weak=True", "-b", "-s")
# inter node
for n in (2, 4, 8, 16):
    qsub("-t", "../firedrake.tpl", "-j", "f_w0.5_", "-w", "1:0", "-r", "-n",
         str(n), "--np", "24", "--", "firedrake_wave", "verbose=True",
         "scale=0.5", "weak=True", "-b", "-s")


os.chdir("..")

# Cahn-Hilliard

os.chdir("cahn_hilliard")

# strong scaling
for variant in ("firedrake", "dolfin"):
    # intra node
    for np in (1, 3, 6, 12, 24):
        qsub("-t", "../%s.tpl" % variant, "-j",
             "%s_s2k_ch10_" % variant[:2], "-w", "2:0", "-r", "-n", "1", "--np",
             str(np), "--", "%s_cahn_hilliard" % variant, "verbose=True",
             "steps=10", "size=2000", "-b", "-s")
    # inter node
    for n in (2, 4, 8, 16, 32, 64):
        qsub("-t", "../%s.tpl" % variant, "-j",
             "%s_s2k_ch10_" % variant[:2], "-w", "2:0", "-r", "-n", str(n), "--np",
             "24", "--", "%s_cahn_hilliard" % variant, "verbose=True",
             "steps=10", "size=2000" "-b", "-s")


# weak scaling
for variant in ("firedrake", "dolfin"):
    # intra node
    for np in (1, 3, 6, 12, 24):
        qsub("-t", "../%s.tpl" % variant, "-j", "%s_w1k_ch10_" %
             variant, "-w", "2:0", "-r", "-n", "1", "--np", str(np),
             "--", "%s_cahn_hilliard" % variant, "verbose=True",
             "steps=10", "size=1000", "weak=True", "-b", "-s")
    # inter node
    for n in (2, 4, 8, 16, 32, 64):
        qsub("-t", "../%s.tpl" % variant, "-j", "%s_w1k_ch10_" %
             variant, "-w", "2:0", "-r", "-n", str(n), "--np", "24", "--",
             "%s_cahn_hilliard" % variant, "verbose=True", "steps=10",
             "size=1000", "weak=True", "-b", "-s")


os.chdir("..")
