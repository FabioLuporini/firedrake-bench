## Poisson

pushd poisson

# strong scaling
for variant in firedrake dolfin; do
  for degree in 1 2 3; do
    # intra node
    for np in 1 3 6 12 24; do
      ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_s80_d${degree}_ -w 1:0 -r -n 1 --np $np -- ${variant}_poisson verbose=True size=80 degree=$degree -b -s
    # inter node
    for n in 2 4 8 16 32 64; do
      ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_s80_d${degree}_ -w 1:0 -r -n $n --np 24 -- ${variant}_poisson verbose=True size=80 degree=$degree -b -s
    done
  done
done

# Poisson weak scaling
for variant in firedrake dolfin; do
  for degree in 1 2 3; do
    # intra node
    for np in 1 3 6 12 24; do
      ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_w1k_d${degree}_ -w 0:30 -r -n 1 --np $np -- ${variant}_poisson verbose=True size=1000 weak=True degree=$degree -b -s
    # inter node
    for n in 2 4 8 16 32 64; do
      ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_w1k_d${degree}_ -w 0:30 -r -n $n --np 24 -- ${variant}_poisson verbose=True size=1000 weak=True degree=$degree -b -s
    done
  done
done

popd

## Wave

pushd wave

# strong scaling, intra node
for np in 1 3 6 12 24; do
  ../qsub.py -t ../firedrake.tpl -j f_s0.072_ -w 1:0 -r -n 1 --np $np -- firedrake_wave verbose=True scale=0.072 -b -s
# inter node
for n in 2 4 8 16 32; do
  ../qsub.py -t ../firedrake.tpl -j f_s0.072_ -w 1:0 -r -n $n --np 24 -- firedrake_wave verbose=True scale=0.072 -b -s
done

# weak scaling, intra node
for np in 1 3 6 12 24; do
  ../qsub.py -t ../firedrake.tpl -j f_w0.5_ -w 1:0 -r -n 1 --np $np -- firedrake_wave verbose=True scale=0.5 weak=True -b -s
# inter node
for n in 2 4 8 16; do
  ../qsub.py -t ../firedrake.tpl -j f_w0.5_ -w 1:0 -r -n $n --np 24 -- firedrake_wave verbose=True scale=0.5 weak=True -b -s
done

popd

## Cahn-Hilliard

pushd cahn_hilliard

# strong scaling
for variant in firedrake dolfin; do
  # intra node
  for np in 1 3 6 12 24; do
    ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_s2k_ch10_ -w 2:0 -r -n 1 --np $np -- ${variant}_cahn_hilliard verbose=True steps=10 size=2000 -b -s
  # inter node
  for n in 2 4 8 16 32 64; do
    ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_s2k_ch10_ -w 2:0 -r -n $n --np 24 -- ${variant}_cahn_hilliard verbose=True steps=10 size=2000 -b -s
  done
done

# weak scaling
for variant in firedrake dolfin; do
  # intra node
  for np in 1 3 6 12 24; do
    ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_w1k_ch10_ -w 2:0 -r -n 1 --np $np -- ${variant}_cahn_hilliard verbose=True steps=10 size=1000 weak=True -b -s
  # inter node
  for n in 2 4 8 16 32 64; do
    ../qsub.py -t ../${variant}.tpl -j ${variant:0:1}_w1k_ch10_ -w 2:0 -r -n $n --np 24 -- ${variant}_cahn_hilliard verbose=True steps=10 size=1000 weak=True -b -s
  done
done

popd
