################################
#!/bin/bash
function pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
        sleep 1
    done
} 

# for j in 8 13 16 21 27 30 34 39
# for j in 41 44 49 52 58 63 66 71
# for j in 76 81 84 90 95 99 102 104
# for j in 109 112 117 105 106 107
for l in 78
do
for k in 1
do
  echo $k
  pwait 8
  source Melt_test_jobscript_clim_var.sh &
    sleep 30
  echo 'test'
done
done
################################

