#!/bin/csh

set t=1
set t_stop=50

cd /home/smith292/sensitivity/data

while ( $t <= $t_stop )

	ln -s /lustre/work/smith292/NWP/ramp_2011_12_03/mem_${t}/wrfout_d01_2011-12-02_18:00:00 wrfout_full_mem.${t}

	@ t++
end
