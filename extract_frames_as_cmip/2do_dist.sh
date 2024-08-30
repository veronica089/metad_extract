#!/bin/bash
initpath="/orozco/projects/E-Dent/VERONICA/DIMER_LARGER/mw_metad/analysis_phase2/mwf_cmipframes"
systems=(wt_Glu0 wt_Glu- mut_Glu0 mut_Glu-)
structures=(0 1 2)
frames=(-15.0 -13.5 -12.0 -10.5 -9.0 -7.5 -6.0 -4.5 -3.0 -1.5 0.0 10.5 9.0 7.5 6.0 4.5 3.0 1.5)
#frames=(-15.0 -13.5)
for sys in ${systems[@]}
   do
     cd $sys
     echo $sys
     ls
     for d1 in ${frames[@]} 
        do
          echo $d1
          for d2 in ${frames[@]}
            do
              echo $d2
	      ls frames_${d1}_${d2}
	      cd frames_${d1}_${d2}
              for st in ${structures[@]}
                do
                 mkdir -p frame$st
                 ls frame$st
                 cd frame$st
                 if [ -f $st.nc ]; then



        	     if [ $sys == 'wt_Glu-' ] || [ $sys == 'wt_Glu0' ]; then
cat << EOF > ptraj_dist.in
parm $initpath/$sys/$sys.prmtop
trajin $st.nc
distance 'ss_116_210' :210@OE1,OE2 :116@NZ out dist_${sys}_st${st}.dat
distance 'ss_515_210' :210@OE1,OE2 :515@NZ out dist_${sys}_st${st}.dat
distance 'ss1_459_210' :210@OE1,OE2 :459@CE out dist_${sys}_st${st}.dat
distance 'ss2_459_210' :210@OE1,OE2 :459@SD out dist_${sys}_st${st}.dat
distance 'ss116_501' :116@NZ :501@OD1,OD2 out dist_${sys}_st${st}.dat
distance 'ss152_153' :152@NZ :153@OE1,OE2 out dist_${sys}_st${st}.dat

go
quit
EOF
                    elif [ $sys == 'mut_Glu-' ] || [ $sys == 'mut_Glu0' ]; then
cat << EOF > ptraj_dist.in
parm $initpath/$sys/$sys.prmtop
trajin $st.nc
distance 'ss_116_210' :210@OE1,OE2 :116@NZ out dist_${sys}_st${st}.dat
distance 'ss_515_210' :210@OE1,OE2 :514@NZ out dist_${sys}_st${st}.dat
distance 'ss1_459_210' :210@OE1,OE2 :459@CE out dist_${sys}_st${st}.dat
distance 'ss2_459_210' :210@OE1,OE2 :459@SD out dist_${sys}_st${st}.dat
distance 'ss116_501' :116@NZ :500@OD1,OD2 out dist_${sys}_st${st}.dat
distance 'ss152_153' :152@NZ :153@OE1,OE2 out dist_${sys}_st${st}.dat
go
quit
EOF
	            fi
                
      cpptraj -i ptraj_dist.in
                 fi  #end if -f st.nc exist
                cd .. #exit from frame$st dir
                done
            cd ..   #exit from frames_d1_d2 dir
            done
        done
   cd ..
   done

