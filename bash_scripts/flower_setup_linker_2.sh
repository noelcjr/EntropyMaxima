#!/bin/bash

# Example run:
#bash setup_linker_2.sh DDDDK DDDDK 2hiu_insulin_no_disulfides/2hiu_1rr_no_disulf.csv "1,1,A,Ndir" "1,2,B,Ndir" 1brs_barnase_barstar/1brs_comp_1rr.csv 2hiu_a_ddddk_b_ddddk.csv

# This script adds two linkers to the centered structure, one in amino acid, model, chain and N-terminal direction.
# This scripts works for insulin only, and it needs to be made more general purpuse by letting the user specify the location
# and direction that the linker is going to be added
if [ -z "$1" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
fi
if [ -z "$2" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
fi
if [ -z "$3" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
fi
if [ -z "$4" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
fi
if [ -z "$5" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
fi
if [ -z "$6" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
fi
if [ -z "$7" ];
then
   echo usage: $0 Linker_String_single_letter_1 Linker_String_single_letter_2 CSV_input_center_structure place_first_link place_second_link CSV_input_rotate_structure output_file_name
   exit
else
   lnk1=$1
   lnk2=$2
   center="../"$3
   plc_lnk1=$4
   plc_lnk2=$5
   rotate="../"$6
   output=$7
   echo $lnk1"_"$lnk2" "$center" "$plc_lnk1" "$plc_lnk2" "$rotate" "$output
   dir=$lnk1"_"$lnk2
   if [ -d $dir ];
   then
      #rm -f -r $dir It is too risky to remove. It might delete a lot of information. Remove manully if needed. 
      echo "Warning: Folder for linker already exist. Script will exit without any changes to the folder corresponding to this linker."
      exit
   fi
   mkdir $dir
   #cp *pdb $dir
   #cp *csv $dir
   cd $dir
   cp $charmm"setup_one.inp" .
   # cp $charmm"minimize.inp" .
   # cp $charmm"gbsw_AB.inp" .
   # The following two lines are already done for DDDDK  so they are comented out.
   # add_residues.py" --apn "1,1,A,Ndir" --res $dir --inp "2hiu.csv" --out "2hiu.csv"
   # add_residues.py" --apn "1,2,B,Ndir" --res $dir --inp "2hiu.csv" --out "2hiu.csv"
   add_residues.py --apn $plc_lnk1 --res $lnk1 --inp $center --out $output
   add_residues.py --apn $plc_lnk2 --res $lnk2 --inp $output --out $output
   
   currentfile=${output%.*v}"_1"   
   reduce -HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000 $currentfile".pdb" > $currentfile"r.pdb"
   pdb.py --prepare --pdbin1 $currentfile"r.pdb" --crdout $currentfile"r.crd" --seqfix yes
   #perl -pi -e 's/generate A first none last none setup/generate A first none last CTER setup/g' setup_one.inp
   #perl -pi -e 's/generate B first none last none setup/generate B first none last CTER setup/g' setup_one.inp
   perl -pi -e 's/INFILE/'$currentfile'r/g' setup_one.inp
   perl -pi -e 's/OUTFILE/'$currentfile'rr/g' setup_one.inp
   charmm_xxlarge < setup_one.inp > setup_one.out
   rm *xplo.psf
   rm *.ic
   rm *out
   rm *SEQ
   rm *INP
   rm $currentfile"r.pdb"
   rm $currentfile".pdb"
   perl -pi -e 's/'$currentfile'rr/OUTFILE/g' setup_one.inp
   perl -pi -e 's/'$currentfile'r/INFILE/g' setup_one.inp
   #perl -pi -e 's/generate A first none last CTER setup/generate A first none last none setup/g' setup_one.inp
   #perl -pi -e 's/generate A first none last CTER setup/generate A first none last none setup/g' setup_one.inp
   pdb.py --fixpdb --pdbin2 $currentfile"rr.pdb"
   
   currentfile=${rotate%.*v}
   cp $currentfile"_1.pdb" .
   currentfile=$(basename $rotate)
   currentfile=${currentfile%.*v}"_1"
   # TODO currentfile has an extra r this is an exception so it is treated specially here
   # currentfile=${currentfile::(-2)}
   echo "currentfile "$currentfile
   reduce -HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000 $currentfile".pdb" > $currentfile"r.pdb"
   pdb.py --prepare --pdbin1 $currentfile"r.pdb" --crdout $currentfile"r.crd" --seqfix yes
   perl -pi -e 's/generate A first none last none setup/generate C first none last none setup/g' setup_one.inp
   perl -pi -e 's/generate B first none last none setup/generate F first none last none setup/g' setup_one.inp
   perl -pi -e 's/A.SEQ/C.SEQ/g' setup_one.inp
   perl -pi -e 's/B.SEQ/F.SEQ/g' setup_one.inp
   perl -pi -e 's/A_FIXRES.INP/C_FIXRES.INP/g' setup_one.inp
   perl -pi -e 's/B_FIXRES.INP/F_FIXRES.INP/g' setup_one.inp
   perl -pi -e 's/INFILE/'$currentfile'r/g' setup_one.inp
   perl -pi -e 's/OUTFILE/'$currentfile'rr/g' setup_one.inp
   charmm_xxlarge < setup_one.inp > "setup_one"$i".out"
   rm *xplo.psf
   rm *.ic
   rm *out
   rm *SEQ
   rm *INP
   rm $currentfile"r.pdb"
   rm $currentfile".pdb"
   rm *crd
   rm *psf
   #perl -pi -e 's/'$currentfile'rr/OUTFILE/g' setup_one.inp
   #perl -pi -e 's/'$currentfile'r/INFILE/g' setup_one.inp
   #perl -pi -e 's/generate A first ACE last none setup/generate A first none last none setup/g' setup_one.inp
   #perl -pi -e 's/generate B first ACE last none setup/generate B first none last none setup/g' setup_one.inp
   pdb.py --fixpdb --pdbin2 $currentfile"rr.pdb"
   cd ..
fi
