#!/bin/bash
#
# Javier Cerezo. May 2011
#
##################################################################################
# DESCRIPTION
# Interface g09-gromacs thruogh external kewword (in g09/g16)
# to be used as any other method.
# Gromacs (MM) can provide Energy, Dipole, Gradients and Hessian. 
# But not polarizability and dipole derivative
#
# G09 passes the coordinates to gromacs through an ASCII
# file. First this file is read and a trr is generated (no
# loss of precision as if a gro/pdb file was generated).
# Then a mdrun step is performed (the top file should be
# provided in command line or through the parameter file). 
# Required data extracted using gromacs tools, such as g_energy 
# (energy) g_dipoles (dipole) and dump (forces from trr since 
# g_traj is broken for small molecules for version<4.5.3). The 
# hessian is obtained after a run with integrator=nm and dump 
# over the hessian.mtx file. 
#
# The most relevant argument to the program is the name of the
# topology file, which can be provide in the route section as
# -p topfilename
#
# A parameter file can also be provided (e.g. parameters.dat)
# specifying the files required for the gromacs runs, including:
# .top (topology), .gro (structure, needed to generate the tpr 
# and to provide support for the output in special UA/QMMM 
# calculations, see below), .mdp files (both SP and NM calcu-
# lations). Separate gro/top files shuold be provide for all
# layers (for ONIOM calculations)
#
# UA/QMMM calculations are supported (no hessian available). 
# A mapping file is used to traslate UA<->AA models (e.g. aa2ua.dat).
# Forces in UA particles are redistributed over the AA equi-
# valent particles. The MM part is indicated in Gaussian with
# a gosht atom. This is traslated to Gromacs with the aa2ua.dat
# file, indicating 0 for the ghost atom.
#
# In the Gaussian input line the method that uses GROMACS
# is specified as:
#
#   External="gromacs_gau.sh [options]"
#
# where [options] may include route section options (-flag value)
# or the name of the parameter file
#
# The interface uses two additionar utilities to transform data from
# gaussian and gromacs formats, that should be included in the package
# along with this script.
#
##################################################################################

#========================================
# GROMACS CONFIGURATION SECTION
#========================================
# # Version 5.x 
# gmxcall="gmx_mpi_d"
# gmxprefix=""
# gmxsufix=""
# dumprefix=""
# addmdpgroup=true
# Version 4.5.x
gmxcall=""
gmxprefix="g_"
gmxsufix="_d"
dumprefix="gmx"
addmdpgroup=false
#
# Only 
export OMP_NUM_THREADS=1 # mcp gmx has problems if threads and MPI are mixed
#========================================
# END OF GROMACS CONFIGURATION SECTION
#========================================

# Defaults (if not defined in the file)
mdpfile='generate'
mdpfile_min='generate'
mdpfile_nm='generate'
topfile='topol.top'
topfile_R='default'
topfile_M='defalt'
topfile_S='default'
grofile='gaussian'
grofile_R='default'
grofile_M='default'
grofile_S='default'
save_intermediate='false'
nproc_gmx='1'
aa2ua="NO"
job_type="SinglePoint"
parfile="none"

##################################################################################
# INPUT DATA:
# First commands are given afeter the program name in External
# Last commans are appended by gaussian in this order:
#  * layer to be calculated (R, M or S)
#  * input file name
#  * output file name
#  * message file name
#  * another file
#  * another file
#
# Now read input except the arguments appended by gaussian
(( nopts = $# - 6 ))
(( i=0 ))
while (( i<nopts )); do
    case $1 in
     -c         ) shift; (( i++ )); grofile=$1       ;;
     -p         ) shift; (( i++ )); topfile=$1       ;;
     -clean     ) save_intermediate=false            ;;
     -noclean   ) save_intermediate=true             ;;
     *          ) parfile=$1                         ;;
    esac
    shift; (( i++ ))
done
# Take gaussian args
layer=$1 ; shift
input=$1 ; shift
output=$1; shift
msg=$1   ; shift
##################################################################################
# Usamos el gaussian ID como identificador del trabajo
# Añadimos _tmp para diferenciarlo de los ficheros de Gaussian
label=${input##/*/}
label=${label%.*}_tmp


##################################################################################
# P R E P A R E   M M    R U N
##################################################################################
echo " " >> gmx_${label}.log
echo " " >> gmx_${label}.log
echo "** NEW GROMACS STEP **" >> gmx_${label}.log
echo " " >> gmx_${label}.log

if [ "$parfile" != "none" ]; then
    echo "Reading parameters from: $parfile" >> gmx_${label}.log
# GET PARAMETERS FROM PROVIDED FILE 

    # Preprocess param to get rid out of commentaries
    egrep -v "^[\ ]*;" $param > procesed.par

    # Get the parameters if defined
    mdpfile_deff=$(egrep "mdpfile[\ ]*=" procesed.par)
    if [ "x$mdpfile_deff" != "x" ]; then
       mdpfile=${mdpfile_deff##*=} 
       mdpfile=${mdpfile%;*}
    fi
    mdpfile_nm_deff=$(egrep "mdpfile_nm[\ ]*=" procesed.par)
    if [ "x$mdpfile_nm_deff" != "x" ]; then
       mdpfile_nm=${mdpfile_nm_deff##*=} 
       mdpfile_nm=${mdpfile_nm%;*}
    fi
    mdpfile_min_def=$(egrep "mdpfile_min[\ ]*=" procesed.par)
    if [ "x$mdpfile_min_def" != "x" ]; then
       mdpfile_min=${mdpfile_min_def##*=}
       mdpfile_min=${mdpfile_min%;*}
    fi
    # If only one layer, we can only define a defauilt top
    topfile_deff=$(egrep "topfile[\ ]*=" procesed.par)
    if [ "x$topfile_deff" != "x" ]; then
        topfile=${topfile_deff##*=}
        topfile=${topfile%;*}
    fi
    topfile_R_deff=$(egrep "topfile_R[\ ]*=" procesed.par)
    if [ "x$topfile_R_deff" != "x" ]; then
        topfile_R=${topfile_R_deff##*=}
        topfile_R=${topfile_R%;*}
    fi
    topfile_M_deff=$(egrep "topfile_M[\ ]*=" procesed.par)
    if [ "x$topfile_M_deff" != "x" ]; then
        topfile_M=${topfile_M_deff##*=}
        topfile_M=${topfile_M%;*}
    fi
    topfile_S_deff=$(egrep "topfile_S[\ ]*=" procesed.par)
    if [ "x$topfile_S_deff" != "x" ]; then
        topfile_S=${topfile_S_deff##*=}
        topfile_S=${topfile_S%;*}
    fi
    grofile_deff=$(egrep "grofile[\ ]*=" procesed.par)
    if [ "x$grofile_deff" != "x" ]; then
        grofile=${grofile_deff##*=}
        grofile=${grofile%;*}
    fi
    grofile_R_deff=$(egrep "grofile_R[\ ]*=" procesed.par)
    if [ "x$grofile_R_deff" != "x" ]; then
        grofile_R=${grofile_R_deff##*=}
        grofile_R=${grofile_R%;*}
    fi
    grofile_M_deff=$(egrep "grofile_M[\ ]*=" procesed.par)
    if [ "x$grofile_M_deff" != "x" ]; then
        grofile_M=${grofile_M_deff##*=}
        grofile_M=${grofile_M%;*}
    fi
    grofile_S_deff=$(egrep "grofile_S[\ ]*=" procesed.par)
    if [ "x$grofile_S_deff" != "x" ]; then
        grofile_S=${grofile_S_deff##*=}
        grofile_S=${grofile_S%;*}
    fi
    save_deff=$(egrep "save_intermediate[\ ]*=" procesed.par)
    if [ "x$save_deff" != "x" ]; then
        save_intermediate=$(echo $save_deff | egrep -i "(true|false)" -o)
    fi
    nproc_deff=$(egrep "nproc_gmx[\ ]*=" procesed.par)
    if [ "x$nproc_deff" != "x" ]; then
        nproc_gmx=$(echo $nproc_deff | egrep "[0-9]+" -o)
    fi
    aa2ua_def=$(egrep "aa2ua_file[\ ]*=" procesed.par)
    if [ "x$aa2ua_def" != "x" ]; then
        aa2ua=${aa2ua_def##*=}
        aa2ua=${aa2ua%;*}
    fi
    job_type_def=$(egrep "job_type[\ ]*=" procesed.par)
    if [ "x$job_type_def" != "x" ]; then
        job_type=${job_type_def##*=}
        job_type=${job_type%;*}
    fi
fi

# Set defaults when no value is given
if [ "$topfile_R" == "default" ]; then topfile_R=$topfile; fi
if [ "$topfile_M" == "default" ]; then topfile_M=$topfile; fi
if [ "$topfile_S" == "default" ]; then topfile_S=$topfile; fi
if [ "$grofile_R" == "default" ]; then grofile_R=$grofile; fi
if [ "$grofile_M" == "default" ]; then grofile_M=$grofile; fi
if [ "$grofile_S" == "default" ]; then grofile_S=$grofile; fi

########################################
# Generate mdp files if needed
########################################
if [ "$mdpfile" == "generate" ]; then
mdpfile="SinglePoint_${label}.mdp"
echo "Generating $mdpfile" >> gmx_${label}.log
cat<<EOF >$mdpfile
integrator               = md
nsteps                   = 0

nstfout                  = 1
nstxout                  = 1

nstlist                  = 0
nstype                   = simple
pbc                      = no
rlist                    = 0.0

coulombtype              = Cut-off
rcoulomb                 = 0.0

vdwtype                  = Cut-off
rvdw                     = 0.0

EOF
if $addmdpgroup; then
    echo "cutoff-scheme            = group" >> $mdpfile
fi
fi
if [ "$mdpfile_min" == "generate" ]; then
mdpfile_min="Minim_${label}.mdp"
echo "Generating $mdpfile_min" >> gmx_${label}.log
cat<<EOF >$mdpfile_min
integrator               = l-bfgs
nsteps                   = 100000
emtol                    = 0.000001

nstenergy                = 1
nstcalcenergy            = 1

nstlist                  = 0
nstype                   = simple
pbc                      = no
rlist                    = 0.0

coulombtype              = Cut-off
rcoulomb                 = 0.0

vdwtype                  = Cut-off
rvdw                     = 0.0

EOF
if $addmdpgroup; then
    echo "cutoff-scheme            = group" >> $mdpfile_min
fi
fi
if [ "$mdpfile_nm" == "generate" ]; then
mdpfile_nm="NormalModes_${label}.mdp"
echo "Generating $mdpfile_nm" >> gmx_${label}.log
cat<<EOF >$mdpfile_nm
integrator               = nm
nsteps                   = 1
emtol                    = 0.0001

nstenergy                = 1
nstcalcenergy            = 1

nstlist                  = 0
nstype                   = simple
pbc                      = no
rlist                    = 0.0

coulombtype              = Cut-off
rcoulomb                 = 0.0

vdwtype                  = Cut-off
rvdw                     = 0.0

EOF
if $addmdpgroup; then
    echo "cutoff-scheme            = group" >> $mdpfile_nm
fi
fi
echo " " >> gmx_${label}.log
    
# Write info to gaussian log file
echo "gmxMM   grofile_R: $grofile_R" >>$msg
echo "gmxMM   grofile_M: $grofile_M" >>$msg
echo "gmxMM   grofile_S: $grofile_S" >>$msg

# Generate trr files if staring from gromacs formats (gro,g96)
if [ "$grofile_R" != "gaussian" ]; then
if [ -f $grofile_R ]; then 
    $gmxcall trjconv$gmxsufix -f $grofile_R -o ${label}_R.trr &>>gmx_${label}.log
fi
fi
if [ "$grofile_M" != "gaussian" ]; then
if [ -f $grofile_M ]; then
    $gmxcall trjconv$gmxsufix -f $grofile_M -o ${label}_M.trr &>>gmx_${label}.log
fi
fi
if [ "$grofile_S" != "gaussian" ]; then
if [ -f $grofile_S ]; then
    $gmxcall trjconv$gmxsufix -f $grofile_S -o ${label}_S.trr &>>gmx_${label}.log
fi
fi


# Set the grofile and topology to use
case $layer in
 R ) grofile=$grofile_R; trrfile=${label}_R.trr; topfile=$topfile_R ;;
 M ) grofile=$grofile_M; trrfile=${label}_M.trr; topfile=$topfile_M ;;
 S ) grofile=$grofile_S; trrfile=${label}_S.trr; topfile=$topfile_S
esac

cat<<EOF >>$msg
gmxMM  GROMACS PARAMETERS
gmxMM   mdpfile = $mdpfile
gmxMM   mdpfile_nm = $mdpfile_nm
gmxMM   mdpfile_min = $mdpfile_min
gmxMM   topfile = $topfile
gmxMM   struct file = $trrfile
gmxMM   save_intermediate = $save_intermediate
gmxMM   nproc_gmx = $nproc_gmx
gmxMM   job_type = $job_type

gmxMM  GAUSSIAN INTERFACE
gmxMM   Layer: $layer
gmxMM   Gaussian input $input
gmxMM   Gaussian output $output
gmxMM   Gaussian message $msg

EOF
cat<<EOF >>gmx_${label}.log
INPUT SUMMARY:

GROMACS PARAMETERS
 mdpfile = $mdpfile
 mdpfile_nm = $mdpfile_nm
 mdpfile_min = $mdpfile_min
 topfile = $topfile
 struct file = $trrfile
 save_intermediate = $save_intermediate
 nproc_gmx = $nproc_gmx
 job_type = $job_type

GAUSSIAN INTERFACE
 Layer: $layer
 Gaussian input $input
 Gaussian output $output
 Gaussian message $msg

EOF

# Get data from G09 to update the trr file (stdout is natoms)
if [ "$aa2ua" == "NO" ]; then
    aa2ua_=""
else
    aa2ua_=$aa2ua
fi
# Set if we are giving a geom in gro format or directly use gaussian one
geomdir=$grofile # if none given, it will be set to use gaussian geom
natoms_gmx_vis=$(gau2gro $input $layer $trrfile $geomdir $aa2ua_)
if (( $? )); then 
    echo "Error in $0 while running gau2gro" >> $msg
    exit 1;
fi

# If we use the geom from gaussian, it has been generated in gau2groas geom.gro
if [ "$geomdir" == "gaussian" ]; then
    grofile=${label}.gro
fi
cat<<EOF >>$msg
gmxMM   Number of gmx atoms mapping gaussian atoms $natoms_gmx_vis

EOF
cat<<EOF >>gmx_${label}.log

Number of gmx atoms mapping gaussian atoms $natoms_gmx_vis

EOF


######################################################################################
#   G R O M A C S    R U N
######################################################################################
echo "=======================" >> gmx_${label}.log
echo "STARTING GROMACS JOBS" >> gmx_${label}.log
echo "=======================" >> gmx_${label}.log
echo " " >> gmx_${label}.log
# Microiteration (only on Real system)
if [ "$job_type" == "MicroOpt" -a "$layer" == "R" ]; then
    echo "--------------------------------------" >> gmx_${label}.log
    echo "Performing microiteration step" >> gmx_${label}.log
    echo "--------------------------------------" >> gmx_${label}.log
    echo " " >> gmx_${label}.log
    # Set model system as freeze group: name MODEL
    $gmxcall make_ndx$gmxsufix -f $grofile_M <<EOF
name 0 MODEL
q

EOF
    # Run minimzation step (TODO: -maxwarn 1 looks pretty unsafe...)
    $gmxcall grompp$gmxsufix -n -c $grofile -p $topfile -t $trrfile -f $mdpfile_min -maxwarn 1 -o topol_${label}.tpr &>>gmx_${label}.log
    if (( $? )); then
        echo "=================="
        echo "An ERROR occurred!"
        echo "=================="
        grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
        exit 1
    fi
    if [ "$nproc_gmx" == "1" ]; then
        $gmxcall mdrun$gmxsufix  -s topol_${label}.tpr -deffnm minim_${label} -o $trrfile &>>gmx_${label}.log 
        if (( $? )); then
            echo "=================="
            echo "An ERROR occurred!"
            echo "=================="
            grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
            exit 1
        fi
    else
        mpiexec -np $nproc_gmx $gmxcall mdrun$gmxsufix  -s topol_${label}.tpr -deffnm minim_${label} -o $trrfile &>>gmx_${label}.log
        if (( $? )); then
            echo "=================="
            echo "An ERROR occurred!"
            echo "=================="
            grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
            exit 1
        fi
    fi
fi

# Single point calculation should be 0-step MD with unconstrained_start=yes
# or 0-step energy minimization 
# -- see GROMACS Mailing List: http://lists.gromacs.org/pipermail/gmx-users/2011-March/059176.html
echo " " >> gmx_${label}.log
echo "------------------------------" >> gmx_${label}.log
echo "Computing energy and dipole" >> gmx_${label}.log
echo "------------------------------" >> gmx_${label}.log
echo " " >> gmx_${label}.log
$gmxcall grompp$gmxsufix -c $grofile -p $topfile -t $trrfile -f $mdpfile -maxwarn 1 -o topol_${label}.tpr 2>>gmx_${label}.log 1>>gmx_${label}.log
if (( $? )); then
    echo "=================="
    echo "An ERROR occurred!"
    echo "=================="
    grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
    exit 1
fi
if [ "$nproc_gmx" == "1" ]; then
    $gmxcall mdrun$gmxsufix -s topol_${label}.tpr -deffnm out_${label} -o $trrfile 2>>gmx_${label}.log 1>>gmx_${label}.log
    if (( $? )); then
        echo "=================="
        echo "An ERROR occurred!"
        echo "=================="
        grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
        exit 1
    fi
else
    mpiexec -np $nproc_gmx $gmxcall mdrun$gmxsufix -s topol_${label}.tpr -deffnm out_${label} -o $trrfile 2>>gmx_${label}.log 1>>gmx_${label}.log
    if (( $? )); then
        echo "=================="
        echo "An ERROR occurred!"
        echo "=================="
        grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
        exit 1
    fi
fi

# ANALYSYS PROGRAMS, TO GET REQUIRED DATA FOR G09:

# First we read from the gaussian input the information in the first line
read natoms Ideriv charge spin  <$input

# 1 Energy and Dipole moment from $gmxcall g_energy
xvg_off='-xvg none'

$gmxcall ${gmxprefix}energy$gmxsufix -f out_${label}.edr -o energy_and_mu_${label}.xvg $xvg_off <<EOF 2>>gmx_${label}.log 1>>gmx_${label}.log
Potential
Mu-X
Mu-Y
Mu-Z
0
EOF
# gmx energy does not compute the dipole in new versions. Using dipoles instead
echo "0" | $gmxcall ${gmxprefix}dipoles$gmxsufix -f $trrfile -s topol_${label}.tpr -xvg none -o Mtot_${label}.xvg 2>>gmx_${label}.log 1>>gmx_${label}.log
# And re-make the old style file
read null mux muy muz null <Mtot_${label}.xvg
read t E null <energy_and_mu_${label}.xvg
echo "$t $E $mux $muy $muz" > energy_and_mu_${label}.xvg
rm dipdist.xvg epsilon.xvg aver.xvg Mtot_${label}.xvg

read null En dx dy dz <energy_and_mu_${label}.xvg
cat <<EOF >>$msg
gmxMM  GROMACS Energy"
gmxMM     $En  (KJ/mol)" 
gmxMM  GROMACS Dipole" 
gmxMM     ( $dx , $dy , $dz )  (Debye)" 

EOF

# We also include the number of atoms in this file taken from [gauss input]-NO: is from update_trr_v2
#read natoms null <$input
echo $natoms >>energy_and_mu_${label}.xvg
echo $natoms_gmx_vis >>energy_and_mu_${label}.xvg

# 2. Forces (=-Gradients) from g_traj
# td m06/6-31G(d) td
if [ "$Ideriv" -ge "1" ]; then
    echo " " >> gmx_${label}.log
    echo "-----------------------" >> gmx_${label}.log
    echo "Computing gradient" >> gmx_${label}.log
    echo "-----------------------" >> gmx_${label}.log
    echo " " >> gmx_${label}.log
    xvg_off='-xvg none'

    gradfile="forces_${label}.xvg"
    $gmxcall ${dumprefix}dump$gmxsufix -f $trrfile 2>>gmx_${label}.log 1>>gmx_${label}.log
    $gmxcall ${dumprefix}dump$gmxsufix -f $trrfile &> forces_${label}_tmp; grep "f\[" forces_${label}_tmp > $gradfile; rm forces_${label}_tmp
else
    gradfile='NO'
fi

# 3. Polarizability and dipole derivatives ??
#    If not available set to 0
    polfile='ZERO'


# 4. Dipole derivatives
#    Get the dipole derivatives from the atomic charges (for non-polarizable FF):
#                           dmu/dxi = qi
ddipfile=ddip_${label}.dat
$gmxcall ${dumprefix}dump$gmxsufix -s topol_${label}.tpr 2> gmx_${label}.log | grep "atom\[" | egrep "q=[ \-]{1}[0-9.e+\-]+" -o | sed "s/q=//" > chr_${label}.dat 
ddip_nonpolar < chr_${label}.dat > $ddipfile


# 5. Force constants (Hessian) -- if frequencies are needed
if [ "$Ideriv" -eq "2" ]; then
    echo " " >> gmx_${label}.log 
    echo "-----------------------" >> gmx_${label}.log
    echo "Computing Hessian" >> gmx_${label}.log
    echo "-----------------------" >> gmx_${label}.log
    echo " " >> gmx_${label}.log
    echo "gmxMM The hessian is also built" >>$msg
    echo >>$msg
    hessfile="hessian_${label}.xvg"
    #     Run job with integrator=nm (normal modes)
    $gmxcall grompp$gmxsufix -c $grofile -p $topfile -t $trrfile -f $mdpfile_nm -maxwarn 1 -o topol_${label}.tpr 2>>gmx_${label}.log 1>>gmx_${label}.log
    if (( $? )); then
        echo "=================="
        echo "An ERROR occurred!"
        echo "=================="
        grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
        exit 1
    fi
    if [ "$nproc_gmx" == "1" ]; then
        $gmxcall mdrun$gmxsufix  -s topol_${label}.tpr -deffnm nm_${label} 2>>gmx_${label}.log 1>>gmx_${label}.log
        if (( $? )); then
            echo "=================="
            echo "An ERROR occurred!"
            echo "=================="
            grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
            exit 1
        fi
    else
        mpiexec -np $nproc_gmx $gmxcall mdrun$gmxsufix  -s topol_${label}.tpr -deffnm nm_${label}.mtx 2>>gmx_${label}.log 1>>gmx_${label}.log
        if (( $? )); then
            echo "=================="
            echo "An ERROR occurred!"
            echo "=================="
            grep "Fatal error" gmx_${label}.log -A 5 -B 5 >> $msg
            exit 1
        fi
    fi

    #     Turn the binary hessian (nm.mtx) into ascii format (dump?)
    $gmxcall ${dumprefix}dump$gmxsufix -mtx nm_${label}.mtx 1>$hessfile 2>/dev/null
    sed -i "1d" $hessfile

else
    hessfile='NO'
fi



########################################################################################
#  U P D A T E    G A U S S I A N
########################################################################################
# Write data into proper format for Gaussian
#       1)                         2)        3)       4)        5)        6)      7)   8)
gro2gau energy_and_mu_${label}.xvg $gradfile $polfile $ddipfile $hessfile $output $msg $aa2ua

# Remove or save intermediate MM files
if $save_intermediate; then
    echo "Keeping all intermediate files from gromacs run" >> $msg
else
    rm *${label}*
fi
# Always delete backups
rm *\# mdout.mdp
