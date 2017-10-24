    /*
    ! J Cerezo, May 2011
    ! ======================================================
    ! This program is part of GAUSSIAN -- GROMACS interface
    ! truough External="gmxMM4gauss.sh parameters.dat" keyword
    ! within Gaussian input file
    !=======================================================
    !
    ! Description
    ! -----------
    ! Read coordinates from Gaussian (optimization) step and
    ! update the corresponding trr file to feed Gromacs using
    ! xdr library
    !
    ! Version info
    ! ------------
    ! v3: Read gosht atoms from g09, using the aa2ua.dat file (0=ghost atom)
    !     This is need to let all MM part to Gromacs while using ONIOM type calulation.  
    ! v4: Added trrfile name as input 
    !
    ! Compilation instructions
    ! -------------------------
    !_ make-o$ LIBFOLDER /home/cerezo/Programas/xdrlib/lib/
    !_ make-o$ LIBNAME xdrfile
    !_ make-o$ INCLUDEFOLDER /home/cerezo/mis_bin/src/include
    !_ make$ gcc update_trr_v2.c -o update_trr_v2
    !make$ gcc xdrfile.c xdrfile_trr.c update_trr_v4.c -o update_trr_v4
    */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <time.h>
#include <float.h>
#include "xdrfile.h"
#include "xdrfile_trr.h" 

/* FACTORS */
#define BOHR_TO_NM 5.29177E-2
#define MAX_ATOMS 1000000

static void _die(char *msg, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die(msg) _die(msg,__LINE__,__FILE__)

static void _die_r(char *msg, int result,int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr,"result = %d\n",result);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die_r(msg,res) _die_r(msg,res,__LINE__,__FILE__)

char *remove_ext(char* mystr) {
    /* 
     * Function to remove extension from a filename
     * From: https://stackoverflow.com/questions/2736753/how-to-remove-extension-from-file-name
     */
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}


int main(int argc, char *argv[]) 
{
    char *gauss_input; /*="gaus_input.dat";*/
    char *condensed_atoms;
    char *trr_file;
    char *layer;
    char *geomdir;
    char *grofile;
    int nua, naa, nghost,ighost;
    int aa2ua;
    FILE *gin, *uain, *grout;
    XDRFILE *xd;
    int result,i,ii,j,k,nframes=1;
    int natoms, natoms_gmx, Ideriv;
    int *Z;
    int step=0;
    float time=1.0;
    float charge, spin, q;
    matrix box;
    rvec *xgauss;
    rvec *xgmx;
    float lambda=0.0;

    /* Read data from command line */
    if(argc<5) {
        printf("ERROR: Required at least %d arguments, found %d\n", 4, argc-1);
        printf(" Uso: %s <gaussian fnm> <layer> <trr fnm> <geomdir> <condensed_atoms fnm (opt)>\n", argv[0]);
        return 1;
     }  
     gauss_input=argv[1];
     layer=argv[2];
     trr_file=argv[3];
     geomdir=argv[4];

         
    /* Open g09 files */
    gin = fopen(gauss_input,"r");
    /*fopen(gauss_msg,"w"); */

    
    /* Set box vectors: 
     * Let them = 0 for the moment ("in vacuo" needs no pbc)
     */
    /*
    It's read from previous trr file
    for(i=0; (i<DIM); i++) {
		for(j=0; (j<DIM); j++) {
			box1[i][j] = 0.0;
        }
    }
    */

    /* =================== * 
     *  Read old trr file  * 
     * =================== */
    xgmx = calloc(MAX_ATOMS,DIM);
    if (strcmp("gaussian", geomdir) != 0) {
        
        xd = xdrfile_open(trr_file,"r");
        if (NULL == xd)
            die("Opening trr file for reading");
        
        /* Read xout */
        result = read_trr_natoms(trr_file,&natoms_gmx);
        if (0 != result)
            die_r("Reading trr file",result);
        result = read_trr(xd,natoms_gmx,&step,&time,&lambda,box,xgmx,NULL,NULL);
        
        xdrfile_close(xd);
    } else {
        /* Set box to 0.0 */
        for(i=0; (i<DIM); i++) {
        for(j=0; (j<DIM); j++) {
            box[i][j] = 0.0;
        }
        }
    } 
  

    /* =============================== * 
     *  Read gaussian input geommetry  * 
     * =============================== */
    /* Header information */
    fscanf(gin,"%d%d%f%f",&natoms, &Ideriv, &charge, &spin); 
    
    /* Store "raw" g09 coordinates in xgauss */
    xgauss = calloc(MAX_ATOMS,DIM);
    Z      = calloc(MAX_ATOMS,sizeof(int));
    if (NULL == xgauss)
        die("Allocating memory for xgauss before reading g09 input geom");

    for(i=0; (i<natoms); i++) {
        fscanf(gin,"%d", &Z[i]);
        for(j=0; (j<DIM); j++) {
            fscanf(gin,"%f", &xgauss[i][j]);
            xgauss[i][j] *= BOHR_TO_NM;
        }
        fscanf(gin,"%f%*[^\n]",&q);

    }

    fclose(gin);
   
    /* Trasform g09 coordinates to UA coordinates if required */  
    nua = 0;
    naa = 0;
    ii  = 0;
    if(argc>5) {
        condensed_atoms=argv[4];  
        uain = fopen(condensed_atoms,"r");    
        fscanf(uain,"%d%d%d",&naa, &nua, &nghost);
 
        for(i=0; (i<nua); i++) {
            fscanf(uain,"%d", &aa2ua);
            if ( aa2ua != 0) {
                for(j=0; (j<DIM); j++) {
                    xgmx[i][j] = xgauss[ii][j];
                }
                ii += aa2ua;
            } else {
                ii++;
            }
            /*
            That should be COM calculation...
            for(j=i; j<aa2ua[i]; j++) {
                x_ua[i] += 
            */
        }
        fclose(uain);
    }
    
    if ( natoms == (naa + nghost) ) {
        natoms = nua;
    } else {
        natoms += nua - naa;
    } 
    printf("%d\n",natoms);
    
    for(i=nua; (i<natoms); i++) {
        for(j=0; (j<DIM); j++) {
            xgmx[i][j] = xgauss[ii][j];
        }
        ii++;
    }
    
    if (strcmp("gaussian", geomdir) != 0) {
        /* If gaussian natoms is greater than gmx natoms, something went wrong... */
        if ( natoms_gmx < natoms ) {
            die("Less atoms in Gromacs than in Gaussian");
        }
        /* Only incorporate gmx gosht environment (extended trr file) if layer is R */
        if ( strcmp("S", layer) == 0 || strcmp("M", layer) == 0 ) {
            natoms_gmx = natoms;
        }
    } else {
        natoms_gmx = natoms;
        /* Write a gro file */
        grofile=remove_ext(gauss_input);
        strcat(grofile,"_tmp.gro");
        grout = fopen(grofile,"w");
        fprintf(grout,"Title\n");
        fprintf(grout,"%5i\n",natoms_gmx);
        for (i=0; (i<natoms_gmx); i++) {
            fprintf(grout,"%5i%5s%-5i%5i",1,"UNK  ",Z[i],i+1);
            for (j=0; (j<DIM); j++) {
                fprintf(grout,"%8.3f",xgauss[i][j]);
            }
            fprintf(grout,"\n");
        }
        fprintf(grout,"%10.5f%10.5f%10.5f\n",0.,0.,0.);
        fclose(grout);
    }

    /* ================ * 
     *  Write trr file  * 
     * ================ */
    /* Open trr to write on */
	xd = xdrfile_open(trr_file,"w");
	if (NULL == xd)
		die("Opening trr file for writing");
    
	/* Write the only frame */
		result = write_trr(xd,natoms_gmx,step,time,lambda,box,xgmx,NULL,NULL);
		if (0 != result)
			die_r("Writing trr file",result);

	xdrfile_close(xd);
    
    return 0;
	
}

