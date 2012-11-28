/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

/*  #define _VERBOSE  */


#include <string.h>
#include <math.h>
#include "main.h"
#include "macros.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
//#include "fatal.h"
#include "xtcio.h"
#include "enxio.h"
#include "assert.h"
#include "smalloc.h"
#include "names.h"
#include "gmxfio.h"
#include "tpxio.h"
#include "trnio.h"
#include "txtdump.h"
#include "vec.h"
#include "statutil.h"
#include "gp_memory.h"
#include "index.h"
#include "confio.h"
#include "pbc.h"
#include "lap.h"
#include "gp_memory.h"
#include "gnrl.h"




/* pointers to two debug output files, only used in debug mode */
#ifdef __DEBUG
static FILE *assignment, *costmatrix;
#endif

/* these variables containing the cost matrix and the solution of the LAP. *
 * They would not neccessary have to be global */
cost **assigncost, *u, *v, lapcost;   
row *colsol;
col *rowsol;

static bool bCom=FALSE;
static bool bPBC=TRUE;
static bool bOrp=FALSE;
static bool bOref=FALSE;
static bool bRef=FALSE;


/* init_permute_stuff - allocate memory for the LAP variables */
void init_permute_stuff(int natoms) {
    assigncost = (cost **)alloc_costmatrix(natoms,natoms); 
    rowsol = (row *)alloc_ints(natoms); 
    colsol = (col *)alloc_ints(natoms); 
    u = (cost *)alloc_costs(natoms); 
    v = (cost *)alloc_costs(natoms);  
    
    if(assigncost == NULL || rowsol == NULL || colsol == NULL
       || u == NULL || v == NULL){
/*	free_costmatrix(assigncost, natoms, natoms);not necessary, will be directly freed, when malloc failed*/
	free(rowsol);
	free(colsol);
	free(u);
	free(v);
	
	printf("cannot allocate memory\n");
	exit(0);
    }   
}

/* init_permute_stuff - free the memory of the LAP variables */
void free_permute_stuff(int dim) {
    free_costmatrix(assigncost,dim,dim);
    free_ints(rowsol);
    free_ints(colsol);
    free_costs(u);
    free_costs(v);
};    

/* permute - given a reference structure 'ref' *
 *           permute a structure 'now' containing 'natoms' *
 *           using the atoms with ids 'oxygens' to build up 
 *           the costmatrix. Its 'result' contains the solution of the LAP *
 *           (i.e., the ids of the relabeled atoms) *
 * *
 *           'permute' is basically a wrapper for the 'lap' function *           
 *           it is called by 'permute_trx' which performs all *
 *           administrative tasks s.a. reading / writing the files *
 *           and relabeling the atoms
 */          

void permute(rvec *ref,rvec *now,int natoms,atom_id* result,atom_id *oxygens,t_pbc * pbc) {
    int i, j;
    
    cost oldconf;  /* costs of the old (unpermuted) and new (permuted) */
    cost permconf; /* configurations */
    rvec dx;
    
    /* compute the cost matrix */
    for (i = 0; i < natoms; i++) {
	for (j = 0; j < natoms; j++) {
	    if (bPBC){
		pbc_dx(pbc,ref[oxygens[i]],now[oxygens[j]], dx);
		assigncost[i][j] = norm2(dx);
	    }
	    else
		assigncost[i][j] = distance2(ref[oxygens[i]],now[oxygens[j]]);
	}
#ifdef __DEBUG	
	fwrite(assigncost[i], sizeof(cost), natoms, costmatrix);
#endif
    }; 
    
#ifdef __DEBUG	
    fflush(costmatrix);
#endif
	
    /*   start with clean vectors */
    for (i=0;i<natoms;i++) { 
	result[i]=0;
	colsol[i]=0;
	u[i]=0;
	v[i]=0;
    };
    
    /* now solve the linear assignment problem to obtain the *
     * optimal permutation for this frame */
    lapcost = lap(natoms, assigncost, result, colsol, u, v);
    
    /* if we are debugging, dump the assignment to a file */
#ifdef __DEBUG
    for(i=0; i<natoms; i++){
	fprintf(assignment, "%d ", result[i]);
    }
    fprintf(assignment, "\n");
    fflush(assignment);
#endif

    /* check the computed assignment for consistency */
    checklap(natoms, assigncost, result, colsol, u, v);     
    
    /*   just to be sure, check whether we obtained a better configuration by *
     *   permuting the atoms */
    oldconf = permconf = 0; 
    for(i=0; i<natoms; i++){ 
	oldconf  += distance2(ref[oxygens[i]], now[oxygens[i]]); 
	permconf += distance2(ref[oxygens[i]], now[oxygens[result[i]]]); 
    } 
    
    if(oldconf < permconf)
	printf("Warning: permutation does not reduce the distance\n");
	
#ifdef __DEBUG    
    /* allow debug file to be overwritten again */
    rewind(costmatrix); 
#endif

    return;
};

/*builds center or mass coordinates if no reference structure was given */
void  BuildxComAr(int lnRes,t_topology ltop,t_trxframe lframe,int *lnAtomsRes,atom_id **lresIndex,rvec *lxCom){
 int ix,i,j;
 rvec x_tmp2;
 real M,mass;
      /* Get center of mass of residues */
   for (i=0; i<lnRes; i++){
		M=0;
		clear_rvec(lxCom[i]);
		for (j=0; j<lnAtomsRes[i]; j++){
		    ix = lresIndex[i][j];
		    mass  = ltop.atoms.atom[ix].m;
		    M += mass;
		    svmul(mass,lframe.x[ix],x_tmp2);
		    rvec_inc(lxCom[i],x_tmp2);
		}
      	for(j=0; (j<DIM); j++){
				    lxCom[i][j] /= M;
		}
   }
}

/*initialise frame */
/*PURPOSE: when in COM mode and no reference struct is supplied*/
t_trxframe initFrame(int step,real time,t_atoms *atoms,rvec *vec,matrix box,int natoms) {
    t_trxframe comframe;
    
    clear_trxframe(&comframe,TRUE);
    comframe.flags = 1;
    comframe.natoms = natoms;
    comframe.bStep = TRUE;
    comframe.step = step;
    comframe.bTime = TRUE;
    comframe.time = time;
    comframe.bAtoms = atoms!=NULL;
    comframe.atoms = atoms;
    comframe.bX = TRUE;
    comframe.x = vec;
    comframe.bV = FALSE;
    comframe.v = NULL;
    comframe.bBox = TRUE;
    copy_mat(box,comframe.box);
    
    return comframe;
}

/* Copy trxframe frame and box info from source to dest */
void copy_trxframe_stuff(t_trxframe * src, t_trxframe * dest)
{
    dest->time = src->time;
    copy_mat(src->box, dest->box);
    dest->t0 = src->t0;
    dest->tpf = src->tpf;
    dest->tppf = src->tppf;
    dest->step = src->step;
    dest->prec = src->prec;
}



/*print frame info up to no atoms of specified frame */
/*PURPOSE: debugging  - mainly*/
void printFrame(t_trxframe comframe,int no){
    int i,j;

    printf("printing frame:\n step %d time %f\n",comframe.step, comframe.time);
    printf(" natoms %d flags %d\n",comframe.natoms, comframe.flags);
    
    for (i=0;i<no;i++){
            printf("\ni %d frame %f %f %f\n",i,comframe.x[i][0],comframe.x[i][1],comframe.x[i][2]);fflush(stdout);
            
    }
    printf("\n printing box vectors\n");

        for (j=0;j<3;j++){
            printf("XX  %f",comframe.box[XX][j]);
            printf("YY  %f",comframe.box[YY][j]);
            printf("ZZ  %f",comframe.box[ZZ][j]);
        }
        printf("\n");
}


/* permute_trx - given the filenames of 'infile', 'indexfile'                       *
 *               and 'outfile' and 'reffile' (if any) and the 'molsize'             *
 *               of the molecules to be relabeled, read all the files               *
 *               and solve the LAP for each frame of 'infile'.                      *
 *               The permuted trajectory is written to 'outfile'                    *
 *               -orp on command line triggers output of CenterOfMass coords to              *
 *               refperm.xtc and -oref will write COM coords of original trajectory *
 *               to ref.xtc                                                         */

void permute_trx(char *infile, char *indexfile, char* outfile, char *reffile, char *resid_file,
				 char *outfile_ref, char *out_refperm_file, int molsize)
{
    int ftp, ftpndx; /* flags holding the file types of infile and indexfile*/
    int trxhandle;   /* file handle for infile*/
    int trxout;      /* file handle for outfile*/

    int natoms;
    int i,j, res;    /*loop variables*/
    int		nRes=0;
    t_topology  top;
	int *nAtomsRes=NULL;
	atom_id **resIndex=NULL;
    int ndist, nref, nsol; /* number of atoms for distance matrix *
			    * calculation, reference structure and solvent group */
    bOref = outfile_ref == NULL ? 0 : 1;
    bOrp = out_refperm_file == NULL ? 0 : 1;
    t_trxframe frame;
    rvec *refstruct, *xref;

    atom_id *dist, *sol;   /* atom ids for distance calculation and solvent group*/
    atom_id *result=NULL;       /* atom ids returned from the LAP solver (i.e., only ndist atoms)*/
    atom_id *all;          /* atom ids of the whole relabeled box (i.e, natoms atoms ) */

    char *distname, *solname, str[STRLEN];

    t_pbc pbc;
    char title[STRLEN];
    int statusCom,statusComPerm,resnr,lastRes=-100, nThisRes=0, iRes, ix, nAtomsPerResMax=0;
    int trxout_ref=0, trxout_refPerm=0, d,m;
    rvec *x_tmp=NULL;
    matrix box;
    atom_id *zeroToNResArray=NULL;
    rvec *xCom=NULL,x_tmp2, *xcomRes, hbox;
    t_trxframe comFrame, comFramePerm;
    real tCom;
    rvec *v;
    t_atoms *atoms,*comAtoms;
	
    int ePBC = 0;

    ftp = fn2ftp(infile);
    ftpndx = fn2ftp(indexfile);
    
    /* let the user choose the distance and the solvent group */
    /* see g_permute user documentation for details */
    printf("\nChoose a group for the distance calculation:\n"); 
    rd_index(indexfile,1,&ndist,&dist,&distname);   
    printf("\nChoose the solvent group:\n"); 
    rd_index(indexfile,1,&nsol,&sol,&solname);   
    
#ifdef __VERBOSE
    printf("read index: ndist = %d\n", ndist);
    for(i=0; i<ndist; i++){
	printf("%d ", dist[i]);
    }
    printf ("\n");
#endif
    
    /* open the output file */
    printf ("writing trajectory to %s\n", outfile);
    trxout=open_trx(outfile,"w"); 
    
    

    /* initialize the vector containing the permuted positions */
    
    /* open the trajectory and read the first frame */
    if (!read_first_frame(&trxhandle,infile,&frame,TRX_READ_X || TRX_READ_V || TRX_READ_F))
	gmx_fatal(FARGS,"Could not read first frame from trajectory %s",infile);
    
    natoms=frame.natoms;

    snew(atoms,natoms);
    snew(v, natoms);
    snew(refstruct,natoms);
    
	   /*   allocate space for the reference structure */

	if (!bCom) {
      snew(result,ndist);
		/* allocate memory for the LAP solver */
   	  init_permute_stuff(ndist);
    }
    
     
    /* allocate the vector containing the indices of all atoms */
    snew(all,natoms);
    /*read in tpr if (bCom)*/   
  
    if (bCom){
        read_tps_conf(resid_file, title,&top, &ePBC, &x_tmp, NULL , box, TRUE); /*must have*/
	    for (i=0; i<ndist; i++)
	        if ( (resnr=top.atoms.atom[dist[i]].resnr) != lastRes){
		        lastRes=resnr;
	      	    nRes++;
		        if (nThisRes>nAtomsPerResMax)
		            nAtomsPerResMax=nThisRes;
		            nThisRes=1;
	        }
	        else{
		        nThisRes++;
            }
        if (!bRef) {
            snew(refstruct,nRes);
        }
        snew(comAtoms,nRes);

        init_t_atoms(comAtoms,nRes,FALSE);
        snew(nAtomsRes,nRes);
        snew(xCom,nRes);
        snew(zeroToNResArray,nRes);
        snew(resIndex,nRes);
        for (i=0; i<nRes; i++){
            snew(resIndex[i],nAtomsPerResMax);
            zeroToNResArray[i]=i;
        }
        lastRes=-100;
        iRes=-1;nAtomsRes[iRes]=0;
	    j=0;
	    for (i=0; i<ndist; i++){
            if ( (resnr=top.atoms.atom[dist[i]].resnr) != lastRes){
                if (iRes >= 0){
                    nAtomsRes[iRes]=j+1;
                }
                iRes++;
                lastRes=resnr;
                j=0;
            }
            else{
                j++;
            }
            resIndex[iRes][j]=dist[i];
        }
        j++;
        nAtomsRes[nRes-1]=j;
    } /*if bCom end*/


    /*     if possible, load the reference structure from a data file *
     *     if no -r is given use structure from first frame           *
     *     and fill it into refstruct                                 */

	if(bRef){
		/* create some dummy variables to store excess information */

		get_stx_coordnum(reffile, &(atoms->nr));
		nref = atoms->nr;
		init_t_atoms(atoms, atoms->nr, TRUE);
		if (nref != natoms && bCom == FALSE){
		  gmx_fatal(FARGS,
				"found %d atoms in the trajectory frames, "
				"but %d atoms in reference structure", natoms, nref);
		} else if (nref != nRes && bCom == TRUE) {
            /*bCom on and different no of COM coords in reffile then residues in structure file*/
              gmx_fatal(FARGS,
                "found %d residues in the structure file, "
                "but %d COM coordinates in reference structure", nRes, nref);
        }
		read_stx_conf(reffile, title, atoms, refstruct, v, &ePBC, box);
		sfree(v); sfree(atoms);
	} else {
	/*     the user did not supply a reference structure */
	/*     use the first frame instead */
       if (!bCom){
           for (i=0;i<natoms;i++) {
           copy_rvec(frame.x[i],refstruct[i]);
           }
       } else {
            printf ("\nNo reference structure provided, building COM coordinates from first frame\n");
		/*we built the com coords for each residue and save it in refstruct*/
/*	in frame are already the coords we need*/	
           BuildxComAr(nRes, top, frame,nAtomsRes,resIndex,refstruct);

	   }
    } /*end of if reffile*/
    


    if (bCom){
		if (bRef){
			if (bOref)
     	      read_first_frame(&statusCom,     reffile, &comFrame,     0);
  		    if (bOrp)
   	            read_first_frame(&statusComPerm, reffile, &comFramePerm, 0);
        } else {
             if (bOref) {
     	        comFrame=initFrame(0,0,comAtoms,xCom,box,nRes);
             }
  		     if (bOrp) {
        	     comFramePerm=initFrame(0,0,comAtoms,xCom,box,nRes);
             }
		}

        snew(result,nRes);
        init_permute_stuff(nRes);

        if (bOref) {
        printf ("\nwriting trajectory of centers of mass to %s\n", outfile_ref);
        trxout_ref=open_trx(outfile_ref,"w");
        }
        if (bOrp) {
            printf ("\nwriting trajectory of permuted centers of mass to %s\n", out_refperm_file);
            trxout_refPerm=open_trx(out_refperm_file,"w");
        }
    } else {
        if ((bOref)){
            printf("\nOutput to %s only when -com is specified\n\n",outfile_ref);
        }
        if ((bOrp)){
            printf("\nOutput to %s only when -com is specified\n\n",out_refperm_file);
        }
    }/*end of if bCom*/

    /*read all the following structures permute them and dump them in the *
     * output file *
     * Note that the first structure is skipped !*/
 do {

	set_pbc(&pbc, ePBC, frame.box);
	fprintf(stderr,"Time %6.3f ps\n",frame.time);
	/* call the LAP solver */
	if (!bCom)
	    permute(refstruct,frame.x,ndist,result,dist,&pbc);
	else{
		BuildxComAr(nRes, top, frame,nAtomsRes,resIndex,xCom);
		if (bOref) {
    /* Write the positions of the center of mass of each residue */
            for (i=0; i<nRes; i++){
                copy_rvec(xCom[i],comFrame.x[i]);
            }

            copy_trxframe_stuff(&frame, &comFrame);
            write_trxframe(trxout_ref,&comFrame);
        }
	    /* Now, permute uses all positions in xCom, therefore
	       zeroToNResArray contains only 0,1,2,...,nRes-1. */
	    permute(refstruct,xCom,nRes,result,zeroToNResArray,&pbc);
	}/*end of if (!bCom)*/
        /* now relabel the atoms in the current frame according to  */
        /* the optimal permutation 'result' */
	
	/* start from the identity permutation */ 
	for(i=0; i<natoms; i++){
	    all[i] = i;
	}
	/*flip the positions of all solvent molecules according to 'result'*
	 * the for loop loops through all solvent atoms  *
	 * j stores the current molecules' index *
	 * res loops through the atoms in a molecule *
	 * so that each solvent molecules is relabeled as a whole */
	for (i=0;i<nsol;i++) {
	    res = i / molsize;
	    j = i % molsize;
	    all[sol[i]]=sol[molsize*result[res]+j];
	}
	/* If the COM is more that half a boxlength away from the reference position
	   shift the whole molecule (and its center of mass)                         */
	if (bCom){
	    for(d=0; d<DIM; d++)
		hbox[d] = 0.5*frame.box[d][d];
	    for (i=0;i<nRes;i++){
		for(m=DIM-1; m>=0; m--)
		    if (hbox[m] > 0){
			while (xCom[result[i]][m]-refstruct[i][m] <= -hbox[m])
			    for(d=0; d<=m; d++){
				for (j=0; j<molsize; j++)
				    frame.x[molsize*result[i]+j][d] += frame.box[m][d];
				    xCom[result[i]][d]              += frame.box[m][d];
			    }
			while (xCom[result[i]][m]-refstruct[i][m] > hbox[m])
			    for(d=0; d<=m; d++){
				for (j=0; j<molsize; j++)
				    frame.x[molsize*result[i]+j][d] -= frame.box[m][d];	
				    xCom[result[i]][d]              -= frame.box[m][d];
			    }
		    }
	    }
	}
	/*we are done: write the relabeled version of the current frame*/
	write_trxframe_indexed(trxout,&frame,natoms,all);
	
	if ((bOrp) && (bCom)){
        for (i=0;i<nRes;i++) {
		    copy_rvec(xCom[result[i]], comFramePerm.x[i]);
        }
        copy_trxframe_stuff(&frame, &comFramePerm);
	    write_trxframe(trxout_refPerm,&comFramePerm);
    }
    }  while (read_next_frame(trxhandle,&frame));
    close_trx(trxout);
    if (!bCom)
        free_permute_stuff(ndist);
    else{
        free_permute_stuff(nRes);
        if (bOrp)
        close_trx(trxout_refPerm);
        if (bOref)
        close_trx(trxout_ref);
/*			for (i=0;i<nRes;i++)
          sfree(refstruct[i]);*/
 }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_permute reads ",
    "a trajectory ([TT].trj[tt]/[TT].trr[tt]/[TT].xtc[tt]) ",
    "output in a readable format and fits all frames against the ",
    "structure file using combinatorial optimization.[PAR]"
  };
  t_filenm fnm[] = {
      { efTRX, "-f", NULL, ffREAD },
      { efTRX, "-o", "permute.xtc",  ffWRITE},
      { efTRX, "-oref", "ref.xtc",  ffOPTWR},
      { efTRX, "-orp", "refperm.xtc",  ffOPTWR},
      { efNDX, "-n", "index.ndx",    ffREAD},
      { efSTX, "-r", "ref.gro",      ffOPTRD  },
      { efSTX, "-s", "frame.pdb",    ffOPTRD  }
  };

#define NFILE asize(fnm)

  /* Command line options */
  static int molsize = 3;
  char *out_file;
  char *in_file;
  char *index_file;
  char *ref_file;
  char *out_refperm_file;
  char *resid_file;
  char *out_file_ref;
  int ftp,ftpin, ftpndx;
  bool bOref,bOrp;
  t_pargs pa[] = {
      { "-m", FALSE, etINT, {&molsize}, "number of atoms of a solvent molecule" },  
      { "-com", FALSE, etBOOL, {&bCom}, "consider center of mass of each residue as reference point for distance calculation."
	"if desired give appropriate reference file with -r, otherwise reference structure is built from COM of residues" },
      { "-rm_pbc", FALSE, etBOOL, {&bPBC}, "calc distance using pbc, default is to remove pbc" },
  };  
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_BE_NICE ,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  /* open files for debug output (mainly crash forensics)*
  *  assignment - contains the latest assignment *
  *  costmatrix - contains the last cost matrix  */
#ifdef __DEBUG
  assignment = fopen("assign.dat", "w");
  costmatrix = fopen("matrix.dat", "w");
#endif /* _DEBUG*/

  in_file=         opt2fn("-f",NFILE,fnm);
  index_file=      opt2fn("-n",NFILE,fnm);
  out_file=        opt2fn("-o",NFILE,fnm);
  out_file_ref=    opt2fn("-oref",NFILE,fnm);
  ref_file=        opt2fn("-r",NFILE,fnm);
  out_refperm_file=opt2fn("-orp",NFILE,fnm);
  resid_file=      opt2fn("-s",NFILE,fnm);
  
 
  bOref=opt2bSet("-oref",NFILE,fnm);
  bOrp =opt2bSet("-orp",NFILE,fnm);
  bRef =opt2bSet("-r",NFILE,fnm);


  if (!fexist(ref_file) && (bRef)){
        gmx_fatal(FARGS,"Reference structure %s not found",ref_file);
  } else if (!bRef) {
        ref_file = NULL;
  }

  out_file_ref = bOref == 1 ? out_file_ref:NULL; /* making sure NULL is set if none of     *
                                                  * the two options -oref / -orp are given */
  out_refperm_file = bOrp == 1 ? out_refperm_file:NULL;

   /* Determine output type */ 
  ftp=fn2ftp(out_file);
  fprintf(stderr,"Will write %s: %s\n",ftp2ext(ftp),ftp2desc(ftp));

if (bCom && (!(resid_file)))
      gmx_fatal(FARGS,"\n\nYou must supply  a structure (-s)\nThe structure file (-s) is used to create the molecule groups using the residue id.\nYou also can provide a reference file, which must contain one atom per molecule you want to permute. Using the first structure is the default if -r is not provided a.\n\n");


  /* all IO has been set up: now we do your job and permute the trajectory */
  if (ftp2bSet(efTRX,NFILE,fnm)) 
      permute_trx(in_file,index_file, out_file, ref_file, resid_file, out_file_ref, out_refperm_file, molsize);
    
  thanx(stderr);

  return 0;
}
