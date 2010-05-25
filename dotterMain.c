/*  File dotterMain.c
 *  Author: esr
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: 
 * Exported functions: 
 * HISTORY:
 * Last edited: Aug 26 15:42 2009 (edgrif)
 * Created: Thu Aug 26 17:17:30 1999 (fw)
 * CVS info:   $Id: dotterMain.c,v 1.6 2010-05-25 15:03:50 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <ctype.h>
#include <wh/regular.h>
#include <wh/graph.h>
#include <wh/gex.h>
#include <SeqTools/blixem_.h>
#include <wh/dotter_.h>

#define UNSET_INT  -1


typedef struct _DotterOptions
  {
    int qoffset;
    int soffset; 
    int selfcall;
    int qlen;
    int slen;
    int revcompq;
    int dotterZoom;
    int install : 1;
    int pixelFacset;
    int seqInSFS;
    
    float memoryLimit;
    
    char *savefile;
    char *loadfile;
    char *FSfilename;
    char *mtxfile;
    
    char *winsize;
    
    char *qname;
    char *sname;
  } DotterOptions;



static void strNamecpy(char *dest, char *src)
{
  char *cp;

  while (*src && *src == ' ')
    src++;

  strcpy(dest, src);

  if ((cp = (char *)strchr(dest, ' ')))
    *cp = 0;
  if ((cp = (char *)strchr(dest, '\n')))
    *cp = 0;

  return ;
}

static char *stringUnprotect(char **textp, char *target)
{
  char *cp, *cpd;
  int count = 0;

 redo:
  cp = *textp;
  cpd = target;
  
  while (*cp)
    {
      if (*cp == '"')
	{
	  cp++;						    /* skip quote */
	  break ;
	}
      else
	cp++ ;
    }

  while (*cp != '"' && *cp)
    {
      if (*cp == '$')
	cp++;

      if (cpd)
	*cpd++ = *cp;
      else
	count++;

      cp++;
    }
  
  if (!target)
    {
      target = messalloc(count+1);
      goto redo;
    }
  
  *cp = '\0' ;
  *textp = cp+1; /* skip quote */

  return target;
}

  
static void addBreakline (MSP **MSPlist, char *name, char *desc, int pos, char seq) {

   MSP   
       *msp;
   char 
       *cp;

   if (!*MSPlist) {
       *MSPlist = (MSP *)messalloc(sizeof(MSP));
       msp = *MSPlist;
   }
   else {
     msp = *MSPlist;
     while(msp->next) msp = msp->next;
     msp->next = (MSP *)messalloc(sizeof(MSP));
     msp = msp->next;
   }

   msp->qname = messalloc(strlen(name)+1);
   strcpy(msp->qname, name);

   msp->desc = messalloc(strlen(desc)+1);
   strcpy(msp->desc, desc);
   if ((cp = (char *)strchr(msp->desc, ' ')))
     *cp = 0;
   if ((cp = (char *)strchr(msp->desc, '\n')))
     *cp = 0;

   msp->qstart = msp->qend = pos;
   *msp->sframe = seq;
   msp->fsColor = DARKGREEN;
   msp->type = BLXMSP_FS_SEG;
   msp->score = 100;
   insertFS(msp, "chain_separator");
}		      


static char* getUsageText()
{

  char *usage;
  static char usageText[] = "\
\n\
 Dotter - Sequence dotplots with image enhancement tools.\n\
\n\
 Reference: Sonnhammer ELL & Durbin R (1995). A dot-matrix program\n\
 with dynamic threshold control suited for genomic DNA and protein\n\
 sequence analysis. Gene 167(2):GC1-10.\n\
\n\
 Usage: dotter [options] <horizontal_sequence> <vertical_sequence>  [X options]\n\
\n\
 Allowed types:                Protein        -      Protein\n\
                               DNA            -      DNA\n\
                               DNA            -      Protein\n\
\n\
 Options:\n\
\n\
 -b <file>      Batch mode, write dotplot to <file>\n\
 -l <file>      Load dotplot from <file>\n\
 -m <float>     Memory usage limit in Mb (default 0.5)\n\
 -z <int>       Set zoom (compression) factor\n\
 -p <int>       Set pixel factor manually (ratio pixelvalue/score)\n\
 -W <int>       Set sliding window size. (K => Karlin/Altschul estimate)\n\
 -M <file>      Read in score matrix from <file> (Blast format; Default: Blosum62).\n\
 -F <file>      Read in sequences and data from <file> (replaces sequencefiles).\n\
 -f <file>      Read feature segments from <file>\n\
 -H             Do not calculate dotplot at startup.\n\
 -R             Reversed Greyramp tool at start.\n\
 -r             Reverse and complement horizontal_sequence (DNA vs Protein)\n\
 -D             Don't display mirror image in self comparisons\n\
 -w             For DNA: horizontal_sequence top strand only (Watson)\n\
 -c             For DNA: horizontal_sequence bottom strand only (Crick)\n\
 -q <int>       Horizontal_sequence offset\n\
 -s <int>       Vertical_sequence offset\n\
\n\
 Some X options:\n\
 -acefont <font> Main font.\n\
 -font    <font> Menu font.\n\
\n\
 See http://www.cgb.ki.se/cgb/groups/sonnhammer/Dotter.html for more info.\n\
\n\
 by Erik.Sonnhammer@cgb.ki.se\n\
 Version ";
    
  extern char *dotterVersion;

  static char *cc_date = 
#if defined(__DATE__)
  __DATE__
#else
  ""
#endif
  ;
  
  usage = messalloc(strlen(usageText) + strlen(dotterVersion) + strlen(cc_date) + 20);
  sprintf(usage, "%s%s, compiled %s\n", usageText, dotterVersion, cc_date);
  
  return usage;
}


int main(int argc, char **argv)
{
  DotterOptions options = {0, 0, 0, UNSET_INT, UNSET_INT, 0, 0, 1, 0, 0, 0.0, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  
  char   
      *qseq=0, *sseq=0, line[MAXLINE+1],
      *cp, *cc, *cq, type, 
      *firstdesc, *qfilename, *sfilename,
      text[MAXLINE+1];

  char opts[] = "D    ";    /* D  display mirror image
                               W  only watson strand
                               C  only crick strand
                               H  only HSPs
                               S  start with swapped greyramptool
                            */
  FILE *qfile, *sfile;
  MSP *MSPlist=0;
  MSP *msp;

  
  int          optc;
  extern int   optind;
  extern char *optarg;
  char        *optstring="b:cDf:F:Hil:M:m:p:q:Rrs:SW:wz:";

  extern char *dotterBinary;

  char *usage = getUsageText();

    
  while ((optc = getopt(argc, argv, optstring)) != EOF)
    {
      switch (optc) 
        {
          case 'b': 
            options.savefile = messalloc(strlen(optarg)+1);
            strcpy(options.savefile, optarg);         break;
          case 'c': opts[1] = 'C';              break;
          case 'D': opts[0] = ' ';              break;
          case 'f': 
            options.FSfilename = messalloc(strlen(optarg)+1);
            strcpy(options.FSfilename, optarg);       break;
          case 'F': 
            options.seqInSFS = 1;        
            options.FSfilename = messalloc(strlen(optarg)+1);
            strcpy(options.FSfilename, optarg);       break;
          case 'H': opts[2] = 'H';              break;
          case 'i': options.install = 0;                break;
          case 'l': 
            options.loadfile = messalloc(strlen(optarg)+1);
            strcpy(options.loadfile, optarg);         break;
          case 'M': 
            options.mtxfile = messalloc(strlen(optarg)+1);
            strcpy(options.mtxfile, optarg);          break;
          case 'm': options.memoryLimit = atof(optarg); break;
          case 'p': options.pixelFacset = atoi(optarg); break;
          case 'q': options.qoffset = atoi(optarg);     break;
          case 'R': opts[3] = 'S';              break;
          case 'r': options.revcompq = 1;               break;
          case 's': options.soffset = atoi(optarg);     break;
          case 'S': 
            options.selfcall = 1;
            options.qname = messalloc(strlen(argv[optind])+1); strcpy(options.qname, argv[optind]);
            options.qlen = atoi(argv[optind+1]);
            options.sname = messalloc(strlen(argv[optind+2])+1); strcpy(options.sname, argv[optind+2]);
            options.slen = atoi(argv[optind+3]);
            dotterBinary = messalloc(strlen(argv[optind+4])+1);
            strcpy(dotterBinary, argv[optind+4]);
                                              break;
          case 'W': 
            options.winsize = messalloc(strlen(optarg)+1);
            strcpy(options.winsize, optarg);          break;
          case 'w': opts[1] = 'W';              break;
          case 'z': options.dotterZoom = atoi(optarg);  break;
          default : fatal("Illegal option");
      }
    }
  
    if (!options.savefile)
      {
	graphInit(&argc, argv);
	gexInit(&argc, argv);
      }

    /* Store X options for zooming in */
    {
	int i, len=0;
	extern char *Xoptions;
	
	for (i = optind+2; i < argc; i++)
	    len += strlen(argv[i])+1;
	
	Xoptions = messalloc(len+1);
	
	for (i = optind+2; i < argc; i++) {
	    strcat(Xoptions, argv[i]);
	    strcat(Xoptions, " ");
	}
    }

    if (options.selfcall) /* Dotter calling dotter */
      {
	qseq = (char *)messalloc(options.qlen+1);
	sseq = (char *)messalloc(options.slen+1);

        int l = fread(qseq, 1, options.qlen, stdin); 
	if (l != options.qlen) 
          {
	    fatal("Only read %d chars to qseq, expected %d", l, options.qlen);
          }
	qseq[options.qlen] = 0;

        l = fread(sseq, 1, options.slen, stdin);
	if (l != options.slen) 
          {
	    fatal("Only read %d chars to sseq, expected %d", l, options.slen);
          }
	sseq[options.slen] = 0;

	/* Read in MSPs */
	while (!feof (stdin))
          { 
	    char *cp;
	    if (!fgets (text, MAXLINE, stdin) || (unsigned char)*text == (unsigned char)EOF)
              {
                break;
              }
	    
	    if (!MSPlist) 
              {
		MSPlist = (MSP *)messalloc(sizeof(MSP));
		msp = MSPlist;
              }
	    else
              {
		msp->next = (MSP *)messalloc(sizeof(MSP));
		msp = msp->next;
              }
	    
	    cp = text;
	    while (*cp == ' ') cp++;
	    msp->type = strtol(cp, &cp, 10);
	    while (*cp == ' ') cp++;
	    msp->score = strtol(cp, &cp, 10);
	    while (*cp == ' ') cp++;
	    msp->fsColor = strtol(cp, &cp, 10);
	    while (*cp == ' ') cp++;
	    msp->qstart = strtol(cp, &cp, 10); 
	    while (*cp == ' ') cp++;
	    msp->qend = strtol(cp, &cp, 10); 
	    while (*cp == ' ') cp++;
	    msp->fs = NULL;
            strtol(cp, &cp, 10);
	    msp->sname = stringUnprotect(&cp, NULL);
	    stringUnprotect(&cp, msp->sframe);
	    msp->qname = stringUnprotect(&cp, NULL);
	    stringUnprotect(&cp, msp->qframe);
	    msp->desc = stringUnprotect(&cp, NULL);


	    /* really horrible hack */
	    if (msp->type == BLXMSP_FS_SEG)
	      {
		insertFS(msp, "chain_separator");
		opts[4] = 'L';
	      }
          }
	fclose(stdin);
      }
    else if (!options.seqInSFS)
      {
	if (argc - optind < 2) 
          {
	    fprintf(stderr, "%s\n", usage); 
	    exit(1);
          }
	else if(!(qfile = fopen(argv[optind], "r"))) 
          {
	    fprintf(stderr,"Cannot open %s\n", argv[optind]); exit(1);
          }
      
	qfilename = argv[optind];
	fseek(qfile, 0, SEEK_END);
	options.qlen = ftell(qfile);
	qseq = (char *)messalloc(options.qlen+1);
	fseek(qfile, 0, SEEK_SET);
      
	if ((cp = (char *)strrchr(argv[optind], '/')))
          {
            options.qname = strnew(cp+1, 0);
          }
	else
          {
            options.qname = strnew(argv[optind], 0);
          }

	if (!(sfile = fopen(argv[optind+1], "r"))) 
          {
	    fprintf(stderr,"Cannot open %s\n", argv[optind+1]); exit(1);
          }
      
	sfilename = argv[optind+1];
	fseek(sfile, 0, SEEK_END);
	options.slen = ftell(sfile);
	sseq = (char *)messalloc(options.slen+1);
	fseek(sfile, 0, SEEK_SET);
      
	if ((cp = (char *)strrchr(argv[optind]+1, '/')))
          {
            options.sname = strnew(cp+1, 0);
          }
	else
          {
            options.sname = strnew(argv[optind+1], 0);
          }

	/* Read in the sequences */
	int l = 0, count = 0;
	cc = qseq;
      
	while (!feof(qfile))
          {
            if (!fgets(line, MAXLINE, qfile))
              {
                break;
              }
        
            if ((cq = (char *)strchr(line, '\n')))
              {
                *cq = 0;
              }

            /* Name headers */
            if ((cq = (char *)strchr(line, '>'))) 
              {
                cq++;
                if (++l == 1) 
                  {
                    options.qname = messalloc(strlen(cq)+1); strNamecpy(options.qname, cq);
                    firstdesc = messalloc(strlen(cq)+1);
                    strcpy(firstdesc, cq);
                  }
                else
                  {
                  /* Multiple sequences - add break lines */
                    if (l == 2) 
                      {
                        opts[4] = 'L';

                        /* Second sequence - add break line to mark first sequence */
                        addBreakline (&MSPlist, qfilename, firstdesc, 0, '1');
                        
                        /* change sequence name to filename */
                        options.qname = messalloc(strlen(qfilename)+1); strcpy(options.qname, qfilename);
                      }
                    
                    addBreakline (&MSPlist, qfilename, cq, count, '1');
                  }
              }
            else 
              {
                for (cq = line; *cq; cq++) if (isalpha(*cq)) 
                  {
                    *cc++ = *cq;
                    count++;
                  }
              }
          }
      
	*cc = 0;
	
	l = 0, count = 0;
	cc = sseq;
      
	while (!feof(sfile))
          {
	    if (!fgets(line, MAXLINE, sfile))
	      break;
	    if ((cq = (char *)strchr(line, '\n')))
	      *cq = 0;

	    /* Name headers */
	    if ((cq = (char *)strchr(line, '>'))) {
		cq++;
		if (++l == 1) {
		    options.sname = messalloc(strlen(cq)+1); strNamecpy(options.sname, cq);
		    firstdesc = messalloc(strlen(cq)+1);
		    strcpy(firstdesc, cq);
		}
		else {
		  /* Multiple sequences - add break lines */

		    if (l == 2) {
			opts[4] = 'L';

		        /* Second sequence - add break line to mark first sequence */
		        addBreakline (&MSPlist, sfilename, firstdesc, 0, '2');
			
			/* change sequence name to filename */
			options.sname = messalloc(strlen(sfilename)+1); strcpy(options.sname, sfilename);
		    }
		    addBreakline (&MSPlist, sfilename, cq, count, '2');
		}
	    }
	    else 
              {
		for (cq = line; *cq; cq++) if (isalpha(*cq)) 
                  {
		    *cc++ = *cq;
		    count++;
                  }
              }
          }
        
	*cc = 0;
      }

    if (options.FSfilename) 
      {
	char dummyopts[32];	/* opts have different meaning in blixem */
	FILE *file;
	
	if (!strcmp(options.FSfilename, "-")) 
          {
	    file = stdin;
          }
	else if(!(file = fopen(options.FSfilename, "r")))
          {
	    messcrash("Cannot open %s\n", options.FSfilename);
          }
	
	parseFS(&MSPlist, file, dummyopts, &qseq, options.qname, &sseq, options.sname, options.qoffset);
      }

    /* Determine sequence types */
    if (Seqtype(qseq) == 'P' && Seqtype(sseq) == 'P') 
      {
	printf("\nDetected sequence types: Protein vs. Protein\n");
	type = 'P';
      }
    else if (Seqtype(qseq) == 'N' && Seqtype(sseq) == 'N') 
      {
	printf("\nDetected sequence types: DNA vs. DNA\n");
	type = 'N';
      }
    else if (Seqtype(qseq) == 'N' && Seqtype(sseq) == 'P') 
      {
	printf("\nDetected sequence types: DNA vs. Protein\n");
	type = 'X';
      }
    else fatal("Illegal sequence types: Protein vs. DNA - turn arguments around!\n\n%s", usage);

    if (options.revcompq) 
      {
	if (type != 'X')
          {
	    fatal("Revcomp'ing horizontal_sequence only needed in DNA vs. Protein");
          }
	else 
          {
	    cc = messalloc(options.qlen+1);
	    revComplement(cc, qseq);
	    messfree(qseq);
	    qseq = cc;
          }
      }

    /* Add -install for private colormaps */
    if (options.install) 
      {
        argvAdd(&argc, &argv, "-install");
      }

    if (!options.savefile)
      {
	dotter(type, opts, options.qname, qseq, options.qoffset, options.sname, sseq, options.soffset, 
	       0, 0, options.savefile, options.loadfile, options.mtxfile, options.memoryLimit, 
               options.dotterZoom, MSPlist, 0, options.winsize, options.pixelFacset) ;

	graphLoop(0) ;

	graphFinish () ;
      }
    else
      {
	/* stop graphical dialog boxes in batch mode. */
	struct messContextStruct nullContext = { NULL, NULL };
	messOutRegister(nullContext);
	messErrorRegister (nullContext);
	messExitRegister(nullContext);
	messCrashRegister (nullContext);
	dotter(type, opts, options.qname, qseq, options.qoffset, options.sname, sseq, options.soffset, 
	       0, 0, options.savefile, options.loadfile, options.mtxfile, options.memoryLimit, 
               options.dotterZoom, MSPlist, 0, options.winsize, options.pixelFacset);
      }

    return (0) ;
}
