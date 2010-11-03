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
 * CVS info:   $Id: dotterMain.c,v 1.21 2010-11-03 13:39:09 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <SeqTools/utilities.h>
#include <SeqTools/blxGff3Parser.h>
#include <SeqTools/blxparser.h>
#include <SeqTools/dotter_.h>
#include <SeqTools/blxmsp.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#define UNSET_INT  -1


static void setDefaultOptions(DotterOptions *options)
{
  options->qoffset = 0;
  options->soffset = 0;
  options->selfcall = FALSE;
  options->qlen = UNSET_INT;
  options->slen = UNSET_INT;
  options->dotterZoom = 0;
  options->install = 1;
  options->pixelFacset = 0;
  options->seqInSFS = 0;
  
  options->memoryLimit = 0.0;
  
  options->savefile = NULL;
  options->loadfile = NULL;
  options->FSfilename = NULL;
  options->mtxfile = NULL;
  options->winsize = NULL;
  options->qname = NULL;
  options->sname = NULL;
  
  options->mirrorImage = TRUE;
  options->watsonOnly = FALSE;
  options->crickOnly = FALSE;
  options->hspsOnly = FALSE;
  options->swapGreyramp = FALSE;
  options->fsEndLinesOn = FALSE;
  options->hspGaps = FALSE;
  options->hozScaleRev = FALSE;
  options->vertScaleRev = FALSE;
}


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


static void addBreakline (MSP **MSPlist, char *name, char *desc, int pos, const int sFrame)
{
   MSP   
       *msp;
   char 
       *cp;

   if (!*MSPlist) {
       *MSPlist = (MSP *)g_malloc(sizeof(MSP));
       msp = *MSPlist;
   }
   else {
     msp = *MSPlist;
     while(msp->next) msp = msp->next;
     msp->next = (MSP *)g_malloc(sizeof(MSP));
     msp = msp->next;
   }

   msp->qname = g_malloc(strlen(name)+1);
   strcpy(msp->qname, name);

   msp->desc = g_malloc(strlen(desc)+1);
   strcpy(msp->desc, desc);
   if ((cp = (char *)strchr(msp->desc, ' ')))
     *cp = 0;
   if ((cp = (char *)strchr(msp->desc, '\n')))
     *cp = 0;

   msp->qRange.min = msp->qRange.max = pos;
//   msp->sFrame = sFrame;
   msp->fsColor = 0;
   msp->type = BLXMSP_FS_SEG;
   msp->score = 100.0;
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
 -v             Reverse and complement vertical_sequence (DNA vs Protein)\n\
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
  
  usage = g_malloc(strlen(usageText) + strlen(dotterVersion) + strlen(cc_date) + 20);
  sprintf(usage, "%s%s, compiled %s\n", usageText, dotterVersion, cc_date);
  
  return usage;
}


int main(int argc, char **argv)
{
  DotterOptions options;
  setDefaultOptions(&options);

  char   
      *qseq=0, *sseq=0, line[MAXLINE+1],
      *curChar, *cc, *cq, type, 
      *firstdesc, *qfilename, *sfilename,
      text[MAXLINE+1];

  FILE *qfile, *sfile;
  MSP *MSPlist = NULL;
  GList *seqList = NULL;
  
  /* MSPlist above is obsolete and should be replaced by featureLists, which contains all the MSPs
   * but in GLists in an array indexed by type. Initialise each GList to NULL. */
  GList* featureLists[BLXMSP_NUM_TYPES];
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    {
      featureLists[typeId] = NULL;
    }
  
  int          optc;
  extern int   optind;
  extern char *optarg;
  char        *optstring="b:cDf:F:Hil:M:m:p:q:Rrs:SvW:wz:";

  static char *dotterBinary = NULL;

  char *usage = getUsageText();

  /* The strand stuff is a bit hacky, because dotter was originally never designed to deal with
   * reverse match seq strands, so match and ref seq strands work in different ways. If the ref seq
   * strand is reversed then the horizontal scale is reversed as well. (Because the 'active' (top) strand
   * in Blixem is always the reverse strand if the display is reversed. Note that for DNA matches both 
   * strands are always shown anyway, so the -r option essentially just swaps which strand is shown at 
   * the top of the alignment tool, and of course reverses the direction of the display.)
   * The vertical scale was never originally designed to be reversed in dotter. I've 
   * added the vertScaleRev flag in case we might want to do this in the future, but this is currently
   * always set to false, even if we have the reverse match seq strand (which is indicated with the -v option). */
  BlxStrand qStrand = BLXSTRAND_FORWARD;
  BlxStrand sStrand = BLXSTRAND_FORWARD;

  while ((optc = getopt(argc, argv, optstring)) != EOF)
    {printf("%c", optc);
      switch (optc) 
        {
          case 'b': 
            options.savefile = g_malloc(strlen(optarg)+1);
            strcpy(options.savefile, optarg);           break;
          case 'c': options.crickOnly = TRUE;           break;
          case 'D': options.mirrorImage = FALSE;        break;
          case 'f': 
            options.FSfilename = g_malloc(strlen(optarg)+1);
            strcpy(options.FSfilename, optarg);         break;
          case 'F': 
            options.seqInSFS = 1;        
            options.FSfilename = g_malloc(strlen(optarg)+1);
            strcpy(options.FSfilename, optarg);         break;
          case 'H': options.hspsOnly = TRUE;            break;
          case 'i': options.install = 0;                break;
          case 'l': 
            options.loadfile = g_malloc(strlen(optarg)+1);
            strcpy(options.loadfile, optarg);           break;
          case 'M': 
            options.mtxfile = g_malloc(strlen(optarg)+1);
            strcpy(options.mtxfile, optarg);            break;
          case 'm': options.memoryLimit = atof(optarg); break;
          case 'p': options.pixelFacset = atoi(optarg); break;
          case 'q': options.qoffset = atoi(optarg);     break;
          case 'R': options.swapGreyramp = TRUE;        break;
          case 'r': 
	    options.hozScaleRev = TRUE;
	    qStrand = BLXSTRAND_REVERSE;		break;
          case 's': options.soffset = atoi(optarg);     break;
          case 'S': 
            options.selfcall = TRUE;
            options.qname = g_malloc(strlen(argv[optind])+1); strcpy(options.qname, argv[optind]);
            options.qlen = atoi(argv[optind+1]);
            options.sname = g_malloc(strlen(argv[optind+2])+1); strcpy(options.sname, argv[optind+2]);
            options.slen = atoi(argv[optind+3]);
            dotterBinary = g_malloc(strlen(argv[optind+4])+1);
            strcpy(dotterBinary, argv[optind+4]);
                                                        break;
          case 'v': sStrand = BLXSTRAND_REVERSE;	break;
          case 'W': 
            options.winsize = g_malloc(strlen(optarg)+1);
            strcpy(options.winsize, optarg);            break;
          case 'w': options.watsonOnly = TRUE;          break;
          case 'z': options.dotterZoom = atoi(optarg);  break;
          default : g_error("Illegal option");
      }
    }

    if (!options.savefile)
      {
	gtk_init(&argc, &argv);
      }

    /* Store X options for zooming in */
    {
	int i, len=0;
	extern char *Xoptions;
	
	for (i = optind+2; i < argc; i++)
	    len += strlen(argv[i])+1;
	
	Xoptions = g_malloc(len+1);
	
	for (i = optind+2; i < argc; i++) {
	    strcat(Xoptions, argv[i]);
	    strcat(Xoptions, " ");
	}
    }

    if (options.selfcall) /* Blixem/Dotter calling dotter */
      {
	qseq = (char *)g_malloc(options.qlen+1);
	sseq = (char *)g_malloc(options.slen+1);

        int l = fread(qseq, 1, options.qlen, stdin); 
	if (l != options.qlen) 
          {
	    g_error("Only read %d chars to qseq, expected %d", l, options.qlen);
          }
	qseq[options.qlen] = 0;

        l = fread(sseq, 1, options.slen, stdin);
	if (l != options.slen) 
          {
	    g_error("Only read %d chars to sseq, expected %d", l, options.slen);
          }
	sseq[options.slen] = 0;

	/* Read in MSPs */
        DEBUG_OUT("Reading MSPs\n");
        MSP *lastMsp = NULL;
        
	while (!feof (stdin))
          { 
            /* read in the blxsequences. they are separated by newlines */
	    if (!fgets (text, MAXLINE, stdin) || (unsigned char)*text == (unsigned char)EOF)
              break;

	    int numMsps = 0;
            BlxSequence *blxSeq = readBlxSequenceFromText(text, &numMsps);
            seqList = g_list_append(seqList, blxSeq);
            
            /* read in the msps for this blx sequence and add them to the msplist */
            int i = 0;
            for ( ; i < numMsps; ++i)
              {
                /* msps are also separated by newlines */
                if (!fgets (text, MAXLINE, stdin) || (unsigned char)*text == (unsigned char)EOF)
                  break;

                MSP *msp = createEmptyMsp(&lastMsp, &MSPlist);
                readMspFromText(msp, text);
                
                /* add it to the relevant feature list. Use prepend because it is quicker */
                featureLists[msp->type] = g_list_prepend(featureLists[msp->type], msp);
                
                /* add the msp to the blxsequence */
                blxSeq->mspList = g_list_append(blxSeq->mspList, msp);
                msp->sSequence = blxSeq;
                
                /* really horrible hack */
                if (msp->type == BLXMSP_FS_SEG)
                  {
                    insertFS(msp, "chain_separator");
                    options.fsEndLinesOn = TRUE;
                  }
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
	qseq = (char *)g_malloc(options.qlen+1);
	fseek(qfile, 0, SEEK_SET);
      
	if ((curChar = (char *)strrchr(argv[optind], '/')))
          {
            options.qname = g_strdup(curChar+1);
          }
	else
          {
            options.qname = g_strdup(argv[optind]);
          }

	if (!(sfile = fopen(argv[optind+1], "r"))) 
          {
	    fprintf(stderr,"Cannot open %s\n", argv[optind+1]); exit(1);
          }
      
	sfilename = argv[optind+1];
	fseek(sfile, 0, SEEK_END);
	options.slen = ftell(sfile);
	sseq = (char *)g_malloc(options.slen+1);
	fseek(sfile, 0, SEEK_SET);
      
	if ((curChar = (char *)strrchr(argv[optind]+1, '/')))
          {
            options.sname = g_strdup(curChar+1);
          }
	else
          {
            options.sname = g_strdup(argv[optind+1]);
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
                    options.qname = g_malloc(strlen(cq)+1); strNamecpy(options.qname, cq);
                    firstdesc = g_malloc(strlen(cq)+1);
                    strcpy(firstdesc, cq);
                  }
                else
                  {
                  /* Multiple sequences - add break lines */
                    if (l == 2) 
                      {
                        options.fsEndLinesOn = TRUE;

                        /* Second sequence - add break line to mark first sequence */
                        addBreakline (&MSPlist, qfilename, firstdesc, 0, 1);
                        
                        /* change sequence name to filename */
                        options.qname = g_malloc(strlen(qfilename)+1); strcpy(options.qname, qfilename);
                      }
                    
                    addBreakline (&MSPlist, qfilename, cq, count, 1);
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
		    options.sname = g_malloc(strlen(cq)+1); strNamecpy(options.sname, cq);
		    firstdesc = g_malloc(strlen(cq)+1);
		    strcpy(firstdesc, cq);
		}
		else {
		  /* Multiple sequences - add break lines */

		    if (l == 2) {
			options.fsEndLinesOn = TRUE;

		        /* Second sequence - add break line to mark first sequence */
		        addBreakline (&MSPlist, sfilename, firstdesc, 0, 2);
			
			/* change sequence name to filename */
			options.sname = g_malloc(strlen(sfilename)+1); strcpy(options.sname, sfilename);
		    }
		    addBreakline (&MSPlist, sfilename, cq, count, 2);
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
	    g_error("Cannot open %s\n", options.FSfilename);
          }
	
        GSList *supportedTypes = blxCreateSupportedGffTypeList();
        IntRange qRange = {UNSET_INT, UNSET_INT};

	parseFS(&MSPlist, file, dummyopts, featureLists, &seqList, supportedTypes, NULL, &qseq, options.qname, &qRange, &sseq, options.sname);
        
        blxDestroyGffTypeList(&supportedTypes);
      }

    /* Determine sequence types */
    if (determineSeqType(qseq) == BLXSEQ_PEPTIDE && determineSeqType(sseq) == BLXSEQ_PEPTIDE) 
      {
	printf("\nDetected sequence types: Protein vs. Protein\n");
	type = 'P';
      }
    else if (determineSeqType(qseq) == BLXSEQ_DNA && determineSeqType(sseq) == BLXSEQ_DNA) 
      {
	printf("\nDetected sequence types: DNA vs. DNA\n");
	type = 'N';
      }
    else if (determineSeqType(qseq) == BLXSEQ_DNA && determineSeqType(sseq) == BLXSEQ_PEPTIDE) 
      {
	printf("\nDetected sequence types: DNA vs. Protein\n");
	type = 'X';
      }
    else
      {
        g_error("Illegal sequence types: Protein vs. DNA - turn arguments around!\n\n%s", usage);
      }
      
    /* Add -install for private colormaps */
    if (options.install) 
      {
        argvAdd(&argc, &argv, "-install");
      }

    if (!options.savefile)
      {
        /* Set the message handlers. (Do this here because we don't want graphical dialog
         * boxes for batch mode.) */
        g_log_set_default_handler(defaultMessageHandler, NULL);
        g_log_set_handler(NULL, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL, popupMessageHandler, NULL);
        
	dotter(type, &options, options.qname, qseq, options.qoffset, qStrand, options.sname, sseq, options.soffset, sStrand,
	       0, 0, options.savefile, options.loadfile, options.mtxfile, options.memoryLimit, 
               options.dotterZoom, MSPlist, seqList, 0, options.winsize, options.pixelFacset) ;

        gtk_main();
      }
    else
      {
        /* Batch mode */
	dotter(type, &options, options.qname, qseq, options.qoffset, qStrand, options.sname, sseq, options.soffset, sStrand,
	       0, 0, options.savefile, options.loadfile, options.mtxfile, options.memoryLimit, 
               options.dotterZoom, MSPlist, seqList, 0, options.winsize, options.pixelFacset);
      }

    return (0) ;
}
