/*  File: blxmain.c
 *  Author: Erik Sonnhammer
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
 * Description: BLXMAIN - Standalone calling program for blxview.
 * Exported functions: 
 * HISTORY:
 * Last edited: May 26 17:13 2009 (edgrif)
 * * Aug 26 16:57 1999 (fw): added this header
 * Created: Thu Aug 26 16:56:45 1999 (fw)
 * CVS info:   $Id: blxmain.c,v 1.14 2010-06-14 11:14:39 gb10 Exp $
 *-------------------------------------------------------------------
 */

#include <SeqTools/blixem_.h>
#include <SeqTools/utilities.h>
#include <wh/regular.h>
#include <wh/graph.h>
#include <wh/gex.h>


/* scrolled message window max messages. */
enum {BLIXEM_MESG_NUM = 50} ;

static void msgPopupsCB(BOOL msg_list) ;


/* Some globals.... */


/* global debug flag for blixem, set TRUE to get debugging output.           */
BOOL blixem_debug_G = FALSE ;


static char *usageText ="\n\
 Blixem - display Blast matches as a multiple alignment.\n\
\n\
 Reference:  Sonnhammer ELL & Durbin R (1994). A workbench for Large Scale\n\
 Sequence Homology Analysis. Comput. Applic. Biosci. 10:301-307.\n\
\n\
 Copyright (c) 2009-2010: Genome Research Ltd.\n\
\n\
\n\
  Usage: blixem [options] <sequencefile> <datafile> [X options] \n\
\n\
 Both <sequencefile> and <datafile> can be substituted by \"-\"\n\
 for reading from stdin (pipe).  If <sequencefile> is piped, the first\n\
 line should contain the sequence name and the second the sequence itself.\n\
\n\
 Options:\n\
 -s <mode>  Sorting mode at startup.\n\
               s = by Score\n\
               i = by Identity\n\
               n = by Name\n\
               p = by Position\n\
 -I         Inverted sorting order\n\
 -a         specify a string giving the names of the alignments, e.g. \"EST_mouse EST_human\" etc.\n\
 -b         Don't start with Big Picture.\n\
 -c <file>  Read configuration options from 'file'.\n\
 -S <#>     Start display at position #.\n\
 -F <file>  Read in query sequence and data from <file> (replaces sequencefile).\n\
 -h         Help and more options.\n\
 -o <optstring>\n\
            Blixem options,  e.g. -o \"MBr\" you'll have to read the source code for details.\n\
 -O <#>     sequence offset.\n\
 -P nodeid<:port>\n\
            Causes Blixem to get sequences from a pfetch server at machine nodeid on the given port (default 22100).\n\
 -r         Remove input files after parsing, used by xace when calling blixem as a\n\
            standalone program.\n\
 -x <file>  Read in extra data (to that in datafile) from from <file>, data may be in a\n\
            different format.\n\
\n\
 Some X options:\n\
 -acefont <font> Main font.\n\
 -font    <font> Menu font.\n\
\n\
 To make the datafile from blast output, run MSPcrunch with option -q.\n\n\
\n\
"BLIXEM_AUTHOR_TEXT"\n\
 Version" ;

static char *help_string = "\n\
 o To pipe MSPcrunch output directly to Blixem, use \"-\"\n\
   as the second parameter ([datafile]).  Example:\n\
\n\
   MSPcrunch -q <my.blast_output> | blixem <my.seq> -\n\
\n\
\n\
 o The BLAST program (blastp, blastn, blastx, tblastn, tblastx)\n\
   is automatically detected from the Blast output header by MSPcrunch\n\
   and is passed on to Blixem in the seqbl format (-q).  If\n\
   your output lacks the header, set the program with these options:\n\
\n\
   -p    Use blastp output (blastx is default)\n\
   -n    Use blastn output\n\
   -t    Use tblastn output\n\
   -l    Use tblastx output\n\n" ;




/* Entry point for blixem standalone program, you should be aware when
 * altering this that blxview.c is also compiled into acedb and that
 * blxview() is called directly by acedb code.
 */
int main(int argc, char **argv)
{
  FILE        *seqfile, *FSfile;
  char seqfilename[1000] = {'\0'};
  char FSfilename[1000] = {'\0'};
  char *refSeq = NULL;
  char *dummyseq = NULL;    /* Needed for blxparser to handle both dotter and blixem */
  char line[MAXLINE+1];
  MSP *mspList = NULL ;
  char opts[BLXOPT_NUM_OPTS]=" MBr Z   ";
  char refSeqName[FULLNAMESIZE+1] = "";
  char dummyseqname[FULLNAMESIZE+1] = "";
  
  int displayStart = 0;
  int qOffset = 0;
  int install = 1;
  gboolean singleFile = FALSE;        /* if true, there is a single file containing both the query 
                                         sequence and the alignment data */
  
  int          optc;
  extern int   optind;
  extern char *optarg;
  char        *optstring="a:bc:F:hIilno:O:pP:rS:s:tx:";
  char *usage;
  BOOL rm_input_files = FALSE ;
  PfetchParams *pfetch = NULL ;
  BOOL xtra_data = FALSE ;
  FILE *xtra_file = NULL ;
  char xtra_filename[1000] = {'\0'} ;
  char *align_types = NULL ;
  char *config_file = NULL ;
  GError *error = NULL ;
 
  /* Set up the GLib message handlers
   *
   * For now, also set up the acedb message package. I'm trying to phase out use of the acedb
   * message package so we can get rid of dependencies on acedb: for now, the GLib message
   * handlers just pass the messages to the acedb package. This is so that we can gradually switch
   * everything over to using GLib calls instead of calling the acedb package directly. Then when
   * we're ready we can just change the GLib handlers to use something else instead of acedb. 
   * 
   * There are two handlers: the default one for all non-critical messages, which will just log
   * output to the console, and one for critical messages and errors, which will display a 
   * pop-up message (the idea being that we don't bother the user unless it's something serious).
   * So, to get a pop-up message use g_critical, and to log a message or warning use g_message, 
   * g_warning, g_debug etc. Note that g_error is always fatal.
   */
  g_log_set_default_handler(defaultMessageHandler, NULL);
  g_log_set_handler(NULL, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL, popupMessageHandler, NULL);

  messSetMsgInfo(argv[0], BLIXEM_TITLE_STRING) ;


  /* Stick version info. into usage string. */
  usage = g_malloc(strlen(usageText) + strlen(blixemVersion) + 10) ;
  sprintf(usage, "%s %s\n", usageText, blixemVersion) ;


  while ((optc = getopt(argc, argv, optstring)) != EOF )
    {
      switch (optc)
	{
	case 'a':
	  align_types = hprintf(0, "%s", optarg) ;
	  break;
	case 'b':
	  opts[2] = 'b';
	  break;
	case 'c': 
	  config_file = strnew(optarg, 0) ;
	  break;
	case 'F': 
	  singleFile = TRUE;        
	  strcpy(FSfilename, optarg);
	  break;
	case 'h': 
	  fprintf(stderr, "%s\n%s\n", usage, help_string) ;
	  exit(EXIT_FAILURE) ;
	  break;
	case 'I':
	  opts[6] = 'I';
	  break;
	case 'i':
	  install = 0;
	  break;
	case 'l':
	  opts[BLXOPT_MODE] = 'L';
	  break;
	case 'n':
	  opts[BLXOPT_MODE] = 'N';
	  break;
	case 'o':
	  {
	    int i ;

	    memset(opts, ' ', strlen(opts)) ;

	    for (i = 0 ; i < strlen(optarg) ; i++)
	      {
		opts[i] = optarg[i] ;
	      }
	  }
	  break;
	case 'O':
	  /* Note that I don't do any checking of the value, this is in line */
	  /* with the way blxview() accepts any value of offset, I just check*/
	  /* its actually an integer.                                        */
	  if (!(utStr2Int(optarg, &qOffset)))
	    g_error("Bad offset argument: \"-%c %s\"\n", optc, optarg) ;
	  break;
	case 'p':
	  opts[BLXOPT_MODE] = 'P';
	  break;
	case 'P':
	  pfetch = g_malloc(sizeof(PfetchParams)) ;
	  pfetch->net_id = strtok(optarg, ":") ;
	  pfetch->port = atoi(strtok(NULL, ":")) ;
	  break;
	case 'r':
	  rm_input_files = TRUE ;
	  break ;
	case 'S': 
	  if (!(displayStart = atoi(optarg)))
	    g_error("Bad diplay start position: %s\n", optarg); 
	  opts[1] = ' ';
	  break;
	case 's': 
	  if (*optarg != 's' && *optarg != 'i' && *optarg != 'n' && *optarg != 'p')
	    {
	      fprintf(stderr,"Bad sorting mode: %s\n", optarg); 
	      exit(EXIT_FAILURE) ;
	    }
	  opts[BLXOPT_SORT_MODE] = *optarg;
	  break;
	case 't':
	  opts[BLXOPT_MODE] = 'T';
	  break;
	case 'x': 
	  xtra_data = TRUE ;
	  strcpy(xtra_filename, optarg);
	  break;

	default : g_error("Illegal option\n");
	}
    }

  if (argc-optind < 2 && !singleFile)
    {
      fprintf(stderr, "%s\n", usage); 

      exit(EXIT_FAILURE);
    }


  /* Add -install for private colormaps */
  if (install)
    argvAdd(&argc, &argv, "-install");

  gtk_init(&argc, &argv);

  /* Old style graph init is still required for calling dotter from blixem */
  graphInit(&argc, argv) ;
  gexInit(&argc, argv);
  gexSetMessPrefs(FALSE, BLIXEM_MESG_NUM, msgPopupsCB) ;  /* for control of popup or mesglist */

  /* Set up program configuration. */
  if (!blxInitConfig(config_file, &error))
    {
      g_error("Config File Error: %s\n", error->message) ;
    }


  if (!singleFile)
    {
      /* The query sequence is in a separate file. */
      
      strcpy(seqfilename, argv[optind]);
      strcpy(FSfilename, argv[optind+1]);

      if (!strcmp(seqfilename, "-"))
	{
	  seqfile = stdin;
	}
      else if(!(seqfile = fopen(seqfilename, "r")))
	{
	  messerror("Cannot open %s\n", seqfilename);
	}
	
      /* Read in query sequence */
      if (seqfile == stdin)
	{
	  /* Same file for query and FS data - sequence on first two lines */
	  Array arr = arrayCreate(5000, char);
	  int i=0;
	    
	  if (!fgets(line, MAXLINE, seqfile))
	    g_error("Error reading seqFile.\n");

	  sscanf(line, "%s", refSeqName);

	  /* Scan the input stream and place the chars in an array */
	  char charFromStream = fgetc(seqfile);
	  while (charFromStream != '\n')
	    {
	      array(arr, i++, char) = charFromStream;
	    }
	  
	  /* Allocate the reference sequence from the array of chars */
	  refSeq = g_malloc(arrayMax(arr)+1);
	  char *refSeqChar = refSeq;
	  
	  for (i = 0; i < arrayMax(arr);)
	    {
	      *refSeqChar++ = arr(arr, i++, char);
	    }
	    
	  arrayDestroy(arr);
	}
      else
	{
	  fseek(seqfile, 0, SEEK_END);
	  refSeq = g_malloc(ftell(seqfile)+1);
	  char *refSeqChar = refSeq;
	  fseek(seqfile, 0, SEEK_SET);
	  
	  while (!feof(seqfile))
	    {
	      if (!fgets(line, MAXLINE, seqfile))
		{
		  break;
		}
		  
	      char *newLineCharPtr = (char *)strchr(line, '\n');
	      if (newLineCharPtr)
		{
		  *newLineCharPtr = 0;
		}
		    
	      if (*line == '>')
		{
		  strncpy(refSeqName, line+1, 255);
		  refSeqName[255]=0;
		  
		  char *spaceCharPtr = (char *)strchr(refSeqName, ' ');
		  if (spaceCharPtr)
		    {
		      *spaceCharPtr = 0;
		    }
		  
		  continue;
		}

	      char *lineChar = NULL;
	      for (lineChar = line; *lineChar; lineChar++) 
		{
		  *refSeqChar++ = *lineChar;
		}
		
	      *refSeqChar = 0;
	    }
	  fclose(seqfile);
	}
    }

  /* Parse the data file containing the homol descriptions.                */
  if (!strcmp(FSfilename, "-"))
    {
      FSfile = stdin;
    }
  else if(!(FSfile = fopen(FSfilename, "r")))
    {
      g_error("Cannot open %s\n", FSfilename);
    }
  
  GList *seqList = NULL; /* parser compiles a list of BlxSequences into this list */
  parseFS(&mspList, FSfile, opts, &seqList, &refSeq, refSeqName, &dummyseq, dummyseqname, qOffset) ;
  
  if (FSfile != stdin)
    {
      fclose(FSfile) ;
    }

  /* There may an additional file containing homol data in an alternative format. */
  if (xtra_data)
    {
      if(!(xtra_file = fopen(xtra_filename, "r")))
	{
	  g_error("Cannot open %s\n", xtra_filename) ;
	}
      
      parseFS(&mspList, xtra_file, opts, &seqList, &refSeq, refSeqName, &dummyseq, dummyseqname, qOffset) ;
      fclose(xtra_file) ;
    }


  /* Remove the input files if requested.                                  */
  if (rm_input_files)
    {
      if(seqfilename[0] != '\0' && unlink(seqfilename) != 0)
	g_warning("Unlink of sequence input file \"%s\" failed: %s\n",
		  seqfilename, messSysErrorText()) ;
      if(FSfilename[0] != '\0' && unlink(FSfilename) != 0)
	g_warning("Unlink of MSP input file \"%s\" failed: %s\n",
		  FSfilename, messSysErrorText()) ;
      if (xtra_filename[0] != '\0' && unlink(xtra_filename) != 0)
	g_warning("Unlink of extra MSP sequence input file \"%s\" failed: %s\n",
		  xtra_filename, messSysErrorText()) ;
    }

  /* Now display the alignments, this call does not return. (Note that
   * TRUE signals blxview() that it is being called from this standalone
   * blixem program instead of as part of acedb. */
  blxview(refSeq, refSeqName, displayStart, qOffset, mspList, seqList, opts, pfetch, align_types, TRUE);
  gtk_main();
    
  /* We should not get here.... */
  return (EXIT_FAILURE) ;
}




/* 
 *                    Internal functions.
 */




/* Enables user to switch between seeing informational messages in popups
 * or in a scrolled window. */
static void msgPopupsCB(BOOL msg_list)
{
  int list_length ;

  if (msg_list)
    list_length = BLIXEM_MESG_NUM ;
  else
    list_length = 0 ;

  gexSetMessPopUps(msg_list, list_length) ;

  return ;
}





