/*  File: belvuMain.c
 *  Author: Gemma Barson, 2011-04-06
 *  Copyright (c) 2011 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * SeqTools is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
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
 * ---------------------------------------------------------------------------
 * This file is part of the SeqTools sequence analysis package, 
 * written by
 *      Gemma Barson      (Sanger Institute, UK)  <gb10@sanger.ac.uk>
 * 
 * based on original code by
 *      Erik Sonnhammer   (SBC, Sweden)           <Erik.Sonnhammer@sbc.su.se>
 * 
 * and utilizing code taken from the AceDB and ZMap packages, written by
 *      Richard Durbin    (Sanger Institute, UK)  <rd@sanger.ac.uk>
 *      Jean Thierry-Mieg (CRBM du CNRS, France)  <mieg@kaa.crbm.cnrs-mop.fr>
 *      Ed Griffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: main() function for Belvu application
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvu_.h"
#include "belvuApp/belvuWindow.h"
#include "belvuApp/belvuTree.h"
#include "seqtoolsUtils/utilities.h"
#include "seqtoolsUtils/blxmsp.h"
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>


/* Usage text. This is a duplicate of the text that is in 
 * doc/User_doc/dotter_usage.txt, so ideally we would get rid of this and use
 * the text from the file instead; for now, we must update both. */

#define USAGE_TEXT "\
\n\
 Belvu - View multiple alignments in good-looking colours.\n\
\n\
 Usage: belvu [options] <multiple_alignment>|- [X options]\n\
\n\
   <multiple_alignment>|- = alignment file or pipe.\n\
\n\
 Options:\n\
  -c          Print Conservation table.\n\
  -l <file>   Load residue color code file.\n\
  -L <file>   Load markup and organism color code file.\n\
  -m <file>   Read file with matching sequence segments.\n\
  -r          Read alignment in 'raw' format (Name sequence).\n\
  -R          Do not parse coordinates when reading alignment.\n\
  -o <format> Write alignment or tree to stdout in this format and exit.\n\
                Valid formats: MSF, Mul(Stockholm), Selex, \n\
                FastaAlign, Fasta, tree.\n\
  -X <cutoff> Print UPGMA-based subfamilies at cutoff <cutoff>.\n\
  -n <cutoff> Make non-redundant to <cutoff> %identity at startup.\n\
  -Q <cutoff> Remove columns more gappy than <cutoff>.\n\
  -q <cutoff> Remove sequences more gappy than <cutoff>.\n\
  -G          Penalize gaps in pairwise comparisons.\n\
  -i          Ignore gaps in conservation calculation.\n\
  -P          Remove partial sequences at startup.\n\
  -C          Don't write coordinates to saved file.\n\
  -z <char>   Separator char between name and coordinates in saved file.\n\
  -a          Show alignment annotations on screen (Stockholm format only).\n\
  -p          Output random model probabilites for HMMER.\n\
              (Based on all residues.)\n\
  -S <order>  Sort sequences in this order.\n\
                a -> alphabetically\n\
                o -> by Swissprot organism, alphabetically\n\
                s -> by score\n\
                n -> by Neighbor-joining tree\n\
                u -> by UPGMA tree\n\
                S -> by similarity to first sequence\n\
                i -> by identity to first sequence\n\
  -s <file>   Read in file of scores.\n\
  -T <method> Tree options:\n\
                i -> Start up showing tree\n\
                I -> Start up showing only tree\n\
                d -> Show distances in tree\n\
                n -> Neighbor-joining\n\
                u -> UPGMA\n\
                c -> Don't color tree by organism\n\
                o -> Don't display sequence coordinates in tree\n\
                b -> Use Scoredist distance correction (default)\n\
                j -> Use Jukes-Cantor distance correction\n\
                k -> Use Kimura distance correction\n\
                s -> Use Storm & Sonnhammer distance correction\n\
                r -> Use uncorrected distances\n\
                p -> Print distance matrix and exit\n\
                R -> Read distance matrix instead of alignment\n\
              (only in combination with Tree routines)\n\
  -b <n>      Apply boostrap analysis with <n> bootstrap samples\n\
  -B          Print out bootstrap trees and exit\n\
              (Negative value -> display bootstrap trees on screen)\n\
  -O <label>  Read organism info after this label (default OS)\n\
  -t <title>  Set window title.\n\
  -g          Draw grid line (for debugging).\n\
  -u          Start up with uncoloured alignment (faster).\n\
  -h, --help  Show this usage information\n\
  --compiled  Show package compile date\n\
  --version   Show package version number\n\
\n\
 Some X options:\n\
  -acefont <font>   Main font.\n\
  -font    <font>   Menu font.\n\
\n\
 Note: X options only work after \"setenv POSIXLY_CORRECT\"\n\
\n\
 setenv BELVU_FETCH to desired sequence fetching program.\n\
 setenv BELVU_FONT_SIZE to specify window font size.\n\n\
"


/* Text to show the authors, version and compile date */
#define FOOTER_TEXT "\
-----\n\
"AUTHOR_TEXT_FULL" \n\
\n\
 Reference: Scoredist: A simple and robust protein sequence distance estimator.\n\
            Erik LL Sonnhammer and Volker Hollich.\n\
            BMC Bioinformatics 6:108 (2005)\n\
\n\
 See http://www.sanger.ac.uk/resources/software/seqtools/ for more info.\n\
\n\
 "BELVU_COPYRIGHT_STRING"\n\
 "BELVU_LICENSE_STRING"\n\
\n\
 Version "BELVU_VERSION_COMPILE"\n\
\n\
"

/* Text to show the version */
#define VERSION_TEXT BELVU_PACKAGE_VERSION"\n"



/* Prints usage info to stderr */
static void showUsageText()
{
  /* Pring usage info followed by authors */
  fprintf(stderr, "%s%s", USAGE_TEXT, FOOTER_TEXT);
}

/* Prints version info to stderr */
static void showVersionInfo()
{
  fprintf(stderr, VERSION_TEXT);  
}

/* Prints compiled date (must go to stdout for our build scripts to work) */
static void showCompiledInfo()
{
  fprintf(stdout, "%s\n", UT_MAKE_COMPILE_DATE());  
}


/* Convert swissprot name suffixes to organisms */
static void suffix2organism(GArray *alignArr, GArray *organismArr) 
{
  int i = 0;
  char *cp = NULL;
  
  for (i = 0; i < alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(alignArr, ALN, i);
      
      if (!alnp->markup && (cp = strchr(alnp->name, '_'))) 
        {
          char *suffix = g_malloc(strlen(cp) + 1);
          strcpy(suffix, cp + 1);
          
          /* Add organism to table of organisms.  This is necessary to make all 
           sequences of the same organism point to the same place and to make a 
           non-redundant list of present organisms */
          alnp->organism = suffix;
          
          /* Only insert it if this organism is not already in the array */
          int ip = 0;
          if (!arrayFind(organismArr, alnp, &ip, (void*)organism_order))
            {
              g_array_append_val(organismArr, *alnp);
              g_array_sort(organismArr, organism_order);
            }
          
          /* Store pointer to unique organism in ALN struct */
          ip = 0;
          if (arrayFind(organismArr, alnp, &ip, (void*)organism_order))
            {
              ALN *alnTmp = &g_array_index(organismArr, ALN, ip);
              alnp->organism = alnTmp->organism;
            }
        }
    }
}


/* This is to read in sequence names and count sequences */
static int treeReadDistancesNames(BelvuContext *bc)
{
  char   *cp = NULL;
  ALN     aln;
  initAln(&aln);
  
  char line[MAXLENGTH + 1];
  
  if (!fgets (line, MAXLENGTH, bc->treeReadDistancesPipe))
    g_error("Error reading distance matrix");
  
  if ((cp = strchr(line, '\n')))
    *cp = 0;
  
  int nseq = 0;
  
  while ((cp = strtok(nseq ? 0 : line, " \t")))
    {
      initAln(&aln);
      
      strncpy(aln.name, cp, MAXNAMESIZE);
      aln.name[MAXNAMESIZE] = 0;
      aln.nr = nseq;
      g_array_insert_val(bc->alignArr, nseq, aln);
      
      nseq++;
      
      /* printf("%d  %s\n", aln.nr, aln.name); */
    }
  
  g_array_sort(bc->alignArr, nrorder);
  
  return nseq;
}


static void readScores(char *filename, BelvuContext *bc)
{
  char line[MAXLENGTH+1], linecp[MAXLENGTH+1], *cp;
  FILE *file;
  int scoreLen;
  
  gboolean found = FALSE;
  gboolean warnings = FALSE;

  ALN aln;
  initAln(&aln);

  if (!(file = fopen(filename, "r")))
    g_error("Cannot open file %s", filename);
  
  while (!feof (file))
    { 
      if (!fgets (line, MAXLENGTH, file)) break;
      strcpy(linecp, line);
      
      initAln(&aln);
      
      if (!(cp = strtok(line, " "))) 
	g_error("Error parsing score file %s.\nLine: %s", filename, linecp);
      
      if (!(sscanf(cp, "%f", &aln.score)))
        g_error("Error parsing score file %s - bad score.\nLine: %s", filename, linecp);
      
      if (!(cp = strtok(0, "/"))) 
	g_error("Error parsing score file %s.\nLine: %s", filename, linecp);
      
      strncpy(aln.name, cp, MAXNAMESIZE);
      aln.name[MAXNAMESIZE] = 0;
      
      if (!(cp = strtok(0, "-"))) 
	g_error("Error parsing score file %s.\nLine: %s", filename, linecp);
      
      if (!(aln.start = atoi(cp)))
        g_error("Error parsing score file %s - no start coordinate.\nLine: %s", filename, linecp);
      
      if (!(cp = strtok(0, "\n"))) 
	g_error("Error parsing score file %s.\nLine: %s", filename, linecp);
      
      if (!(aln.end = atoi(cp)))
        g_error("Error parsing score file %s - no end coordinate.\nLine: %s", filename, linecp);
      
      int idx = 0;
      if (!alignFind(bc->alignArr, &aln, &idx)) 
        {
	  warnings = TRUE;
          /* printf("Warning: %s/%d-%d (score %.1f) not found in alignment\n", 
           aln.name, aln.start, aln.end, aln.score);*/
        }
      else
        {
	  found = TRUE;
	
          g_array_index(bc->alignArr, ALN, idx).score = aln.score;

          char *scoreStr = blxprintf("%.1f", aln.score);
          scoreLen = strlen(scoreStr);
          g_free(scoreStr);

          if (scoreLen > bc->maxScoreLen) 
            bc->maxScoreLen = scoreLen;
        }
    }

  fclose(file);
  
  if (found)
    {
      bc->displayScores = TRUE;
    
      if (warnings)
        g_warning("Some sequences in the scores file were not found in the alignment.\n");
    }
  else
    {
      g_critical("Error reading scores file: no sequences in the scores file were found in the alignment.\n");
    }
}


int main(int argc, char **argv)
{
  FILE    
  *file, *pipe;
  
  char    
  *scoreFile = 0,
  *readMatchFile = 0,
  *colorCodesFile = 0,
  *markupColorCodesFile = 0,
  *output_format = 0,
  *optargc;
  
  int     
  i,
//  pw, ph,	 /* pixel width and height */
  output_probs = 0,
  init_tree = 0,
  only_tree = 0,
  show_ann = 0;
  
  double   
  makeNRinit = 0.0,
  init_rmEmptyColumns = 0.0,
  init_rmGappySeqs = 0.0;
  
  
  
  /* Initialise GTK - do this early so we can show any errors in a pop-up dialog */
  gtk_init(&argc, &argv);
  
  
  /* Set up the GLib message handlers
   * 
   * There are two handlers: the default one for all non-critical messages, which will just log
   * output to the console, and one for critical messages and errors, which will display a 
   * pop-up message (the idea being that we don't bother the user unless it's something serious).
   * 
   * All errors and warnings will be sent to stderr, as will info messages (g_message_info).
   * Program output destined for stdout should use g_message (or printf is fine too).
   * g_debug will direct to stdout as well.
   * 
   * In summary:
   *   g_error: pop-up error message (always fatal)
   *   g_critical: pop-up error message
   *   g_warning: error message logged to stderr
   *   g_message_info: program info message sent to stderr
   *   g_message: program output message set to stdout
   *   (g_debug: not sure of usage scenarios but directs to stdout)
   *   
   */
  BlxMessageData msgData;
  msgData.titlePrefix = g_strdup("Belvu - ");
  msgData.parent = NULL;
  msgData.statusBar = NULL;

  g_log_set_default_handler(defaultMessageHandler, &msgData);
  g_log_set_handler(NULL, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION, 
                    popupMessageHandler, &msgData);


  /* Initialise up the context with default values */
  BelvuContext *bc = createBelvuContext();


  /* Set up tree defaults */
  char treePickString[50];
  strcpy(treePickString, SWAPstr);
  
  setTreeScaleCorr(bc, bc->treeMethod);
  
  
  char *OrganismLabel = "OS";
    
  gboolean verbose = FALSE;
  gboolean gridOn = FALSE;
  gboolean init_rmPartial = FALSE;
  
  static gboolean showHelp = FALSE;
  static gboolean showCompiled = FALSE;
  static gboolean showVersion = FALSE;
  
  /* Get the input args. We allow long args, so we need to create a long_options array */
  static struct option long_options[] =
    {
      {"compiled",		no_argument,        &showCompiled, 1},
      {"version",	        no_argument,        &showVersion, 1},

      {"help",                  no_argument,        0, 'h'},
      {0, 0, 0, 0}
    };
  
  char        *optstring="aBb:CcGghil:L:m:n:O:o:PpQ:q:RrS:s:T:t:uX:z:";
  extern int   optind;
  extern char *optarg;
  int          optionIndex; /* getopt_long stores the index into the option struct here */
  int          optc;        /* the current option gets stored here */
  
  while ((optc = getopt_long(argc, argv, optstring, long_options, &optionIndex)) != EOF)
    {
      switch (optc) 
        {
          case 0:
            /* we get here if getopt_long set a flag; nothing else to do */
            break; 
            
          case 'a': show_ann = 1;                       break;
          case 'B': bc->outputBootstrapTrees = TRUE;    break;
          case 'b': bc->treebootstraps = atoi(optarg);  break;
          case 'C': bc->saveCoordsOn = FALSE;           break;
          case 'c': verbose = TRUE;                     break;
          case 'G': bc->penalize_gaps = TRUE;           break;
          case 'g': gridOn = TRUE;                      break;
          case 'l': 
            colorCodesFile = g_malloc(strlen(optarg)+1);
            strcpy(colorCodesFile, optarg);             break;
          case 'h': showHelp = TRUE;                    break;
          case 'i': bc->ignoreGapsOn = TRUE;            break;
          case 'L': 
            markupColorCodesFile = g_malloc(strlen(optarg)+1);
            strcpy(markupColorCodesFile, optarg);       break;
          case 'm': 
            readMatchFile = g_malloc(strlen(optarg)+1);
            strcpy(readMatchFile, optarg);              break;
          case 'n':  makeNRinit = atof(optarg);         break;
          case 'O': strncpy(OrganismLabel, optarg, 2);  break;
          case 'o': 
            output_format = g_malloc(strlen(optarg)+1);
            strcpy(output_format, optarg);              break;
          case 'P': init_rmPartial = TRUE;              break;
          case 'Q': init_rmEmptyColumns = atof(optarg); break;
          case 'q': init_rmGappySeqs = atof(optarg);    break;
          case 'p': output_probs = 1;                   break;
          case 'R': bc->stripCoordTokensOn = bc->saveCoordsOn = FALSE;    break;
          case 'r': bc->IN_FORMAT = RAW;                break;
            
          case 'S': 
            switch (*optarg)
            {
              case 'a': bc->sortType = BELVU_SORT_ALPHA;    break;
              case 'o': bc->sortType = BELVU_SORT_ORGANISM; break;
              case 's': bc->sortType = BELVU_SORT_SCORE;    break;
              case 'S': bc->sortType = BELVU_SORT_SIM;      break;
              case 'i': bc->sortType = BELVU_SORT_ID;       break;
              case 'n': 
                bc->treeMethod = NJ;
                bc->sortType = BELVU_SORT_TREE;             break;
              case 'u': 
                bc->treeMethod = UPGMA; 
                bc->sortType = BELVU_SORT_TREE;             break;
              default : g_error("Illegal sorting order: %s", optarg);
            };
            break;
            
          case 's': 
            scoreFile = g_malloc(strlen(optarg)+1);
            strcpy(scoreFile, optarg);                      break;
            
          case 'T': 
          for (optargc = optarg; *optargc; optargc++) 
            {
              switch (*optargc)
                {
                  case 'n': 
                    strcpy(bc->treeMethodString, NJstr);
                    bc->treeMethod = NJ;          break;
                  case 'u': 
                    strcpy(bc->treeMethodString, UPGMAstr);
                    bc->treeMethod = UPGMA;       break;
                  case 'c':
                    bc->treeColorsOn = FALSE;         break;
                  case 'd':
                    bc->treeShowBranchlen = TRUE;    break;
                  case 'I':
                    only_tree=1;
                  case 'i':
                    init_tree = 1;            break;
                  case 'j':
                    strcpy(bc->treeDistString, JUKESCANTORstr);
                    bc->treeDistCorr = JUKESCANTOR;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'k':
                    strcpy(bc->treeDistString, KIMURAstr);
                    bc->treeDistCorr = KIMURA;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'o':
                    bc->treeCoordsOn = FALSE;         break;
                  case 'p':
                    bc->treePrintDistances = TRUE;  
                    init_tree = 1;            break;
                  case 'R':
                    bc->treeReadDistancesOn = TRUE;  break;
                  case 's':
                    strcpy(bc->treeDistString, STORMSONNstr);
                    bc->treeDistCorr = STORMSONN;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'b':
                    strcpy(bc->treeDistString, SCOREDISTstr);
                    bc->treeDistCorr = SCOREDIST;
                    setTreeScaleCorr(bc, bc->treeMethod);         break;
                  case 'r':
                    strcpy(bc->treeDistString, UNCORRstr);
                    bc->treeDistCorr = UNCORR;
                    setTreeScale(bc, 1.0);          break;
                  default : g_error("Illegal sorting order: %s", optarg);
                }
            } 
            break;
            
          case 't': strncpy(bc->Title, optarg, 255);		  break;
          case 'u': bc->displayColors = FALSE;                    break;
          case 'X': bc->mksubfamilies_cutoff = atof(optarg);      break;
          case 'z': bc->saveSeparator = *optarg;		  break;
          default : g_error("Illegal option");                    break;
        }
    }
  
  if (showVersion)
    {
      showVersionInfo();
      exit (EXIT_FAILURE);
    }

  if (showCompiled)
    {
      showCompiledInfo();
      exit (EXIT_FAILURE);
    }
  
  
  if (argc-optind < 1 || showHelp) 
    { 
      showUsageText();
      exit(1);
    }
  
  if (!strcmp(argv[optind], "-")) 
    {
      pipe = stdin;
      if (!*bc->Title) strcpy(bc->Title, "stdin");
    }
  else 
    {
      if (!(pipe = fopen(argv[optind], "r")))
	g_error("Cannot open file %s", argv[optind]);
      if (!*bc->Title) 
	strncpy(bc->Title, argv[optind], 255);
      
      bc->fileName = g_path_get_basename(argv[optind]);
      bc->dirName = g_path_get_dirname(argv[optind]);
    }
  
  int nseq = 0;
  if (bc->treeReadDistancesOn) 
    {
      /* Should this really be either or?  Problem: cannot read organism info when reading tree */
      bc->treeReadDistancesPipe = pipe;
      nseq = treeReadDistancesNames(bc);
      
      init_tree = 1;
      only_tree = 1;
      bc->treeCoordsOn = FALSE;
    }
  else
    {
      readFile(bc, pipe);
    }
  
  if (bc->organismArr->len == 0)
    suffix2organism(bc->alignArr, bc->organismArr);
  
  setOrganismColors(bc->organismArr);
  
  if (scoreFile) 
    readScores(scoreFile, bc);
  
  doSort(bc, bc->sortType);
  
  if (!bc->matchFooter && readMatchFile) 
    {
      if (!(file = fopen(readMatchFile, "r"))) 
        g_error("Cannot open file %s", readMatchFile);
      
      readMatch(bc, file);
      fclose(file);
    }
  else if (bc->matchFooter) 
    {	 
      readMatch(bc, pipe);
      fclose(pipe);
    }
  
  if (!bc->treeReadDistancesOn) 
    {
      checkAlignment(bc);
      setConsSchemeColors(bc);
    }
  
  if (verbose)
    {
      /* Print conservation statistics */
      int i, j, max, consensus = 0 ;
      double totcons = 0.0;
    
      printf("\nColumn Consensus        Identity       Conservation\n");
      printf  ("------ ---------  -------------------  ------------\n");
    
      for (i = 0; i < bc->maxLen; ++i)
        {
          max = 0;
          
          for (j = 1; j < 21; j++)
            {
              if (bc->conservCount[j][i] > max)
                {
                  max = bc->conservCount[j][i];
                  consensus = j;
                }
            }
          
          printf("%4d       %c      %4d/%-4d = %5.1f %%  %4.1f\n", 
                 i+1, b2aIndex(consensus), max, nseq, (double)max/nseq*100, bc->conservation[i]);
          totcons += bc->conservation[i];
        }
      
      printf ("\nAverage conservation = %.1f\n", totcons/(bc->maxLen*1.0));
      
      exit(0);
    }
  
  setResidueSchemeColors(bc);
  
  if (colorCodesFile) 
    {
      if (!(file = fopen(colorCodesFile, "r"))) 
        g_error("Cannot open file %s", colorCodesFile);
      
      readResidueColorScheme(bc, file, getColorArray());

      bc->residueScheme = BELVU_SCHEME_CUSTOM;
      bc->schemeType = BELVU_SCHEME_TYPE_RESIDUE;
      bc->colorByResIdOn = FALSE;
    }
  
  initMarkupColors();
  initCustomColors();
  
  if (markupColorCodesFile) 
    {
      if (!(file = fopen(markupColorCodesFile, "r"))) 
        g_error("Cannot open file %s", markupColorCodesFile);

      readResidueColorScheme(bc, file, getMarkupColorArray());

      bc->residueScheme = BELVU_SCHEME_CUSTOM;
      bc->schemeType = BELVU_SCHEME_TYPE_RESIDUE;
    }
  
  if (makeNRinit)
    mkNonRedundant(bc, makeNRinit);
  
  if (init_rmPartial)
    rmPartialSeqs(bc);
  
  if (init_rmEmptyColumns)
    rmEmptyColumns(bc, init_rmEmptyColumns/100.0);
  
  if (init_rmGappySeqs) {
    rmGappySeqs(bc, init_rmGappySeqs);
    rmFinaliseGapRemoval(bc);
  }
  
  if (output_format) 
    {
      if (!strcasecmp(output_format, "Stockholm") ||
          !strcasecmp(output_format, "Mul") ||
          !strcasecmp(output_format, "Selex"))
        {
          writeMul(bc, stdout);
        }
      else if (!strcasecmp(output_format, "MSF"))
        {
          writeMSF(bc, stdout);
        }
      else if (!strcasecmp(output_format, "FastaAlign")) 
        {
	  bc->saveFormat = BELVU_FILE_ALIGNED_FASTA;
          writeFasta(bc, stdout);
        }
      else if (!strcasecmp(output_format, "Fasta")) 
        {
	  bc->saveFormat = BELVU_FILE_UNALIGNED_FASTA;
          writeFasta(bc, stdout);
        }
      else if (!strcasecmp(output_format, "tree")) 
        {
          Tree *treeStruct = g_malloc(sizeof(Tree));
          separateMarkupLines(bc);
          treeStruct->head = treeMake(bc, TRUE);
          saveTreeNH(treeStruct->head, treeStruct->head, stdout);
          printf(";\n");
        }
      else
        {
          g_error("Illegal output format: %s", output_format);
        }

      exit(0);
    }
  
  if (bc->outputBootstrapTrees && bc->treebootstraps > 0) 
    {
      treeBootstrap(bc);
      exit(0);
    } 
  
  if (output_probs) 
    {
      outputProbs(bc, stdout);
      exit(0);
    }
  
  if (bc->mksubfamilies_cutoff) 
    {
      mksubfamilies(bc, bc->mksubfamilies_cutoff);
      exit(0);
    }

  fprintf(stderr, "\n%d sequences, max len = %d\n", bc->alignArr->len, bc->maxLen);
  
  /* Try to get 8x13 font for menus, if not set on command line */
  for ( i=0; i < argc; i++)
    {
      if (!strcmp(argv[i], "-font")) 
        break;
    }
  
  if (i == argc) 
    {
      argvAdd(&argc, &argv, "-font");
      argvAdd(&argc, &argv, "8x13");
    }
  
  // to do: get font size?
//  graphScreenSize(&screenWidth, &screenHeight, &fontwidth, &fontheight, &pw, &ph);


  if (show_ann) 
    showAnnotationWindow(bc);

  /* Calculate screen width of alignment */
  double cw = 1 + bc->maxNameLen + 1;      /* character width of initial alignment */
  
  if (bc->maxStartLen) 
    cw += bc->maxStartLen + 1;
  
  if (bc->maxEndLen)   
    cw += bc->maxEndLen + 1;
  
  if (bc->maxScoreLen) 
    cw += bc->maxScoreLen + 1;

  int scrollbarWidth = 20; /* to do: calculate scrollbar width properly */
  cw += bc->maxLen + scrollbarWidth + 2;
  
//  colorMenu = menuInitialise ("color", (MENUSPEC*)colorMENU);
//  colorEditingMenu = menuInitialise ("color", (MENUSPEC*)colorEditingMENU);
//  sortMenu = menuInitialise ("sort", (MENUSPEC*)sortMENU);
//  editMenu = menuInitialise ("edit", (MENUSPEC*)editMENU);
//  showColorMenu =  menuInitialise ("", (MENUSPEC*)showColorMENU);
//  saveMenu =  menuInitialise ("", (MENUSPEC*)saveMENU);
//  belvuMenu = menuInitialise ("belvu", (MENUSPEC*)mainMenu);
//  treeGUIMenu = menuInitialise ("", (MENUSPEC*)treeGUIMENU);
//  treeDistMenu = menuInitialise ("", (MENUSPEC*)treeDistMENU);
//  treePickMenu = menuInitialise ("", (MENUSPEC*)treePickMENU);
//  if (!displayScores) {
//    menuSetFlags(menuItem(sortMenu, "Sort by score"), MENUFLAG_DISABLED);
//    menuSetFlags(menuItem(editMenu, "Remove sequences below given score"), MENUFLAG_DISABLED);
//    menuSetFlags(menuItem(belvuMenu, "Print score and coords of line"), MENUFLAG_DISABLED);
//    menuSetFlags(menuItem(belvuMenu, "Output score and coords of line"), MENUFLAG_DISABLED);
//  }
//  menuSetFlags(menuItem(colorMenu, thresholdStr), MENUFLAG_DISABLED);


  if (init_tree)
    {
      createAndShowBelvuTree(bc);
    
      if (only_tree)
        {
          /*ACEOUT out = aceOutCreateToFile("t", "w", 0);
            graphGIF(treeGraph, out, 0);*/
          
          //graphLoop(FALSE);
          //graphFinish();
          
          exit(0);
        }
    }
  
  if (bc->outputBootstrapTrees && bc->treebootstraps < 0)
    {	
      /* Display [treebootstraps] bootstrap trees */
      bc->treebootstrapsDisplay = TRUE;
      
      bc->treebootstraps = -bc->treebootstraps;
      
      treeBootstrap(bc);
    }

  if (!colorCodesFile) 
    {
      bc->schemeType = BELVU_SCHEME_TYPE_CONS;
      bc->consScheme = BELVU_SCHEME_BLOSUM;
      setConsSchemeColors(bc);
    }
  
  /* Create the main belvu graph display of aligned sequences. */
  if (createBelvuWindow(bc, &msgData))
    {
      gtk_main();
    }
  
  return(0) ;
}
