/*  File: belvu.c
 *  Author: Erik Sonnhammer
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
 * You should have received a copy of the GNU General Public License2 * along with this program; if not, write to the Free Software
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
 * Description: Logic functions for the belvu applications
 *----------------------------------------------------------------------------
 */



/* 
   belvu.c - Pretty colour viewing of multiple alignments

-------------------------------------------------------------
|  File: belvu.c                                            |
|  Author: Erik.Sonnhammer@sbc.su.se                        |
|  Copyright (C) E Sonnhammer                               |
-------------------------------------------------------------

    Pending:

        Internal bootstrapping.
	Rationale for basic bootstrapping:
	     0. struct Subtree {nodep, seqstring}.  nodep points to node in orig tree, seqnamestring a sorted concatenation of seqence ordinal numbers
	     1. Traverse original tree from root to 2-seq nodes, insert Subtree structs at each node into a list (sorted on seqstring)
	     2. For each Boostrap tree {
	         Traverse tree, for each internal node except root: if seqstring found in Subtreelist, increment node's boostrap counter
	     }
	        Traverse bootstrap and original tree
	Rationale for group bootstrapping:
		   For each sequence, store array of sequences in merging order in bootstrap tree.
		   For each node in original tree, check if group of sequences 
		   corresponds to the first n sequences in bootstraptree array.
	Rationale for exact bootstrapping:
		   Traverse bootstrap tree from leaf to root:
		      oritree nodes = {same, different, 0}
		      if issequence(parent->other_child) {
		         if (parent->other_child == boottree(parent->other_child))
			    parent = same
			 else
			    parent = different
		      if isseen(parent->other_child) {
		         if (parent->other_child == boottree(parent->other_child))
			    parent = same
			 else
			    parent = diferent
		      if same, continue to parent

	Clean up the parsing code with str2aln calls.

	read in other trees with bootstraps
	use confidence cutoff in find orthologs
	
        make alignment collapsing easier to use, see above.

        Keyboard edit of one sequence.  How to find the right place ? !!

        Undo edits.

	Would like to abort parsing if two sequence lines accidentally have same name -
	How? - this is used in selex format...

	Koonin ideas:
	    Consensus pattern lines (Default 0% and 100% lines, maybe 90% too).
	    Add/remove any nr of lines cutoffs 0 - 100%.
	    Tool to define groups of residues and corresponding symbol.  These
	    are shown on pattern lines (priority to smaller groups).

	Low priority:

	   "Add sequences with matching segments" for more than one sequence.
	    Will probably be very deleterious to the alignment.

	    Ability to choose both foreground and background color - makes it twice
	    as slow - see "Black/white for printing.

*/


/*  Description of color_by_similarity algorithm in setConsSchemeColors():

    ( Corresponds to summing all pairwise scores)

    for each residue i {
        for each residue j {
	    if (i == j) 
	        score(i) += (count(i)-1)*count(j)*matrix(i,j)
	    else
	        score(i) += count(i)*count(j)*matrix(i,j)
	    }
	}

	if (ignore gaps)
	    n = nresidues(pos)
	else 
	    n = nsequences
		    
	if (n == 1)
	    id = 0.0
	else
	    id = score/(n*(n-1))
    }
		
*/


#include <belvuApp/belvu_.h>
#include <belvuApp/belvuWindow.h>
#include <belvuApp/belvuTree.h>
#include <belvuApp/belvuConsPlot.h>

#include <stdarg.h>
/*#include <stdlib.h> / * Needed for RAND_MAX but clashes with other stuff */
#include <sys/types.h>
#include <sys/time.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <math.h>


#define FETCH_PROG_ENV_VAR        "BELVU_FETCH"       /* environment variable used to specify the fetch program */
#define FETCH_URL_ENV_VAR         "BELVU_FETCH_WWW"   /* environment variable used to specify the WWW-fetch URL */
#define DEFAULT_FETCH_PROG        "efetch"            /* default program for fetching sequences */
#define DEFAULT_FETCH_URL         "http://www.sanger.ac.uk/cgi-bin/otter/52/nph-pfetch?request=%s"
//#define DEFAULT_FETCH_URL         "http://www.sanger.ac.uk/cgi-bin/seq-query?%s"  /* default url for fetching sequences in WWW-fetch mode */



/*  BLOSUM62 930809

#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209

  Note: to use with a2b[], always subtract 1 from the values !!!!

  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X  \* */ 
static int BLOSUM62[24][24] = {
  {4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
  {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
  {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
  {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
  {0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
  {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
  {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
  { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
  {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
  {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
  {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
  {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
  {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
  {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
  {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
  { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
  { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
  {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
  {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
  { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
  {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
  {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
  { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
  {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -1}
 };


/* ASCII-to-binary translation table
  Note: to use with BLOSUM62[], always subtract 1 from the values !!!! */
#undef NA
#define NA 0
static int a2b[] =
  {
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA, 1,NA, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
    15, 6, 2,16,17,NA,20,18,NA,19,NA,NA,NA,NA,NA,NA,
    NA, 1,NA, 5, 4, 7,14, 8, 9,10,NA,12,11,13, 3,NA,
    15, 6, 2,16,17,NA,20,18,NA,19,NA,NA,NA,NA,NA,NA,

    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA 
  };














#ifdef OLD_BELVU_CODE

/* Global variables *******************/

static int
                matchFooter = 0,
                maxfgColor = BLACK,
                midfgColor = BLACK,
                lowfgColor = BLACK,
                maxbgColor = CYAN,
                midbgColor = MIDBLUE,
                lowbgColor = LIGHTGRAY,
                oldmaxfgColor,
                oldmidfgColor,
                oldlowfgColor,
                oldmaxbgColor,
                oldmidbgColor,
                oldlowbgColor,
                printColorsOn = 0,
                wrapLinelen,
                ignoreGapsOn = 0,
                treeOrderNr,	/* Tree node number */
                conservationWindow=1,
                rmEmptyColumnsOn=1,
                saveCoordsOn=1,
                treeMethod,
/*                treeMethod = NJ, my code....should this be initialised ? */
                maxTreeWidth = 0,
                treeColorsOn = 1,
                treeShowBranchlen = 0,
                treeShowOrganism = 1,
                treeCoordsOn = 1,
                treebootstraps = 0, /* Number of bootstrap trees to be made */
	outputBootstrapTrees = 0,  /* Output the individual bootstrap trees */
                treebootstrapsDisplay = 0,  /* Display bootstrap trees on screen */
                stripCoordTokensOn = 1,
                lowercaseOn = 0,
                treePrintDistances = 0,
                treeReadDistancesON = 0,
                penalize_gaps = 0;

static double    Xpos,		/* x screen position of alignment window */
                Aspect,		/* Aspect ratio of coordinate system (width/height) */
                VSCRWID,	/* Width of Vertical scroll bar - depends on Aspect */
                HscrollStart,	/* Start of legal Horizontal scrollbar slider area */
                HscrollEnd,     /* End of -''- */
                VscrollStart,	/* as above */
                VscrollEnd,	/* as above */
                HsliderLen,     /* Length of Horizontal slider */
                VsliderLen,     /* Length of Vertical slider */
                SliderOffset,
                colorByResIdCutoff = 20.0, /* Colour by residue + id cutoff */
                lowIdCutoff = 0.4,	/* %id cutoff for lowest colour */
                midIdCutoff = 0.6,	/* %id cutoff for medium colour */
                maxIdCutoff = 0.8,	/* %id cutoff for maximum colour */
                lowSimCutoff = 0.5,	/* %id cutoff for lowest colour */
                midSimCutoff = 1.5,	/* %id cutoff for medium colour */
                maxSimCutoff = 3.0,	/* %id cutoff for maximum colour */
                oldMax,
                oldMid,
                oldLow,
                screenWidth,
                screenHeight,
                fontwidth, fontheight,
                tree_y, treeLinewidth = 0.3,
                conservationLinewidth = 0.2,
               *conservation,	        /* The max conservation in each column [0..maxLen] */
                conservationRange,
                conservationXScale=1,
                conservationYScale=2,
                oldlinew,
                treeBestBalance,
                treeBestBalance_subtrees,
                treeDistCorr,
                treePickMode,
                treeScale,
                mksubfamilies_cutoff = 0.0;

static char     belvuVersion[] = "2.35",
                line[MAXLENGTH+1],  /* General purpose array */
                Title[256] = "",    /* Window title */
                *ruler=0,	    /* Residue ruler on top */
                stats[256],	    /* Status bar string */
                dirName[DIR_BUFFER_SIZE+1], /* Default directory for file browser */
                fileName[FIL_BUFFER_SIZE+1],
                maxText[11],	/* Cutoff text arrays for boxes */
                midText[11],
                lowText[11],
                treeMethodString[50],
                treeDistString[50],
                treePickString[50],
                saveSeparator='/',
                gapChar='.',
                *OrganismLabel = "OS";

static Graph    belvuGraph = GRAPH_NULL, showColorGraph = GRAPH_NULL, treeGraph = GRAPH_NULL, 
                conservationGraph = GRAPH_NULL, saveGraph = GRAPH_NULL, annGraph = GRAPH_NULL, 
                treeGUIGraph = GRAPH_NULL, gapCharGraph = GRAPH_NULL, organismsGraph = GRAPH_NULL ;

static Stack    stack ;
static Stack    AnnStack ;
static Array    Align, MarkupAlign, organismArr, bootstrapGroups;

static ALN      aln,	     /* General purpose ALN */
               *alnp,        /* General purpose ALNP */
               *Highlight_alnp=0;

static treeStruct *treestruct;	/* General purpose treeStruct */
static treeNode   *treeHead,	/* Global current treehead */
                  *treeBestBalancedNode;
static FILE    *treeReadDistancesPipe;

ACEIN ace_in = NULL ;

static int stripCoordTokens(char *cp);
static void readColorCodesMenu(void);
static void readColorCodes(FILE *file, int *colorarr);
static void showColorCodesRedraw(void);
static void saveColorCodes(void);
static void belvuRedraw(void);
static void colorSchemeStandard(void);
static void colorByResidue(void);
static void colorSchemeGibson(void);
static void colorSchemeCys(void);
static void colorSchemeEmpty(void);
static void colorId(void);
static void colorIdSim(void);
static void colorSim(void);
static void colorRect(void);
static void xmosaicFetch(void);
static void saveAlignment(void);
static void saveRedraw(void);
static void saveMSFSelect(void);
static void saveMulSelect(void);
static void saveFastaAlnSelect(void);
static void saveFastaSelect(void);
static void treeUPGMAselect(void);
static void treeNJselect(void);
static void treeSWAPselect(void);
static void treeROOTselect(void);
static void treeUNCORRselect(void);
static void treeKIMURAselect(void);
static void treeJUKESCANTORselect(void);
static void treeSTORMSONNselect(void);
static void treeScoredistselect(void);
static void treeBootstrap(void);
static void wrapAlignmentWindow(void);
static void mkNonRedundant(double cutoff);
static void mkNonRedundantPrompt(void);
static void rmPartialSeqs(void);
static void rmGappySeqs(double cutoff);
static void rmGappySeqsPrompt(void);
static void readLabels(void);
static void selectGaps(void);
static void listIdentity(void);
static double identity(char *s1, char *s2);
static void rmPicked(void);
static void rmPickLeft(void);
static void Help(void);
static void rmColumnPrompt(void);
static void rmColumnCutoff(void);
static void rmColumnLeft(void);
static void rmColumnRight(void);
static void rmEmptyColumns(double cutoff);
static void rmEmptyColumnsInteract(void);
static void rmScore(void);
static void alphaSort(void);
static void simSort(void);
static void idSort(void);
static void organismSort(void);
static void scoreSortBatch(void);
static void scoreSort(void);
static void treeSettings(void);
static void treeDisplay(void);
static void treeSortBatch(void);
static void treeSort(void);
static void printScore(void);
static void belvuMenuDestroy(void);
static void belvuDestroy(void);
static void readMatch(FILE *fil);
int  findCommand (char *command, char **retp);
static void colorByResId(void);
static void printColors (void);
static void ignoreGaps (void);
static void finishShowColorOK(void);
static void finishShowColorCancel(void);
static void initResidueColors(void);
static void initMarkupColors(void);
static int  isAlign(char c);
static int  isGap(char c);
static void treePrintNH_init(void);
static void treeFindOrthologs(void);
static void treeSetLinewidth(char *cp);
static void treeSetScale(char *cp);
static void conservationPlot(void);
static void conservationSetWindow(void);
static void conservationSetLinewidth(void);
static void conservationSetScale(void);
static void rmEmptyColumnsToggle(void);
static void markup(void);
static void hide(void);
static void unhide(void);
static double score(char *s1, char *s2);
static void arrayOrder(void);
static void arrayOrder10(void);
static void highlightScoreSort(char mode);
static int bg2fgColor (int bkcol);
static void alncpy(ALN *dest, ALN *src);
static int alignFind( Array arr, ALN *obj, int *ip);
static void centerHighlighted(void);
static void str2aln(char *src, ALN *alnp);
static treeNode *treeReroot(treeNode *node);
static void treeSwapNode(treeNode *node);
static void treeDraw(treeNode *treeHead);
static treeNode *treecpy(treeNode *node);
static void treeTraverse(treeNode *node, void (*func)());
static void treeBalanceByWeight(treeNode *node1, treeNode *node2, double *llen, double *rlen);
static void treeReadDistances(double **pairmtx);
static double treeSize3way(treeNode *node, treeNode *fromnode);
static void showOrganisms(void);
static void parseMulLine(char *line, ALN *aln);
static void lowercase(void);
static void rmFinaliseGapRemoval(void);
static void graphButtonCheck(char* text, void (*func)(void), double x, double *y, int On);


myGraphDestroyDeclare(belvuDestroy) ;
myGraphDestroyDeclare(treeDestroy) ;

static int oldColor[256];



#define xmosaicFetchStr        "Use efetch via WWW"

static MENU belvuMenu;

static MENUOPT mainMenu[] = 
{   
  {belvuMenuDestroy,    "Quit"},
  {Help,                "Help"},
  {wrapAlignmentWindow, "Wrap-around alignment window (pretty print)"},
  {graphPrint,          "Print this window"},
  {menuSpacer,          ""},
  {treeDisplay,         "Make tree from current alignment"},
  {treeSettings,        "Tree options"},
  {menuSpacer,          ""},
  {conservationPlot,    "Plot conservation profile"},
  {menuSpacer,          ""},
  {saveAlignment,       "Save alignment"},
  {saveRedraw,          "Save alignment as..."},
  {printScore,          "Output score and coords of line"},
    /* readMatchPrompt,     "Add sequence with matching segments", */   
    /* Too difficult to use interactively - rely solely on matchFooters instead */
  {menuSpacer,          ""},
  {xmosaicFetch,        xmosaicFetchStr},
  {listIdentity,        "Compare all vs. all, list identity"},
  {graphCleanUp,        "Clean up windows"},
  {0, 0}
};

static MENUOPT treeMenu[] = 
{   
  {treeDestroy,         "Quit"},
  {graphPrint,          "Print"},
  {treeSettings,        "Tree options"},
  {treePrintNH_init,    "Save tree in New Hampshire format"},
  {treeFindOrthologs,   "Find putative orthologs"},
  {showOrganisms,       "Show current organisms"},
  {0, 0}
};

static MENUOPT conservationMenu[] = 
{   
  {graphDestroy,        "Quit"},
  {graphPrint,          "Print"},
  {conservationSetWindow, "Set window size"},
  {conservationSetScale,  "Set plot scale factors"},
  {conservationSetLinewidth,    "Set line width of tree"},
  {0, 0}
};

static MENUOPT wrapMenu[] = 
{   
  {graphDestroy,        "Quit"},
  {graphPrint,          "Print"},
  {wrapAlignmentWindow, "New Wrap-around alignment window"},    
  {0, 0}
};



#define colorByResidueStr      "Colour scheme: By residue, Last palette"
#define colorSchemeStandardStr "Colour scheme: By residue, Erik's"
#define colorSchemeGibsonStr   "Colour scheme: By residue, Toby's"
#define colorSchemeCysStr      "Colour scheme: By residue, Cys/Gly/Pro"
#define colorSchemeEmptyStr    "Colour scheme: By residue, Clean slate"
#define colorSimStr            "Colour scheme: Average similarity by Blosum62"
#define colorIdStr             "Colour scheme: Percent identity only"
#define colorIdSimStr          "Colour scheme: Percent identity + Blosum62"
#define colorRectStr           "Display colors (faster without)"
#define printColorsStr         "Use gray shades (for printing)"
#define ignoreGapsStr          "Ignore gaps in conservation calculation"
#define thresholdStr           "Only colour residues above %id threshold"
static MENU colorMenu;
static MENUOPT colorMENU[] = 
{   
  {colorSchemeStandard, colorSchemeStandardStr}, 
  {colorSchemeGibson,   colorSchemeGibsonStr},   
  {colorSchemeCys,      colorSchemeCysStr},      
  {colorSchemeEmpty,    colorSchemeEmptyStr},
  {colorByResidue,      colorByResidueStr}, 
  {colorByResId,        thresholdStr},
/*    setColorCodes,       "Set colours manually",*/
  {saveColorCodes,      "Save current colour scheme"},
  {readColorCodesMenu,  "Read colour scheme from file"},
  {menuSpacer,          ""},
  {colorSim,            colorSimStr},            
  {colorId,             colorIdStr},             
  {colorIdSim,          colorIdSimStr},
  {ignoreGaps,          ignoreGapsStr},
  {menuSpacer,          ""},
  {printColors,         printColorsStr},
  {markup,              "Exclude highlighted from calculations on/off"},
  {colorRect,           colorRectStr},
  {lowercase,           "Highlight lowercase characters"},
  {showColorCodesRedraw, "Open window to edit colour scheme"},
  {0, 0}
};

static MENU colorEditingMenu;
static MENUOPT colorEditingMENU[] = 
{   
  {showColorCodesRedraw, "Close Colour Codes window first"},
  {0, 0}
};

#define rmEmptyColumnsStr "Automatically remove all-gap columns after sequence removals"
static MENU editMenu;
static MENUOPT editMENU[] = 
{   
  {rmPicked,            "Remove highlighted line"},
  {rmPickLeft,          "Remove many sequences"},
  {rmGappySeqsPrompt,   "Remove gappy sequences"},
  {rmPartial,           "Remove partial sequences"},
  {mkNonRedundantPrompt,"Make non-redundant"},
  {rmOutliers,          "Remove outliers"},
  {rmScore,             "Remove sequences below given score"},
  {menuSpacer,          ""},
  {rmColumnPrompt,      "Remove columns"},
  {rmColumnLeft,        "<- Remove columns to the left of cursor (inclusive)"},
  {rmColumnRight,       "Remove columns to the right of cursor (inclusive) -> "},
  {rmColumnCutoff,      "Remove columns according to conservation"},
  {rmEmptyColumnsInteract, "Remove gappy columns"},
  {rmEmptyColumnsToggle, rmEmptyColumnsStr},
  {menuSpacer,          ""},
  {menuSpacer,          ""},
  {readLabels,          "Read labels of highlighted sequence and spread them"},
  {selectGaps,          "Select gap character"},
  {hide,                "Hide highlighted line"},
  {unhide,              "Unhide all hidden lines"},
  {0, 0}
};

static MENU showColorMenu;
static MENUOPT showColorMENU[] = 
{   
  {finishShowColorCancel,   "CANCEL"},
  {finishShowColorOK,       "OK"},
  {0, 0}
};

static MENU sortMenu;
static MENUOPT sortMENU[] = 
{   
  {scoreSort,           "Sort by score"},
  {alphaSort,           "Sort alphabetically"},
  {organismSort,        "Sort by swissprot organism"},
  {treeSort,            "Sort by tree order"},
  {simSort,             "Sort by similarity to highlighted sequence"},
  {idSort,              "Sort by identity to highlighted sequence"},
  {0, 0}
};

static MENU saveMenu;
static MENUOPT saveMENU[] = 
{   
  {saveMSFSelect,           MSFStr},
  {saveMulSelect,           MulStr},
  {saveFastaAlnSelect,      FastaAlnStr},
  {saveFastaSelect,         FastaStr},
  {0, 0}
};

static MENU treeGUIMenu;
static MENUOPT treeGUIMENU[] = 
{   
  {treeUPGMAselect,         UPGMAstr},
  {treeNJselect,            NJstr},
  {0, 0}
};


static MENU treePickMenu;
static MENUOPT treePickMENU[] = 
{   
  {treeSWAPselect, SWAPstr},
  {treeROOTselect, ROOTstr},
  {0, 0}
};

static MENU treeDistMenu;
static MENUOPT treeDistMENU[] = 
{   
  {treeUNCORRselect,         UNCORRstr},
  {treeJUKESCANTORselect,    JUKESCANTORstr},
  {treeKIMURAselect,         KIMURAstr},
  {treeSTORMSONNselect,      STORMSONNstr},
  {treeScoredistselect,      Scorediststr},
  {0, 0}
};

static void fatal(char *format, ...)
{
    va_list  ap;

    printf("\nFATAL ERROR: "); 

    va_start(ap, format);
    vprintf(format, ap);
    va_end(ap);

    printf("\n"); 
    exit(-1);
}

static void scrRight() { AlignXstart++; belvuRedraw(); }
static void scrLeft()  { AlignXstart--; belvuRedraw(); }
static void scrUp()    { AlignYstart--; belvuRedraw(); }
static void scrDown()  { AlignYstart++; belvuRedraw(); }
static void scrRightBig() { AlignXstart += Alignwid; belvuRedraw(); }
static void scrLeftBig()  { AlignXstart -= Alignwid; belvuRedraw(); }
static void scrUpBig()    { AlignYstart -= Alignheig; belvuRedraw(); }
static void scrDownBig()  { AlignYstart += Alignheig; belvuRedraw(); }


/* Unset the scroll bar dragging functions. */
static void unregister(void)
{
  graphRegister(LEFT_DRAG, NULL);
  graphRegister(LEFT_UP, NULL);

  graphRegister(MIDDLE_DRAG, NULL);
  graphRegister(MIDDLE_UP, NULL);

  return ;
}



static int x2col(double *x)
{
    /* Round off to half units */
    *x = floor(*x)+0.5;

    if (*x < HscrollStart - VSCRWID+1) *x = HscrollStart - VSCRWID+1;
    else if (*x > HscrollStart - VSCRWID + Alignwid)
	*x = HscrollStart - VSCRWID + Alignwid;

    return (int) (*x - (HscrollStart - VSCRWID + 0.5) + AlignXstart + 1);
}



static void showOrganisms(void)
{
    int i;

    if (!graphActivate(organismsGraph)) {
	organismsGraph = graphCreate (TEXT_FULL_SCROLL, "Organisms", 0, 0, 0.3, 0.5);

	graphTextFormat(FIXED_WIDTH);
	/* graphRegister(PICK, savePick); */
    }
    graphPop();
    graphClear();
    graphTextBounds(50, arrayMax(organismArr)+ 2);

    for (i = 0; i < arrayMax(organismArr); i++) {
	graphColor(arrp(organismArr, i, ALN)->color);
	graphText(arrp(organismArr, i, ALN)->organism, 1, 1+i);
    }
    graphRedraw();
}



/* Action to picking a box in a tree window
*/
static void treeboxPick (int box, double x, double y, int modifier_unused)
{
    treeNode *treenodePicked;

    if (!box)
      return;

    if (!graphAssFind(assVoid(1), &treestruct))
      {
	messout("Can not find treeStruct");
	return;
      };

    if (!graphAssFind(assVoid(100+box), &treenodePicked))
      {
	messout("Picked box %d not associated", box);
	return;
      }

    /* Invert box */
    graphBoxDraw(treestruct->currentPickedBox, BLACK, WHITE);
    graphBoxDraw(box, WHITE, BLACK);
    treestruct->currentPickedBox = box;

    if (box > treestruct->lastNodeBox) { /* sequence box */

	/* Highlight it in the alignment */
	bc->selectedAln = treenodePicked->aln;
	centerHighlighted();
	belvuRedraw();
    }
    else { /* Tree node box picked - swap or reroot */
	if (treePickMode == NODESWAP) {
	    treeSwapNode(treenodePicked);
	    treeDraw(treeHead);
	}
	else {
	    if (treenodePicked->parent->parent) { /* Do not reroot with the same root */
		/* treestruct->head = treecpy(treestruct->head); */
		treeDraw(treeReroot(treenodePicked));
	    }
	}
    }
}


static void keyboard(int key, int unused)
{
    switch (key)
    {
/*    case UP_KEY: scrUp();             break;
    case DOWN_KEY: scrDown();         break;
    case LEFT_KEY: scrLeft();         break;
    case RIGHT_KEY: scrRight();       break;
*/	
    case PAGE_UP_KEY: scrUpBig();     break; /* Doesn't work */
    case PAGE_DOWN_KEY: scrDownBig(); break; /* Doesn't work */

    case UP_KEY: scrUpBig();             break;
    case DOWN_KEY: scrDownBig();         break;
    case LEFT_KEY: scrLeftBig();         break;
    case RIGHT_KEY: scrRightBig();       break;
	
    case INSERT_KEY: AlignXstart = 0; belvuRedraw();      break;
    case DELETE_KEY: AlignXstart = maxLen; belvuRedraw(); break;
    case HOME_KEY:   AlignYstart = 0; belvuRedraw();      break;
    case END_KEY:    AlignYstart = nseq; belvuRedraw();   break;
    }
}


static void colorId(void)
{
    colorScheme = COLORID;
    bc->color_by_conserv = 1;
    color_by_similarity = 0;
    bc->id_blosum = FALSE;
    colorCons();
}
static void colorIdSim(void)
{
    colorScheme = COLORIDSIM;
    bc->color_by_conserv = 1;
    color_by_similarity = 0;
    bc->id_blosum = TRUE;
    colorCons();
}

static void graphScrollBar(double x1, double y1, double x2, double y2)
{
    int
	oldColor = graphColor(LIGHTBLUE);

    graphFillRectangle(x1, y1, x2, y2);
    graphColor(BLACK);
    graphRectangle(x1, y1, x2, y2);

    graphColor(oldColor);
}    


static void graphTriangle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    Array temp = arrayCreate(6, double) ;
    int oldColor = graphColor(LIGHTBLUE);
    
    array(temp, 0, double) = x1;
    array(temp, 1, double) = y1;
    array(temp, 2, double) = x2;
    array(temp, 3, double) = y2;
    array(temp, 4, double) = x3;
    array(temp, 5, double) = y3;
    graphPolygon(temp);
    arrayDestroy(temp);

    graphColor(BLACK);
    graphLine(x1, y1, x2, y2);
    graphLine(x2, y2, x3, y3);
    graphLine(x3, y3, x1, y1);

    graphColor(oldColor);
}


static void belvuMenuDestroy()
{
    if (graphActive() == belvuGraph && 
	!saved && 
	messQuery ("Alignment was modified - save ?")) 
    {
	saveMul();
    }

    belvuDestroy() ;

    return ;
}



static void lowercase(void)
{
    lowercaseOn = !lowercaseOn;
    belvuRedraw();
}
	
static void belvuRedraw(void)
{
  char  seq[MAXLINE+1], ichar[10], widBylen[99], ch[2] = " ";
  int   i, j, box, nx, ny, seqlen, statsStart, bkcol;
  double Ypos=3, Ysave, x1, x2, y1, y2;

  /* Batch usage */
  if (!graphActivate(belvuGraph))
    return ;

  graphClear();
  graphTextFormat(FIXED_WIDTH);
  graphFitBounds(&nx, &ny);


  /* Development grid lines */
  if (gridOn)
    {
      graphColor(ORANGE);

      for (i = 0; i <= ny; i++)
	graphLine(0, i, nx, i);

      for (i=0; i <= nx; i++)
	graphLine(i, 0, i, ny);

      graphColor(BLACK);
    }

  Xpos = maxNameLen+1 + 1;
  if (1 /* coordsOn */ )
    Xpos += maxStartLen+1 + maxEndLen+1;
  if (displayScores)
    Xpos += maxScoreLen+1;


  Alignwid = nx-0.5 - VSCRWID - Xpos;
  if (Alignwid > maxLen)
    {
      Alignwid = maxLen;
      AlignXstart = 0;
    }
  if (AlignXstart < 0)
    AlignXstart = 0;
  if (AlignXstart + Alignwid > maxLen)
    AlignXstart = maxLen - Alignwid;
  seqlen = (Alignwid > MAXLINE ? MAXLINE : Alignwid);


  /* RULER *******/
  {
    if (!ruler)
      {
	ruler = g_malloc(maxLen+1);
	memset(ruler, '-', maxLen);
	for (i = 10 ; i <= maxLen; i += 10)
	  {
	    sprintf(ichar, "%d", i);
	    strncpy(ruler + i - strlen(ichar), ichar, strlen(ichar));
	  }
      }

    strncpy(seq, ruler+AlignXstart, seqlen);
    seq[seqlen] = 0;
    graphText(seq, Xpos, Ypos);
    Ypos += 2;
  }
    
  /* Matching sequences ********************/
  {
    /* No, bad idea - have to rewrite all boxPick functions, and
       horizontal scrolling... */
	
    /* graphLine(0.5, Ypos-0.5, nx-0.5, Ypos-0.5);
	
    graphText("Query goes here", 1, Ypos);
    graphTextFormat (PLAIN_FORMAT);
    Ypos += 2;
    */
  }

  Alignheig = ny - 0.5 - HSCRWID - Ypos - 0.5;
  if (Alignheig > nseq)
    {
      Alignheig = nseq;
      AlignYstart = 0;
    }
  if (AlignYstart < 0)
    AlignYstart = 0;
  if (AlignYstart + Alignheig > nseq)
    AlignYstart = nseq - Alignheig;

  /* Alignment ***********/
  Ysave = Ypos;
  for (j = AlignYstart; j < AlignYstart+Alignheig && j < nseq; j++)
    {
      alnp = arrp(Align, j, ALN);
      if (alnp->hide) {
	/* Need dummy box to get box-counting right. Ugly but cheap */
	graphBoxStart();
	graphBoxEnd();
	continue;
      }

      graphBoxStart();

	for (i = AlignXstart; i < AlignXstart+Alignwid && i < maxLen; i++) {

	  if (!isGap(alnp->seq[i]) && alnp->seq[i] != ' ' && colorRectangles) {
	    bkcol = findResidueBGcolor(alnp, i);
	    graphColor(bkcol);
	    graphFillRectangle(Xpos+i-AlignXstart, Ypos, Xpos+i-AlignXstart+1, Ypos+1);
	  }
	  else 
	    bkcol = NOCOLOR;

	  /* Foreground color */
	  if (bc->color_by_conserv) {
	    graphColor(bg2fgColor(bkcol));

	    *ch = alnp->seq[i];
	    graphText(ch, Xpos+i-AlignXstart, Ypos);
	  }
	}

      graphColor(BLACK);
      if (!bc->color_by_conserv) {
	strncpy(seq, alnp->seq+AlignXstart, seqlen);
	seq[seqlen] = 0;
	graphText(seq, Xpos, Ypos);
      }
	
      graphBoxEnd(); 
      Ypos++;
    }
  graphColor(BLACK);

  /* Names, coordinates and scores *******************************/
  Ypos = Ysave;
  for (j = AlignYstart; j < AlignYstart+Alignheig && j < nseq; j++)
    {
      alnp = arrp(Align, j, ALN);
      if (alnp->hide) {
	/* Need dummy box to get box-counting right. Ugly but cheap */
	graphBoxStart();
	graphBoxEnd();
	continue;
      }

      box = graphBoxStart();
      graphText(alnp->name, 1, Ypos);
      graphBoxEnd(); 
      graphBoxDraw(box, BLACK, alnp->color);

      if (alnp->markup != GC /* && coordsOn */ ) {
	graphText(messprintf("%*d", maxStartLen, alnp->start), maxNameLen+2, Ypos);
	graphText(messprintf("%*d", maxEndLen, alnp->end), maxNameLen+maxStartLen+3, Ypos);
      }
      if (displayScores && !alnp->markup) {
	graphText(messprintf("%*.1f", maxScoreLen, alnp->score), maxNameLen+maxStartLen+maxEndLen+4, Ypos);
      }  
    
      /* Fetch command */
      /* Sort out database here since short PIR ones will be seen as
       * Swissprot otherwise - this overrides Efetch' default
       * parsing.  I do this since normally belvu has complete seq
       * names - no need to guess them.
       * /
       if (strchr(alnp->name, ':')) *db = 0;
       else if (strchr(alnp->name, '_')) strcpy(db, "SW:");
       else if (strchr(alnp->name, '.')) strcpy(db, "WP:");
       else strcpy(db, "PIR:");
       sprintf(alnp->fetch, "%s%s", db, alnp->name);
      */
      sprintf(alnp->fetch, "%s", alnp->name);

      Ypos++;
    }


  /* Scrollbars ************/

  HscrollStart = Xpos - 0.5 + VSCRWID;
  HscrollEnd = nx - 0.5 - 2*VSCRWID;
  VscrollStart = Ysave - 0.5 + HSCRWID;
  VscrollEnd = ny - 0.5 - 2*HSCRWID;

  x1 = HscrollStart + ((double)AlignXstart/maxLen)*(HscrollEnd-HscrollStart);
  x2 = HscrollStart + ((double)(AlignXstart+Alignwid)/maxLen)*(HscrollEnd-HscrollStart);
  HsliderLen = x2-x1;

  y1 = VscrollStart + ((double)AlignYstart/nseq)*(VscrollEnd-VscrollStart);
  y2 = VscrollStart + ((double)(AlignYstart+Alignheig)/nseq)*(VscrollEnd-VscrollStart);
  VsliderLen = y2-y1;

  /* Horizontal scrollbar */



  /* patch up frame of scroll boxes */
  graphRectangle(HscrollStart-VSCRWID, ny-0.5-HSCRWID,
		 HscrollEnd+VSCRWID, ny-0.5);
  graphRectangle(nx-0.5-VSCRWID, VscrollStart-HSCRWID,
		 nx-0.5, VscrollEnd+HSCRWID);


  /* Buttons *********/
  statsStart = 1;
  graphButton("File", Help, statsStart, 0.9);  /* Uses main menu, haha */
  statsStart += 6;
  editButtonBox = graphButton("Edit", rmPickLeft, statsStart, 0.9);
  statsStart += 6;
  colorButtonBox = graphButton("Colour", showColorCodesRedraw, statsStart, 0.9);
  statsStart += 8;
  sortButtonBox = graphButton("Sort", alphaSort, statsStart, 0.9);
  statsStart += 6;

  /* Status bar ************/
  sprintf(widBylen, "(%dx%d)", nseq, maxLen);
  graphText(widBylen, 1, 3);
  /* statsStart += strlen(widBylen)+1; */
  graphText("Picked:", statsStart, 1);
  statsStart += 8;
  statsBox = graphBoxStart();
  graphTextPtr(stats, statsStart, 1, nx-statsStart-1);
  graphBoxEnd();
  graphBoxDraw(statsBox, BLACK, LIGHTBLUE);

  /* Frame lines *************/
  graphLine(0.5, 0.5, nx-0.5, 0.5);
  graphLine(0.5, 0.5, 0.5, ny-0.5);
  graphLine(nx-0.5, 0.5, nx-0.5, ny-0.5);
  graphLine(0.5, ny-0.5, nx-0.5, ny-0.5);
  graphLine(0.5, 2.5, nx-0.5, 2.5);
  graphLine(0.5, Ysave-0.5, nx-0.5, Ysave-0.5);

  Xpos = 0.5 + maxNameLen + 1; 
  graphLine(Xpos, Ysave-0.5, Xpos, ny-0.5);
  Xpos += maxStartLen + 1;
  graphLine(Xpos, Ysave-0.5, Xpos, ny-0.5);
  Xpos += maxEndLen + 1;
  graphLine(Xpos, Ysave-0.5, Xpos, ny-0.5);
  if (displayScores) {
    Xpos += maxScoreLen + 1;
    graphLine(Xpos, Ysave-0.5, Xpos, ny-0.5);
  }
  

  /* Highlighted sequence **********************/
  highlight(1);

  setMenuCheckmarks();

  graphRedraw() ;

  return ;
}



static void graphButtonCheck(char* text, void (*func)(void), double x, double *y, int On)
{
    char buttont[1024];

    strcpy(buttont, messprintf("%s %s", On ? "*" : " ", text));

    graphButton(buttont, func, x, *y);
    *y += 2;
}

static void treeUPGMAselect(void){
    strcpy(treeMethodString, UPGMAstr);
    treeMethod = UPGMA;
    treeScale = 1.0;
    treeSettings();
}
static void treeNJselect(void){
    strcpy(treeMethodString, NJstr);
    treeMethod = NJ;
    if (treeDistCorr == UNCORR)
	treeScale = 1.0;
    else 
	treeScale = 0.3;
    treeSettings();
}



static void treeSWAPselect(void) {
    strcpy(treePickString, SWAPstr);
    treePickMode = NODESWAP;
    treeSettings();
}
static void    treeROOTselect(void) {
    strcpy(treePickString, ROOTstr);
    treePickMode = NODEROOT;
    treeSettings();
}


static void treeUNCORRselect(void){
    strcpy(treeDistString, UNCORRstr);
    treeDistCorr = UNCORR;
    treeScale = 1.0;
    treeSettings();
}

static void treeJUKESCANTORselect(void){
    strcpy(treeDistString, JUKESCANTORstr);
    treeDistCorr = JUKESCANTOR;
    setTreeScaleCorr(bc);
    treeSettings();
}
static void treeKIMURAselect(void){
    strcpy(treeDistString, KIMURAstr);
    treeDistCorr = KIMURA;
    setTreeScaleCorr(bc);
    treeSettings();
}
static void treeSTORMSONNselect(void){
    strcpy(treeDistString, STORMSONNstr);
    treeDistCorr = STORMSONN;
    setTreeScaleCorr(bc);
    treeSettings();
}
static void treeScoredistselect(void){
    strcpy(treeDistString, Scorediststr);
    treeDistCorr = Scoredist;
    setTreeScaleCorr(bc);
    treeSettings();
}

static void saveCoordsToggle(void) {
    saveCoordsOn = !saveCoordsOn;
    saveRedraw();
}

static void saveSeparatorToggle(void) {
    if (saveSeparator == '/')
	saveSeparator = '=';
    else 
	saveSeparator = '/';
    saveRedraw();
}

static void saveOK(void) {
    
    saveAlignment();

    graphDestroy();
    graphActivate(belvuGraph);
}

static void saveRedraw(void)
{
    double 
	x = 1.0,
	y = 1.0;

    if (!graphActivate(saveGraph)) {
	saveGraph = graphCreate (TEXT_FIT, "Save As", 0, 0, 0.5, 0.25);

	graphTextFormat(FIXED_WIDTH);
	/* graphRegister(PICK, savePick); */
    }
    graphPop();
    graphClear();

    graphText("Format:", x, y);
    graphNewBoxMenu(graphButton(saveFormat, saveRedraw, x+9, y), saveMenu);
    y += 2.0;

    graphButtonCheck("Save with coordinates", saveCoordsToggle, x, &y, saveCoordsOn);

    if (saveCoordsOn) {
	char saveSeparatort[10];

	graphText("Separator char between name and coordinates:", x, y);
	strcpy(saveSeparatort, messprintf(" %c ", saveSeparator));
	graphButton(saveSeparatort, saveSeparatorToggle, x+45, y);
	y += 1.0;
	graphText("(Use = for GCG)", x, y);
	y += 1.0;
    }

    y += 2;
    
    graphButton("OK", saveOK, x, y);
    graphButton("Cancel", graphDestroy, x+4, y);

    graphRedraw();
}


static void treeShowBranchlenToggle(void) {
    treeShowBranchlen = !treeShowBranchlen;
    treeSettings();
}


static void treeShowOrganismToggle(void) {
    treeShowOrganism = !treeShowOrganism;
    treeSettings();
}


static void treeSettings(void)
{
    double 
	x = 1.0,
	y = 1.0;
    static char
	linewidthstr[256],
	treeScalestr[256];

    if (!graphActivate(treeGUIGraph)) {
	treeGUIGraph = graphCreate (TEXT_FIT, "tree GUI", 0, 0, 0.6, 0.35);

	graphTextFormat(FIXED_WIDTH);
	/* graphRegister(PICK, savePick); */
    }
    graphPop();
    graphClear();

    graphText("Select tree building method:", x, y);
    graphNewBoxMenu(graphButton(treeMethodString, treeSettings, x+30, y), treeGUIMenu);
    y += 2.0;

    graphText("Select distance correction method:", x, y);
    graphNewBoxMenu(graphButton(treeDistString, treeSettings, x+35, y), treeDistMenu);
    y += 3.0;

    graphText("Display options:", x, y);
    y += 2.0;

    sprintf(treeScalestr, "%.2f", treeScale);
    graphText("Tree Scale:", x, y);
    graphTextEntry (treeScalestr, 5, x+12, y, (void*)(treeSetScale));
    y += 2.0;

    sprintf(linewidthstr, "%.2f", treeLinewidth);
    graphText("Line width:", x, y);
    graphTextEntry (linewidthstr, 5, x+12, y, (void*)(treeSetLinewidth));
    y += 2.0;

    graphText("Action when picking a node:", x, y);
    graphNewBoxMenu(graphButton(treePickString, treeSettings, x+30, y), treePickMenu);
    y += 2.0;

    graphButtonCheck("Display branch lengths", treeShowBranchlenToggle, x, &y, treeShowBranchlen);

    graphButtonCheck("Display organism", treeShowOrganismToggle, x, &y, treeShowOrganism);

    /*graphButtonCheck("Sample trees in intervals", treeSampleToggle, x, &y, treeBootstrapOn);
    if (treeBootstrapOn) {
              char saveSeparatort[10];

	graphText("Number of bootstrapsSeparator char between name and coordinates:", x, y);
	strcpy(saveSeparatort, messprintf(" %c ", saveSeparator));
	graphButton(saveSeparatort, saveSeparatorToggle, x+45, y);
	y += 1.0;
	graphText("(Use = for GCG)", x, y);
	y += 1.0;
    }*/

    y += 2;
    graphButton("Make tree", treeDisplay, x, y);

    y += 2;
    graphButton("Quit", graphDestroy, x, y);

    graphRedraw();
}



static void treeSetLinewidth(char *cp)
{
    treeLinewidth = atof(cp);
    treeSettings();
}


static void treeSetScale(char *cp)
{
    treeScale = atof(cp);
    treeSettings();
}


/* Note: this treecpy does not reallocate any strings, assuming they will 
   always exist and never change.  Serious problems await if this is
   not true
*/
static treeNode *treecpy(treeNode *node)
{
    treeNode *newnode;

    if (!node) return 0;

    newnode = g_malloc(sizeof(treeNode));

    newnode->dist      = node->dist;
    newnode->branchlen = node->branchlen;	
    newnode->boot      = node->boot;
    newnode->name      = node->name;
    newnode->organism  = node->organism;
    newnode->aln       = node->aln;
    newnode->box       = node->box;
    newnode->color     = node->color;

    newnode->left      = treecpy(node->left);
    newnode->right     = treecpy(node->right);

    return newnode;
}


static void simSort(void) 
{
  highlightScoreSort('P');
}

static void idSort(void) 
{
  highlightScoreSort('I');
}


void argvAdd(int *argc, char ***argv, char *s)
{
    char **v;
    int i;

    v = (char **)malloc(sizeof(char *)*(*argc+2));

    /* Copy argv */
    for (i = 0; i < (*argc); i++)
	v[i] = (*argv)[i];

    /* Add s */
    v[*argc] = (char *)malloc(strlen(s)+1);
    strcpy(v[*argc], s);

    (*argc)++;

    /* Terminator - ANSI C standard */
    v[*argc] = 0;

    /* free(*argv); */   /* Too risky */
    
    *argv = v;
}


static void finishShowColorOK(void)
{
  graphActivate(showColorGraph);
  graphDestroy();
  showColorGraph = GRAPH_NULL ;
  
  belvuRedraw();

  return ;
}

static void finishShowColorCancel(void)
{
    int i;
    
    if (bc->color_by_conserv) {
	bc->maxfgColor = oldbc->maxfgColor;
	bc->midfgColor = oldbc->midfgColor;
	bc->lowfgColor = oldbc->lowfgColor;
	bc->maxbgColor = oldbc->maxbgColor;
	bc->midbgColor = oldbc->midbgColor;
	bc->lowbgColor = oldbc->lowbgColor;

	if (color_by_similarity) {
	    bc->maxSimCutoff = oldMax;
	    bc->midSimCutoff = oldMid;
	    bc->lowSimCutoff = oldLow;
	}
	else {
	    bc->maxIdCutoff = oldMax;
	    bc->midIdCutoff = oldMid;
	    bc->lowIdCutoff = oldLow;
	}
    }
    else
	for (i=65; i<96; i++) color[i] = oldColor[i];
    
    setConsSchemeColors();

    finishShowColorOK();
}

static void setSim(char *cp) {
    bc->maxSimCutoff = atof(maxText);
    bc->midSimCutoff = atof(midText);
    bc->lowSimCutoff = atof(lowText);
    showColorCodesRedraw();
}
static void setId(char *cp) {
    bc->maxIdCutoff = atof(maxText);
    bc->midIdCutoff = atof(midText);
    bc->lowIdCutoff = atof(lowText);
    showColorCodesRedraw();
}

static void belvuConfColourMenu(KEY key, int box)
{
  int *colour;

  if (graphAssFind(assVoid(box+2000), &colour))
    { 
      *colour = key;
      *(colour+32) = key;
      graphBoxDraw(box, BLACK, *colour);
      showColorCodesRedraw();
    } 
  else
    {
      messout("Cannot find color box %d", box);
    }

  return ;
}

static void belvuConfColour(int *colour, int init, double x, double y, int width, char *text)
{ 
    int box;
    if (text)
	graphText(text, x+width+1, y);

    box = graphBoxStart();
    graphRectangle(x, y, x+width, y+1);
    graphBoxEnd();
    graphBoxFreeMenu(box, belvuConfColourMenu, graphColors);
    *colour = init;
    graphBoxDraw(box, BLACK, init);
    graphAssRemove (assVoid(box+2000)) ;
    graphAssociate(assVoid(box+2000), colour);
}

static void showConservCode(char *text, int *fgcolor, int *bgcolor, double y, char *id, int box, void (*fn)(char*))
{
    graphText(text, 1, y); 
    graphTextScrollEntry (id, 10, 5, 6, y, fn);

    belvuConfColour(fgcolor, *fgcolor, 12, y, 1, 0);
    belvuConfColour(bgcolor, *bgcolor, 13, y, 3, messprintf("%s", colorNames[*bgcolor]));
}


static void paintMessage(int *i) 
{
    /* graphText("All residues with a positive", 1, 1.5*i++);
       graphText("BLOSUM62 score to a coloured", 1, 1.5*i++);
       graphText("residue get the same colour.", 1, 1.5*i++);
       */
    
    graphText("Similar residues according", 1, (*i) * 1.5);
    graphText("to BLOSUM62 are coloured", 1, (*i) * 1.5 + 1);
    graphText("as the most conserved one.", 1, (*i) * 1.5 +2);
    (*i) += 3;
}

	
static void showColorCodesRedraw(void)
{
    int i, j, p;
    char *res = g_malloc(2);
    int box;

    *res = ' ' ;

    if (!graphActivate(showColorGraph)) {
	showColorGraph = graphCreate (TEXT_FIT, "Colour Codes", 0, 0, 0.35, 0.5);

	if (bc->color_by_conserv) {
	    oldbc->maxfgColor = bc->maxfgColor;
	    oldbc->midfgColor = bc->midfgColor;
	    oldbc->lowfgColor = bc->lowfgColor;
	    oldbc->maxbgColor = bc->maxbgColor;
	    oldbc->midbgColor = bc->midbgColor;
	    oldbc->lowbgColor = bc->lowbgColor;
	    
	    if (color_by_similarity) {
		oldMax = bc->maxSimCutoff;
		oldMid = bc->midSimCutoff;
		oldLow = bc->lowSimCutoff;
	    }
	    else {
		oldMax = bc->maxIdCutoff;
		oldMid = bc->midIdCutoff;
		oldLow = bc->lowIdCutoff;
	    }
	}
	else
	    for (i=65; i<96; i++) oldColor[i] = color[i];
    
	graphTextFormat(FIXED_WIDTH);
	graphNewMenu(showColorMenu);
	graphRegister(PICK, settingsPick);
    }
    graphPop();
    graphClear();

    if (bc->color_by_conserv) {
	i = 1; 

	if (color_by_similarity) {

	    graphText("Colouring by average", 1, 1.5*i);
	    graphText("BLOSUM62 score.", 1, 1.5*i+1);
	    i += 2;
	    
	    paintMessage(&i);
	    
	    sprintf(maxText, "%.2f", bc->maxSimCutoff);
	    showConservCode("Max:", &bc->maxfgColor, &bc->maxbgColor, 1.5*i++, maxText, maxBox, setSim);
	    sprintf(midText, "%.2f", bc->midSimCutoff);
	    showConservCode("Mid:", &bc->midfgColor, &bc->midbgColor, 1.5*i++, midText, midBox, setSim);
	    sprintf(lowText, "%.2f", bc->lowSimCutoff);
	    showConservCode("Low:", &bc->lowfgColor, &bc->lowbgColor, 1.5*i++, lowText, lowBox, setSim);
	}
	else {
	    oldMax = bc->maxIdCutoff;
	    oldMid = bc->midIdCutoff;
	    oldLow = bc->lowIdCutoff;

	    graphText("Colouring by %identity.", 1, 1.5*i++);
	    i++;

	    if (bc->consScheme == BELVU_SCHEME_ID_BLOSUM) 
	      paintMessage(&i);

	    sprintf(maxText, "%.1f", bc->maxIdCutoff);
	    showConservCode("Max:", &bc->maxfgColor, &bc->maxbgColor, 1.5*i++, maxText, maxBox, setId);
	    sprintf(midText, "%.1f", bc->midIdCutoff);
	    showConservCode("Mid:", &bc->midfgColor, &bc->midbgColor, 1.5*i++, midText, midBox, setId);
	    sprintf(lowText, "%.1f", bc->lowIdCutoff);
	    showConservCode("Low:", &bc->lowfgColor, &bc->lowbgColor, 1.5*i++, lowText, lowBox, setId);
	}
    
	i += 2;
	graphText("Note: foreground colors are only valid for different background colors.", 1, 1.5*i++);
	graphText("(Hit return after entering new cutoffs.)", 1, 1.5*i++);
	graphButton("OK", finishShowColorOK, 1.5, i*1.5+1);
	graphButton("Cancel", finishShowColorCancel, 5.5, i*1.5+1);

	setConsSchemeColors();
    }
    else {
	int used, col=1, row=1;

	for (i=1; i <= 20; i++) {
	    *res = b2a[i];
	    if (i == 11) {
		row = 1;
		col = 20;
	    }
	    graphText(res, col, row*1.5);
	    belvuConfColour(color+*res, color[(unsigned char)(*res)], col+2, row*1.5, 3, 
			    messprintf("%s", colorNames[color[(unsigned char)(*res)]]));
	    row++;
	}
	row++;
	graphText("Groups:", 1, row*1.5);
	for (i=0; i < NUM_TRUECOLORS; i++) {
	    used = 0;
	    for (j=65; j<96; j++) if (color[j] == i) used = 1;
	    if (used) {
		box = graphBoxStart();
		graphText(colorNames[i], 1, ++row*1.5);
		graphBoxEnd();
		graphBoxDraw(box, BLACK, i);
		graphText(":", 11, row*1.5);
		for (p=13, j=65; j<96; j++)
		    if (color[j] == i) {
			graphColor(BLACK);
			*res = j;
			graphText(res, p, row*1.5);
			p++;
		    }
	    }
	}

	graphButton("OK", finishShowColorOK, 1.5, ++row*1.5+1);
	graphButton("Cancel", finishShowColorCancel, 5.5, row*1.5+1);

	if (colorByResId(bc)) setConsSchemeColors();
    }
    graphColor(BLACK);
    belvuRedraw();

    graphActivate(showColorGraph);
    graphRedraw();
    /* graphLoop(TRUE);*/
}


static void colorByResId(void)
{
  ACEIN ace_in ;

  if ((ace_in = messPrompt ("Only colour residues above %identity: ", 
			     messprintf("%.1f", bc->colorByResIdCutoff), 
			    "fz", 0)))
    {
      aceInDouble(ace_in, &bc->colorByResIdCutoff);
      aceInDestroy(ace_in);
  
      color_by_similarity = 0;
  
      bc->colorByResIdOn = 1;
      setConsSchemeColors();

      belvuRedraw();
    }

  return ;
}


static void colorRes(void)
{
    menuUnsetFlags(menuItem(colorMenu, thresholdStr), MENUFLAG_DISABLED);
    menuSetFlags(menuItem(colorMenu, printColorsStr), MENUFLAG_DISABLED);
    menuSetFlags(menuItem(colorMenu, ignoreGapsStr), MENUFLAG_DISABLED);
    bc->color_by_conserv = bc->colorByResIdOn = 0;
    belvuRedraw();
}


static void colorByResidue(void)
{
    colorScheme = COLORBYRESIDUE;
    colorRes();
}


/* call an external shell command and print output in a text_scroll window
 *
 * This is a replacement for the old graph based text window, it has the advantage
 * that it uses gtk directly and provides cut/paste/scrolling but...it has the
 * disadvantage that it will use more memory as it collects all the output into
 * one string and then this is _copied_ into the text widget.
 * 
 * If this proves to be a problem I expect there is a way to feed the text to the
 * text widget a line a time. */
void externalCommand (char *command)
{
#if !defined(MACINTOSH)
  FILE *pipe ;
  char text[MAXLINE+1], *cp ;
  int line=0, len, maxlen=0;
  static Stack stack ;
  Graph old ;
  GString *str_text ;

  old = graphActive() ;

  str_text = g_string_new(NULL) ;

  stack = stackReCreate (stack, 50) ;

  pipe = popen (command, "r") ;
  while (!feof (pipe))
    { 
      if (!fgets (text, MAXLINE, pipe))
	break;

      len = strlen (text) ;
      if (len)
	{ 
	  if (text[len-1] == '\n') 
	    text[len-1] = '\0';
	  pushText (stack, text) ;
	  line++;
	  if (len > maxlen)
	    maxlen = len;
	}
    }
  pclose (pipe);

  stackCursor(stack, 0) ;

  while ((cp = stackNextText (stack)))
    {
      g_string_append_printf(str_text, "%s\n", cp) ;
    }

  gexTextEditorNew(command, str_text->str, 0,
		   NULL, NULL, NULL,
		   FALSE, FALSE, TRUE) ;

  g_string_free(str_text, TRUE) ;

  graphActivate (old) ;

#endif
  return ;
}






static void xmosaicFetch(void)
{
    xmosaic = (xmosaic ? 0 : 1);
    belvuRedraw();
}


static void wrapAlignmentWindow(void)
{
    /*  Column makeup:
	space, maxNameLen, 2xspace, maxEndLen, space, maxLen, maxEndLen
	i.e. 4 spaces
	*/

    int 
	paragraph=0, alnstart, alnend, alnlen, line=1, 
	i, j, oldpos, bkcol, empty, totCollapsed=0,
	collapseRes = 0, collapsePos = 0, collapseOn = 0;
    char 
	collapseStr[10],
	ch[2] = " ";
    static char 
	*seq=0, title[256];
    static int  
	*pos=0,			/* Current residue position of sequence j */
	linelen=80;
    double
	xpos;
    ACEIN wrap_in;

    if (!pos) strcpy(title, Title);
    if (!(wrap_in = messPrompt("Give line width and title (separated by a space):", 
		     messprintf("%d %s", linelen, title), "it", 0))) return;
    aceInInt(wrap_in, &linelen);
    strncpy(title, aceInWord(wrap_in), 255);
    title[255] = 0;
    aceInDestroy(wrap_in);
    wrap_in = NULL ;

    seq = g_malloc(linelen+1);
    wrapLinelen = maxNameLen + maxEndLen + maxEndLen + 5 +
		   (linelen > maxLen ? maxLen : linelen);

    if (!pos) pos = (int *)g_malloc(nseq*sizeof(int *));
    for (j = 0; j < nseq; j++) pos[j] = arrp(Align, j, ALN)->start;

    graphCreate(TEXT_FULL_SCROLL, Title, 0, 0, screenWidth, screenHeight);
    graphTextFormat(FIXED_WIDTH);
    graphRegister(LEFT_DOWN, leftDown);
    graphMenu(wrapMenu);

    if (*title) {
	graphText(title, 1, 1);
	line += 2;
    }

    while (paragraph*linelen +totCollapsed < maxLen)
    {
	for (j = 0; j < nseq; j++)
	{
	    alnstart = paragraph*linelen +totCollapsed; 
/*printf("alnstart= %d  totCollapsed= %d  collapsePos= %d \n", alnstart, totCollapsed, collapsePos);*/
	    alnlen = ( (paragraph+1)*linelen +totCollapsed < maxLen ? linelen : maxLen-alnstart );
	    alnend = alnstart + alnlen;

	    alnp = arrp(Align, j, ALN);
	    if (alnp->hide) continue;

	    for (empty=1, i = alnstart; i < alnend; i++) {
		if (!isGap(alnp->seq[i]) && alnp->seq[i] != ' ') 
		    empty = 0;
	    }

	    if (!empty) {

		for (collapsePos = 0, oldpos = pos[j], i = alnstart; i < alnend; i++) {	

		    xpos = maxNameLen+maxEndLen+4+i-alnstart -collapsePos;

		    if (alnp->seq[i] == '[') {
			/* Experimental - collapse block:
			   "[" -> Collapse start
			   "]" -> Collapse end

			   e.g.
			   
			   1 S[..........]FKSJFE
			   2 L[RVLKLKYKNS]KDEJHF
			   3 P[RL......DK]FKEJFJ

			          |
				  V

			   1 S[  0]FKSJFE
			   2 L[ 10]KDEJHF
			   3 P[  4]FKEJFJ
			   
			   Minimum collapse = 5 chars.
			   The system depends on absolute coherence in format, 
			   very strange results will be generated otherwise.
			   Edges need to be avoided manually.
			*/
			collapseOn = 1;
			collapsePos += 1;
			collapseRes = 0;
			continue;
		    }
		    if (alnp->seq[i] == ']') {
			collapseOn = 0;
			collapsePos -= 4;
			alnend -= 3;

			sprintf(collapseStr, "[%3d]", collapseRes);

/*printf("\n%s\n", collapseStr);*/

			graphText(collapseStr, xpos, line);
			continue;
		    }
		    if (collapseOn) {
			collapsePos++;
			if (!isGap(alnp->seq[i])) {
			    collapseRes++;
			    pos[j]++;
			}
			alnend++;
			continue;
		    }


		    if (!isGap(alnp->seq[i]) && alnp->seq[i] != ' ') {
			bkcol = findResidueBGcolor(alnp, i);
			graphColor(bkcol);
			graphFillRectangle(xpos, line, xpos+1, line+1);
			pos[j]++;
		    }
		    else
			bkcol = WHITE;


		    /* Foreground color */
		    if (colorByConservation(bc))
			graphColor(bg2fgColor(bkcol));
		    else
			graphColor(BLACK);
		    
		    *ch = alnp->seq[i];
		    graphText(ch, xpos, line);

		    graphColor(BLACK);
		}
		
		/*if (!(bc->color_by_conserv && printColorsOn)) {
		    strncpy(seq, alnp->seq+alnstart, alnend-alnstart);
		    seq[alnend-alnstart] = 0;
		    graphText(seq, maxNameLen+maxEndLen+4, line);
		}*/

		graphText(alnp->name, 1, line);

		if (!alnp->markup) {
		    graphText(messprintf("%*d", maxEndLen, oldpos), maxNameLen+3, line);
		    if (alnend == maxLen) 
			graphText(messprintf("%-d", pos[j]-1), maxNameLen+maxEndLen+alnlen+5, line);
		}

		line++;
	    }
	}
	paragraph++;
	line++;
	totCollapsed += collapsePos;
    }
    g_free(seq);
    graphTextBounds(wrapLinelen+4, line + 2);
    graphRedraw();
}


static void rmGappySeqsPrompt()
{
    static double cutoff = 50.0 ;

    if (!(ace_in = messPrompt ("Remove sequences with more (or equal) gaps than (%) :", 
			       messprintf("%.0f", cutoff), 
			       "fz", 0)))
	return;

    aceInDouble(ace_in, &cutoff);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    rmGappySeqs(cutoff);
    rmFinaliseGapRemoval();
}


static void gapCharOK (void) {

    int i,j;

    for (i = 0; i < nseq; i++) {
	alnp = &g_array_index(bc->alignArr, ALN, i);
	for (j = 0; j < maxLen; j++) {
	    if (isGap(alnp->seq[j])) alnp->seq[j] = gapChar;
	}
    }
    graphDestroy();
    belvuRedraw();
}


static void gapCharToggle(void) {
    if (gapChar == '.')
	gapChar = '-';
    else 
	gapChar = '.';
    selectGaps();
}


static void selectGaps(void) {

    char text[10];

    if (!graphActivate(gapCharGraph)) {
	gapCharGraph = graphCreate (TEXT_FIT, "Select gap char", 0, 0, 0.5, 0.25);

	graphTextFormat(FIXED_WIDTH);
    }
    graphPop();
    graphClear();

    graphText("Gap character:", 1, 2);
    strcpy(text, messprintf(" %c ", gapChar));
    graphButton(text, gapCharToggle, 16, 2);

    graphButton("OK", gapCharOK, 1, 4);
    graphButton("Cancel", graphDestroy, 1, 5.5);

    graphRedraw();
}


static void rmPicked(void)
{
    if (!bc->selectedAln) {
	messout ("Pick a sequence first!");
	return;
    }

    g_message("Removed %s/%d-%d.  ", bc->selectedAln->name, bc->selectedAln->start, bc->selectedAln->end);

    nseq--;
    arrayRemove(Align, bc->selectedAln, (void*)nrorder);
    saved = 0;
    arrayOrder();
    bc->selectedAln = 0;

    g_message("%d seqs left.\n\n", nseq);

    rmFinaliseGapRemoval();
}


static void rmPickLeftOff(void)
{
    rmPickLeftOn = 0;
}

static void rmPickLeft(void)
{
    rmPickLeftOn = 1;
    graphRegister(MESSAGE_DESTROY, rmPickLeftOff);
    graphMessage("Remove sequences by double clicking.  Remove this window when finished.");

}


static void rmEmptyColumnsToggle(void)
{
    bc->rmEmptyColumnsOn = !bc->rmEmptyColumnsOn;
    belvuRedraw();
}


static void rmEmptyColumnsInteract(void)
{
    static double 
	cutoff = 50.0;

    if (!(ace_in = messPrompt ("Remove columns with higher fraction of gaps than (%): ", 
		      messprintf("%.0f", cutoff), 
			       "fz", 0)))
	return;

    aceInDouble(ace_in, &cutoff);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    rmEmptyColumns(cutoff/100.0);
    rmFinaliseColumnRemoval();
}


static void rmColumnPrompt(void)
{
    int 
	from=1, to = maxLen;
    
    if (!(ace_in = messPrompt ("Remove columns From: and To: (inclusive)", 
		      messprintf("%d %d", from, to), 
			       "iiz", 0)))
	return;

    aceInInt(ace_in, &from);
    aceInInt(ace_in, &to);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    if (from < 1 || to > maxLen || to < from) {
	messout("Bad coordinates: %d-%d.  The range is: 1-%d", from, to, maxLen);
	return;
    }

    rmColumn(from, to);
    rmFinaliseColumnRemoval();
}


static void rmColumnLeft(void)
{
    int 
	from=1, to=selectedCol;
    
    if (!selectedCol) {
	messout("Pick a column first");
	return;
    }

    rmColumn(from, to);

    AlignXstart=0;

    rmFinaliseColumnRemoval();
}

static void rmColumnRight(void)
{
    int 
	from=selectedCol, to=maxLen;
    
    if (!selectedCol) {
	messout("Pick a column first");
	return;
    }

    rmColumn(from, to);

    rmFinaliseColumnRemoval();
}


/* Copied from master copy in dotter.c */
/* Bourne shells (popen) sometimes behave erratically - use csh instead
 */
int findCommand (char *command, char **retp)
{
#if !defined(NO_POPEN)
    FILE *pipe;
    static char retval[1025], *cp, csh[]="/bin/csh";

    /* Check if csh exists */
    if (access(csh, X_OK)) {
	messout("Could not find %s", csh);
	return 0;
    }

    if (!(pipe = (FILE *)popen(messprintf("%s -cf \"which %s\"", csh, command), "r"))) {
	return 0;
    }

    while (!feof(pipe))
	fgets(retval, 1024, pipe);
    retval[1024] = 0;
    pclose(pipe);

    if ((cp = strchr(retval, '\n')))
      *cp = 0;
    if (retp) *retp = retval;
    
    /* Check if whatever "which" returned is an existing and executable file */
    if (!access(retval, F_OK) && !access(retval, X_OK))
        return 1;
    else
	return 0;
#endif
}


static void printColors (void)
{
    static int 
	oldmaxbg, oldmidbg, oldlowbg,
	oldmaxfg, oldmidfg, oldlowfg;

    if (!printColorsOn) {
	oldmaxbg = bc->maxbgColor; bc->maxbgColor = BLACK; 
	oldmaxfg = bc->maxfgColor; bc->maxfgColor = WHITE;

	oldmidbg = bc->midbgColor; bc->midbgColor = GRAY;
	oldmidfg = bc->midfgColor; bc->midfgColor = BLACK;

	oldlowbg = bc->lowbgColor; bc->lowbgColor = LIGHTGRAY;
	oldlowfg = bc->lowfgColor; bc->lowfgColor = BLACK;
	printColorsOn = 1;
    }
    else {
	bc->maxbgColor = oldmaxbg;
	bc->maxfgColor = oldmaxfg;
	bc->midbgColor = oldmidbg;
	bc->midfgColor = oldmidfg;
	bc->lowbgColor = oldlowbg;
	bc->lowfgColor = oldlowfg;
	printColorsOn = 0;
    }
    setConsSchemeColors();
    belvuRedraw();
}


/* Not used - the idea was to find balancing point of a tree by walking to
   neighbors and checking if they are better.  However, the best balanced
   node may be beyond less balanced nodes.  (This is because the subtree weights
   are averaged branchlengths and not the sum of the total branchlengths).
*/
static void treeBalanceTreeRecurse(treeNode *node, double *bestbal, treeNode **bestNode)
{
    double newbal, lweight, rweight;

    /* Probe left subtree */

    if (node->left) {
	lweight = treeSize3way(node, node->left) - node->branchlen;
	rweight = treeSize3way(node->left, node) - node->left->branchlen;
	newbal = fabsf(lweight - rweight);
	
	/*
	printf("Subtree weights = %.1f  %.1f\n", lweight, rweight);
	printf("oldbal = %.1f   newbal = %.1f\n", *bestbal, newbal);
	*/

	if (newbal < *bestbal) { /* better balance */

	    printf("Better balance: %.1f (oldbest = %.1f)\n", newbal, *bestbal);
	    
	    *bestbal = newbal;
	    *bestNode = node;
	    
	    treeBalanceTreeRecurse(node->right, bestbal, bestNode);
	}
    }

    if (node->right) {

	lweight = treeSize3way(node, node->right) - node->branchlen;
	rweight = treeSize3way(node->right, node) - node->right->branchlen;
	newbal = fabsf(lweight - rweight);
	
	/*
	printf("Subtree weights = %.1f  %.1f\n", lweight, rweight);
	printf("oldbal = %.1f   newbal = %.1f\n", *bestbal, newbal);
	*/
	
	if (newbal < *bestbal) { /* better balance */
	    
	    printf("Better bal: %f\n", newbal);
	    
	    *bestbal = newbal;
	    *bestNode = node;
	    
	    treeBalanceTreeRecurse(node->right, bestbal, bestNode);
	}
    }
}


/* Declare common destroy funcs.... */
myGraphDestroy(belvuDestroy, belvuGraph)

myGraphDestroy(treeDestroy, treeGraph)

#endif /* old code */



/* These values define the defaults for the thresholds when coloring by
 * conservation; the first three are when coloring by %ID and the last 
 * three when coloring by similarity (i.e. BLOSUM62) */
#define DEFAULT_LOW_ID_CUTOFF           0.4
#define DEFAULT_MID_ID_CUTOFF           0.6
#define DEFAULT_MAX_ID_CUTOFF           0.8
#define DEFAULT_LOW_SIM_CUTOFF          0.5
#define DEFAULT_MID_SIM_CUTOFF          1.5
#define DEFAULT_MAX_SIM_CUTOFF          3.0

/* These values define the default color IDs for the conservation colors */
#define DEFAULT_MAX_FG_COLOR            BLACK
#define DEFAULT_MID_FG_COLOR            BLACK
#define DEFAULT_LOW_FG_COLOR            BLACK
#define DEFAULT_MAX_BG_COLOR            CYAN
#define DEFAULT_MID_BG_COLOR            MIDBLUE
#define DEFAULT_LOW_BG_COLOR            LIGHTGRAY  
#define DEFAULT_MAX_FG_PRINT_COLOR      WHITE
#define DEFAULT_MID_FG_PRINT_COLOR      BLACK
#define DEFAULT_LOW_FG_PRINT_COLOR      BLACK
#define DEFAULT_MAX_BG_PRINT_COLOR      BLACK
#define DEFAULT_MID_BG_PRINT_COLOR      GRAY
#define DEFAULT_LOW_BG_PRINT_COLOR      LIGHTGRAY  

/* Global variables */

/* Color names (must be in same order as Color enum) */
static char *colorNames[NUM_TRUECOLORS] = {
"WHITE", 
"BLACK", 
"LIGHTGRAY", 
"DARKGRAY",
"RED", 
"GREEN", 
"BLUE",
"YELLOW", 
"CYAN", 
"MAGENTA",
"LIGHTRED", 
"LIGHTGREEN", 
"LIGHTBLUE",
"DARKRED", 
"DARKGREEN", 
"DARKBLUE",
"PALERED", 
"PALEGREEN", 
"PALEBLUE",
"PALEYELLOW", 
"PALECYAN", 
"PALEMAGENTA",
"BROWN", 
"ORANGE", 
"PALEORANGE",
"PURPLE", 
"VIOLET", 
"PALEVIOLET",
"GRAY", 
"PALEGRAY",
"CERISE", 
"MIDBLUE"
};


/* Define the old-style acedb colors as hex values that can be parsed by GDK */
static char* colorTable[NUM_TRUECOLORS]= {
"#ffffff", 	   /* WHITE           */
"#000000",	   /* BLACK           */
"#c8c8c8",	   /* LIGHTGRAY       */
"#646464", 	   /* DARKGRAY        */
"#ff0000",         /* RED             */
"#00ff00",         /* GREEN           */
"#0000ff",	   /* BLUE            */
"#ffff00",	   /* YELLOW          */
"#00ffff",	   /* CYAN            */
"#ff00ff", 	   /* MAGENTA         */
"#ffa0a0",	   /* LIGHTRED        */
"#a0ffa0",	   /* LIGHTGREEN      */
"#a0c8ff",	   /* LIGHTBLUE       */
"#af0000", 	   /* DARKRED         */
"#00af00",	   /* DARKGREEN       */
"#0000af",         /* DARKBLUE        */
"#ffe6d2",	   /* PALERED         */ 
"#d2ffd2", 	   /* PALEGREEN       */
"#d2ebff",	   /* PALEBLUE        */
"#ffffc8",	   /* PALEYELLOW      */
"#c8ffff",	   /* PALECYAN        */
"#ffc8ff",	   /* PALEMAGENTA     */
"#a05000",	   /* BROWN           */
"#ff8000",	   /* ORANGE          */
"#ffdc6e",         /* PALEORANGE      */
"#c000ff",	   /* PURPLE          */
"#c9aaff",	   /* VIOLET          */
"#ebd7ff",	   /* PALEVIOLET      */
"#969696",	   /* GRAY            */
"#ebebeb",	   /* PALEGRAY        */
"#ff0080",	   /* CERISE          */
"#56b2de"	   /* MIDBLUE         */
} ;


/* File format names (must be in same order as BelvuFileFormat enum)  */
static char *fileFormatNames[BELVU_NUM_FILE_FORMATS] = {
  "Stockholm (Pfam/HMMER)",
  "MSF",
  "Aligned Fasta",
  "Unaligned Fasta"
};


/* Current residue colors */
static int color[] = {
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
        NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};


/* Saved custom residue colors */
static int customColor[] = {
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,

NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};


/* Residue colors */
static int markupColor[] = {
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,

        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,
        BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG,BG
};


static char b2a[] = ".ARNDCQEGHILKMFPSTWYV" ;

static int a2b_sean[] =
  {
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA, 1,NA, 2, 3, 4, 5, 6, 7, 8,NA, 9,10,11,12,NA,
    13,14,15,16,17,NA,18,19,NA,20,NA,NA,NA,NA,NA,NA,
    NA, 1,NA, 2, 3, 4, 5, 6, 7, 8,NA, 9,10,11,12,NA,
    13,14,15,16,17,NA,18,19,NA,20,NA,NA,NA,NA,NA,NA,

    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
  };


 

/* Local function declarations */
static void		   alncpy(ALN *dest, ALN *src);
static double		   score(char *s1, char *s2, const gboolean penalize_gaps);
static void		   initConservMtx(BelvuContext *bc);
static int		   countResidueFreqs(BelvuContext *bc);
static int                 stripCoordTokens(char *cp, BelvuContext *bc);
int*                       getConsColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground);
static void                rmFinalise(BelvuContext *bc) ;


/***********************************************************
 *		          Sorting			   *
 ***********************************************************/

/* Sort comparision function to sort ALNs by name */
gint alphaorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = (const ALN*)xIn;
  const ALN *y = (const ALN*)yIn;

  int retval=0;
  
  if (!(retval = strcmp(x->name, y->name)))
    {
      if (x->start == y->start)
	{
	  if (x->end < y->end)
	    retval = -1;
	  else if (x->end > y->end)
	    retval = 1;
	}
      else if (x->start < y->start)
	retval = -1;
      else if (x->start > y->start)
	retval = 1;
    }
  
  /* printf("Comparing %10s %4d %4d with %10s %4d %4d = %d\n", 
   x->name, x->start, x->end, y->name, y->start, y->end, retval); */
  
  return retval;
}


/* Sort comparision function to sort ALNs by organism */
gint organism_order(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = (const ALN*)xIn;
  const ALN *y = (const ALN*)yIn;

  if (!x->organism && !y->organism)
    return 0;
  else if (!x->organism)
    return -1;
  else if (!y->organism)
    return 1;
  else
    return strcmp(x->organism, y->organism);
}


/* Sort comparision function to sort ALNs by organism */
gint organismorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = (const ALN*)xIn;
  const ALN *y = (const ALN*)yIn;

  int retval=0;
  char *p1 = strchr(x->name, '_'), 
       *p2 = strchr(y->name, '_');

  if (!p1 && !p2) return alphaorder(x, y);
  if (!p1) return 1;
  if (!p2) return -1;
  
  if (!(retval = strcmp(p1, p2))) 
    return alphaorder(x, y);

  return retval;
}


/* Sort comparison function to sort ALNs by score */
gint scoreorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = (const ALN*)xIn;
  const ALN *y = (const ALN*)yIn;

  if (x->score < y->score)
    return -1;
  else if (x->score > y->score)
    return 1;
  else return 0;
}

static gint scoreorderRev(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = (const ALN*)xIn;
  const ALN *y = (const ALN*)yIn;

    if (x->score > y->score)
	return -1;
    else if (x->score < y->score)
	return 1;
    else return 0;
}


/* Sort comparison function to sort ALNs by nr */
gint nrorder(gconstpointer xIn, gconstpointer yIn)
{
  const ALN *x = (const ALN*)xIn;
  const ALN *y = (const ALN*)yIn;

  if (x->nr < y->nr)
    return -1;
  else if (x->nr > y->nr)
    return 1;
  else return 0;
}


void scoreSort(BelvuContext *bc)
{
  if (!bc->displayScores) 
    { 
      g_critical("No scores available.\n");
      return;
    }

  ALN aln;
  initAln(&aln);

  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  g_array_sort(bc->alignArr, scoreorder);
  
  if (bc->selectedAln) 
    {
      int ip = 0;
    
      if (!alignFind(bc->alignArr, &aln, &ip)) 
	{
	  g_critical("Program error: cannot find back highlighted sequence after sort.\n");
	  bc->selectedAln = NULL;
	}
      else
	{
	  bc->selectedAln = &g_array_index(bc->alignArr, ALN, ip);
	}
  }
  
  arrayOrder(bc->alignArr);
}


static void alphaSort(BelvuContext *bc)
{
  ALN aln;
  initAln(&aln);
  
  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  g_array_sort(bc->alignArr, alphaorder);
  
  if (bc->selectedAln) 
    { 
      int ip = 0;
      if (!alignFind(bc->alignArr, &aln, &ip)) 
	{
	  g_critical("Program error: cannot find back highlighted seq after sort.\n");
	  bc->selectedAln = NULL;
	}
      else
	{
	  bc->selectedAln = &g_array_index(bc->alignArr, ALN, ip);
	}
    }
  
  arrayOrder(bc->alignArr);
}


static void organismSort(BelvuContext *bc)
{
  ALN aln;
  initAln(&aln);

  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  g_array_sort(bc->alignArr, organismorder);
  
  if (bc->selectedAln) 
    {
      int ip = 0;
      if (!alignFind(bc->alignArr, &aln, &ip)) 
	{
	  g_critical("Program error: cannot find back highlighted seq after sort.\n");
	  bc->selectedAln = NULL;
	}
      else
	{
	  bc->selectedAln = &g_array_index(bc->alignArr, ALN, ip);
	}
    }
  
  arrayOrder(bc->alignArr);
}


void highlightScoreSort(char mode, BelvuContext *bc)
{
  if (!bc->selectedAln) 
    {
      g_critical("Please highlight a sequence first");
      return;
    }
  
  if (bc->selectedAln->markup) 
    {
      g_critical("Please do not highlight a markup line");
      return;
    }
  
  ALN aln;
  initAln(&aln);
  alncpy(&aln, bc->selectedAln);
  
  separateMarkupLines(bc);
  
  /* if (displayScores) {
   if (!(graphQuery("This will erase the current scores, do you want to continue?")))
   return;
   }*/
  
  bc->displayScores = TRUE;
  
  /* Calculate score relative to highlighted sequence */
  int i = 0;
  
  for (i = 0; i < bc->alignArr->len; ++i)
    {
      if (mode == 'P')
	{
	  ALN *curAln = &g_array_index(bc->alignArr, ALN, i);
	  curAln->score = score(alnGetSeq(&aln), alnGetSeq(curAln), bc->penalize_gaps);
	}
      else if (mode == 'I')
	{
	  ALN *curAln = &g_array_index(bc->alignArr, ALN, i);
	  curAln->score = identity(alnGetSeq(&aln), alnGetSeq(curAln), bc->penalize_gaps);
	}
      
      char *scoreStr = blxprintf("%.1f", g_array_index(bc->alignArr,ALN, i).score);
      int len = strlen(scoreStr);
      
      if (len > bc->maxScoreLen)
	bc->maxScoreLen = len;
      
      g_free(scoreStr);
    }
  
  g_array_sort(bc->alignArr, scoreorderRev);
  
  reInsertMarkupLines(bc);
  
  if (!alignFind(bc->alignArr, &aln, &i)) 
    {
      g_critical("Program error: cannot find back highlighted seq after sort.\n");
      bc->selectedAln = NULL;
    }
  else
    {
      bc->selectedAln = &g_array_index(bc->alignArr, ALN, i);
    }
  
  arrayOrder(bc->alignArr);
  bc->alignYStart = 0;
}



void doSort(BelvuContext *bc, const BelvuSortType sortType)
{
  g_array_sort(bc->alignArr, nrorder);
  
  switch(sortType) 
  {
    case BELVU_SORT_ALPHA :	  alphaSort(bc);                        break;
    case BELVU_SORT_ORGANISM :	  organismSort(bc);                     break;
    case BELVU_SORT_SCORE :	  scoreSort(bc);                        break;
    case BELVU_SORT_TREE  :	  treeSort(bc);                         break;
    case BELVU_SORT_CONS  :	  /* sort by nrorder - already done */  break; 
    
    case BELVU_SORT_SIM :
      bc->selectedAln = &g_array_index(bc->alignArr, ALN, 0);
      highlightScoreSort('P', bc); 
      break;
    
    case BELVU_SORT_ID : 
      bc->selectedAln = &g_array_index(bc->alignArr, ALN, 0);
      highlightScoreSort('I', bc); break;
    
    case BELVU_UNSORTED : break;
    
    default: 
      g_warning("Initial sort order '%d' not recognised.\n", sortType);
      break;
  }
}

/***********************************************************
 *		          Trees				   *
 ***********************************************************/

void setTreeScaleCorr(BelvuContext *bc, const int treeMethod) 
{
  if (treeMethod == UPGMA)
      bc->treeScale = 1.0;
  else if (treeMethod == NJ)
      bc->treeScale = 0.3;
}


void setTreeScale(BelvuContext *bc, const double newScale) 
{
  bc->treeScale = newScale;
}


static int treeOrder(TreeNode *node, const int treeOrderNrIn) 
{
  int treeOrderNr = treeOrderNrIn;
  
  if (node) 
    {
      treeOrderNr = treeOrder(node->left, treeOrderNr);
      
      if (node->aln)
        node->aln->nr = treeOrderNr++;
      
      treeOrderNr = treeOrder(node->right, treeOrderNr);    
    }

  return treeOrderNr;
}


void treeSortBatch(BelvuContext *bc)
{
  separateMarkupLines(bc);
  
  if (!bc->treeHead)
    {
      TreeNode *treeHead = treeMake(bc, FALSE);
      setTreeHead(bc, treeHead);
    }
  
  treeOrder(bc->treeHead, 1); /* Set nr field according to tree order */
  
  g_array_sort(bc->alignArr, nrorder);
  
  reInsertMarkupLines(bc);
}


void treeTraverse(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *treeNode)) 
{
  if (!node) 
    return;
  
  treeTraverse(bc, node->left, func);
  func(bc, node);
  treeTraverse(bc, node->right, func);
}


/* General purpose routine to convert a string to ALN struct.
   Note: only fields Name, Start, End are filled!
 */
static void str2aln(BelvuContext *bc, char *src, ALN *alnp) 
{
  char *tmp = g_strdup(src);
  stripCoordTokens(tmp, bc);

  if (sscanf(tmp, "%s%d%d", alnp->name, &alnp->start, &alnp->end) != 3)
    {
      g_critical("Name to field conversion failed for %s (%s).\n", src, tmp);
      return;
    }
  
  if (strlen(alnp->name) > MAXNAMESIZE)
    g_error("buffer overrun in %s !!!!!!!!!\n", "str2aln") ;
  
  g_free(tmp);
}


/* Find back the ALN pointers lost in treeSortBatch */
static void treeFindAln(BelvuContext *bc, TreeNode *node) 
{
  if (!node->name) 
    return;

  ALN aln;
  initAln(&aln);
  
  str2aln(bc, node->name, &aln);

  int ip = 0;
  
  if (!alignFind(bc->alignArr, &aln, &ip)) 
    {
      g_critical("Program error: cannot find back seq after sort.\n");
    }
  else
    {
      node->aln = &g_array_index(bc->alignArr, ALN, ip);
    }
}


void treeSort(BelvuContext *bc)
{
  ALN aln;
  initAln(&aln);
  
  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  treeSortBatch(bc);
  
  /* After treeSortBatch the treeNodes point to the wrong ALN,
   because arraySort only copies contents and not entire
   structs. This can be fixed by either a linear string search or
   by remaking the tree... */
  treeTraverse(bc, bc->treeHead, treeFindAln);
  
  if (bc->selectedAln) 
    {
      int ip = 0;
      if (!alignFind(bc->alignArr, &aln, &ip)) 
        {
          g_critical("Program error: cannot find back highlighted seq after sort.\n");
          bc->selectedAln = 0;
        }
      else
        {
          bc->selectedAln = &g_array_index(bc->alignArr, ALN, ip);
        }
    }
  
  /* Show the tree window (create it if necessary) */
  if (bc->belvuTree)
    gtk_window_present(GTK_WINDOW(bc->belvuTree));
  else
    createBelvuTreeWindow(bc, bc->treeHead);
}


static void treeTraverseLRfirst(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *node)) 
{
  if (!node) 
    return;
  
  treeTraverseLRfirst(bc, node->left, func);
  treeTraverseLRfirst(bc, node->right, func);
  func(bc, node);
}


static void subfamilyTrav(BelvuContext *bc, TreeNode *node) 
{
  static double dist = 0.0;
  static int 
  newgroup = 1,
  groupnr = 0;
  
  if (!node) 
    return;
  
  if (node->name) 
    {
    dist = node->branchlen;
    
    if (newgroup) 
      {
      printf("\nGroup nr %d:\n", ++groupnr);
      newgroup = 0;
      }
    
    printf("%s\n", node->name);
    }
  else
    {
    /* internal node */
    dist += node->branchlen;
    }
  
  if ( bc->mksubfamilies_cutoff > (100.0-dist) )
    { 
      newgroup = 1; 
    }
  
  /* printf("abs=%.1f  branch=%.1f\n", dist, node->branchlen); */
}


void mksubfamilies(BelvuContext *bc, double cutoff)
{
  Tree *treeStruct = g_malloc(sizeof(Tree));
  
  separateMarkupLines(bc);
  
  strcpy(bc->treeMethodString, UPGMAstr);
  bc->treeMethod = UPGMA;
  
  treeStruct->head = treeMake(bc, FALSE);
  
  treeTraverseLRfirst(bc, treeStruct->head, subfamilyTrav);
}


/***********************************************************
 *		          Alignments			   *
 ***********************************************************/

static void readFastaAlnFinalise(BelvuContext *bc, ALN *aln, char *seq)
{
  aln->sequenceStr = g_string_new(seq);
  
  if (bc->maxLen) 
    {
      if (alnGetSeqLen(aln) != bc->maxLen) 
        g_error("Differing sequence lengths: %d %d", bc->maxLen, alnGetSeqLen(aln));
    }
  else
    {
      bc->maxLen = alnGetSeqLen(aln);
    }

  int ip = 0;
  if (arrayFind(bc->alignArr, &aln, &ip, (void*)alphaorder))
    {
      g_error("Sequence name occurs more than once: %s%c%d-%d", 
              aln->name, bc->saveSeparator, aln->start, aln->end);
    }
  else
    {
      aln->nr = bc->alignArr->len + 1;
      g_array_append_val(bc->alignArr, *aln);
      g_array_sort(bc->alignArr, alphaorder);
    }
}


/* initialise an ALN with empty values */
void initAln(ALN *alnp)
{
  alnp->name[0] = '\0';
  alnp->start = 0;
  alnp->end = 0;
  alnp->sequenceStr = NULL;
  alnp->nr = 0;
  alnp->fetch[0] = '\0';
  alnp->score = 0.0;
  alnp->color = 0;
  alnp->markup = 0;
  alnp->hide = FALSE;
  alnp->nocolor = FALSE;
  alnp->organism = NULL; 
  alnp->startColIdx = 0;
}

/* Read in fasta sequences from a file and create a sequence in the given
 * alignments array for each of the fasta sequences. */
static void readFastaAln(BelvuContext *bc, FILE *pipe)
{
  char line[MAXLENGTH+1];
  GString *currentSeq = NULL;

  ALN currentAln;
  initAln(&currentAln);

  while (!feof (pipe))
    { 
      if (!fgets (line, MAXLENGTH, pipe))
	break;

      char *cp = strchr(line, '\n');
      
      /* Cut off the newline char at the end, if it has one */
      if (cp)
	*cp = 0;

      if (*line == '>')
        {
          if (currentSeq)
            {
              /* Finish off the previous sequence */
              readFastaAlnFinalise(bc, &currentAln, currentSeq->str);
              g_string_free(currentSeq, TRUE);
              currentSeq = NULL;
            }

          /* Start a new sequence */
          currentSeq = g_string_new(NULL);
          
          /* Parse the new line. Note that this resets the ALN struct for us. */
          parseMulLine(bc, line + 1, &currentAln);
        }
      else if (currentSeq)
        {
          /* Part-way through reading a sequnce; append the current line to it */
          g_string_append(currentSeq, line);
        }
    }
  
  if (currentSeq)
    {
      readFastaAlnFinalise(bc, &currentAln, currentSeq->str);
      g_string_free(currentSeq, TRUE);
      currentSeq = NULL;
    }

  bc->saveFormat = BELVU_FILE_ALIGNED_FASTA;
  
  return ;
}


static void alncpy(ALN *dest, ALN *src)
{
  strncpy(dest->name, src->name, MAXNAMESIZE);
  dest->start = src->start;
  dest->end = src->end;
  /* dest->sequenceStr = src->sequenceStr ? g_string_new(src->sequenceStr->str) : NULL; */
  dest->sequenceStr = src->sequenceStr;
  dest->nr = src->nr;			
  strncpy(dest->fetch, src->fetch, MAXNAMESIZE+10);
  dest->score = src->score;
  dest->color = src->color;
}


/***********************************************************
 *		          Parsing			   *
 ***********************************************************/

/* Convenience routine for converting "name/start-end" to "name start end".
 Used by parsing routines.
 
 Return 1 if both tokens were found, otherwise 0.
 */
static int stripCoordTokens(char *cp, BelvuContext *bc)
{
  char *end = cp;
  
  while (*end && !isspace(*end)) end++;
  
  if ((cp = strchr(cp, bc->saveSeparator)) && cp < end) {
    *cp = ' ';
    if ((cp = strchr(cp, '-')) && cp < end) {
      *cp = ' ';
      bc->IN_FORMAT = MUL;
      return 1;
    }
  }

  bc->IN_FORMAT = RAW;
  return 0;
}


/*
 Parse name, start and end of a Mul format line
 
 Convenience routine, part of readMul and other parsers
 */
void parseMulLine(BelvuContext *bc, char *line, ALN *aln)
{
  char line2[MAXLENGTH+1], *cp=line2, *cq, GRfeat[MAXNAMESIZE+1];
  
  strncpy(cp, line, MAXLENGTH);
  
  initAln(aln);
  
  if (!strncmp(cp, "#=GC", 4))
    {
      aln->markup = GC;
      cp += 5;
    }

  if (!strncmp(cp, "#=RF", 4)) 
    {
      aln->markup = GC;
    } 

  if (!strncmp(cp, "#=GR", 4)) 
    {
      aln->markup = GR;
      cp += 5;
    }
  
  if (bc->stripCoordTokensOn) 
    stripCoordTokens(cp, bc);
  
  /* Name */
  strncpy(aln->name, cp, MAXNAMESIZE);
  aln->name[MAXNAMESIZE] = 0;

  if ((cq = strchr(aln->name, ' '))) 
    *cq = 0;
  
  /* Add Start and End coords */
  if (bc->stripCoordTokensOn && (bc->IN_FORMAT != RAW) )
    sscanf(strchr(cp, ' ') + 1, "%d%d", &aln->start, &aln->end);
  
  if (aln->markup == GR) 
    {
      /* Add GR markup names to name 
       
         #=GR O83071 192 246 SA 999887756453524252..55152525....36463774777.....948472782969685958
         ->name = O83071_SA
      */
      if (bc->IN_FORMAT == MUL) 
        {
          sscanf(cp + strlen(aln->name), "%d%d%s", &aln->start, &aln->end, GRfeat);
        }
      else
        {
          sscanf(cp + strlen(aln->name), "%s", GRfeat);
        }
      
      if (strlen(aln->name)+strlen(GRfeat)+2 > MAXNAMESIZE)
        g_critical("Too long name or/and feature name");

      strcat(aln->name, " ");
      strncat(aln->name, GRfeat, MAXNAMESIZE-strlen(aln->name));
      
      /* printf("%s, %d chars\n", aln->name, strlen(aln->name)); fflush(stdout); */
    }
}


/***********************************************************
 *		          Colours			   *
 ***********************************************************/

/* Set the default colors of organisms to something somewhat intelligent */
void setOrganismColors(GArray *organismArr) 
{
    static int treeColors[16] = {
      RED, 
      BLUE,
      DARKGREEN, 
      ORANGE, 
      MAGENTA,
      BROWN, 
      PURPLE, 
      CYAN, 
      VIOLET, 
      MIDBLUE,
      CERISE, 
      LIGHTBLUE,
      DARKRED, 
      GREEN, 
      DARKBLUE,
      GRAY
  };
  
  int i = 0;
  int color = 0;
  
  for (i = 0; i < organismArr->len; ++i) 
    {
      color = treeColors[i % 16];
      /*if (i > 15) color = treeColors[15];*/
     
      g_array_index(organismArr, ALN, i).color = color;
    }
}



/* Return the markup color for the given char */
int getMarkupColor(const char inputChar)
{
  return markupColor[(unsigned char)(inputChar)];
}

/* Return the conservation color for the given char at the given index */
int getConservColor(BelvuContext *bc, const char inputChar, const int idx)
{
  return bc->colorMap[a2b[(unsigned char)(inputChar)]][idx];
}

/* Return the color from the colors array for the given char */
int getColor(const char inputChar)
{
  return color[(unsigned char)(inputChar)];
}

/* Set the color in the colors array for the given char */
void setColor(const char inputChar, const int colorNum)
{
  color[(unsigned char)(toupper(inputChar))] = colorNum;
  color[(unsigned char)(tolower(inputChar))] = colorNum;
}

/* Returns the 'color' array */
int* getColorArray()
{
  return color;
}

/* Returns the 'markupColor' array */
int* getMarkupColorArray()
{
  return markupColor;
}


/* Convert one of the old acedb-style color numbers to a hex string */
static const char* convertColorNumToStr(const int colorNum)
{
  const char *result = colorTable[WHITE];
  
  if (colorNum >= 0 && colorNum < NUM_TRUECOLORS)
    result = colorTable[colorNum];
  
  return result;
}


/* Convert an old-style ACEDB color number to a GdkColor */
void convertColorNumToGdkColor(const int colorNum, const gboolean isSelected, GdkColor *result)
{
  const char *colorStr = convertColorNumToStr(colorNum);
  getColorFromString(colorStr, result, NULL);
  
  /* If an item is selected, we use a slightly different shade of the same color. */
  if (isSelected)
    getSelectionColor(result, result);
}


//static void colorCons(BelvuContext *bc)
//{
//  setConsSchemeColors(bc);
//  
//  //  menuSetFlags(menuItem(colorMenu, thresholdStr), MENUFLAG_DISABLED);
//  //  menuUnsetFlags(menuItem(colorMenu, printColorsStr), MENUFLAG_DISABLED);
//  //  menuUnsetFlags(menuItem(colorMenu, ignoreGapsStr), MENUFLAG_DISABLED);
//  
//  bc->colorByResIdOn = FALSE;
//  //  belvuRedraw();
//}


/* Returns true if we're coloring by conservation */
gboolean colorByConservation(BelvuContext *bc)
{
  return (bc->schemeType == BELVU_SCHEME_TYPE_CONS);
}

/* Returns true if we're coloring by residue */
gboolean colorByResidue(BelvuContext *bc)
{
  return (bc->schemeType == BELVU_SCHEME_TYPE_RESIDUE);
}

/* Returns true if we're coloring by similarity */
gboolean colorBySimilarity(BelvuContext *bc)
{
  return (colorByConservation(bc) && bc->consScheme == BELVU_SCHEME_BLOSUM);
}


/* Returns true if we should only color residues above a set threshold */
gboolean colorByResId(BelvuContext *bc)
{
  /* This only applies in color-by-residue mode */
  return (colorByResidue(bc) && bc->colorByResIdOn);
}


/* Utility to clear all residue colors to white */
static void clearResidueColors(BelvuContext *bc)
{
  color['A'] = color['a'] = WHITE;
  color['B'] = color['b'] = NOCOLOR;
  color['C'] = color['c'] = WHITE;
  color['D'] = color['d'] = WHITE;
  color['E'] = color['e'] = WHITE;
  color['F'] = color['f'] = WHITE;
  color['G'] = color['g'] = WHITE;
  color['H'] = color['h'] = WHITE;
  color['I'] = color['i'] = WHITE;
  color['J'] = color['j'] = WHITE;
  color['K'] = color['k'] = WHITE;
  color['L'] = color['l'] = WHITE;
  color['M'] = color['m'] = WHITE;
  color['N'] = color['n'] = WHITE;
  color['O'] = color['o'] = WHITE;
  color['P'] = color['p'] = WHITE;
  color['Q'] = color['q'] = WHITE;
  color['R'] = color['r'] = WHITE;
  color['S'] = color['s'] = WHITE;
  color['T'] = color['t'] = WHITE;
  color['V'] = color['v'] = WHITE;
  color['U'] = color['u'] = WHITE;
  color['W'] = color['w'] = WHITE;
  color['X'] = color['x'] = WHITE;
  color['Y'] = color['y'] = WHITE;
  color['Z'] = color['z'] = NOCOLOR;
}


/* Utility to clear all custom colors to white. Should be called during initialisation. */
void initCustomColors()
{
  customColor['A'] = customColor['a'] = WHITE;
  customColor['B'] = customColor['b'] = NOCOLOR;
  customColor['C'] = customColor['c'] = WHITE;
  customColor['D'] = customColor['d'] = WHITE;
  customColor['E'] = customColor['e'] = WHITE;
  customColor['F'] = customColor['f'] = WHITE;
  customColor['G'] = customColor['g'] = WHITE;
  customColor['H'] = customColor['h'] = WHITE;
  customColor['I'] = customColor['i'] = WHITE;
  customColor['J'] = customColor['j'] = WHITE;
  customColor['K'] = customColor['k'] = WHITE;
  customColor['L'] = customColor['l'] = WHITE;
  customColor['M'] = customColor['m'] = WHITE;
  customColor['N'] = customColor['n'] = WHITE;
  customColor['O'] = customColor['o'] = WHITE;
  customColor['P'] = customColor['p'] = WHITE;
  customColor['Q'] = customColor['q'] = WHITE;
  customColor['R'] = customColor['r'] = WHITE;
  customColor['S'] = customColor['s'] = WHITE;
  customColor['T'] = customColor['t'] = WHITE;
  customColor['V'] = customColor['v'] = WHITE;
  customColor['U'] = customColor['u'] = WHITE;
  customColor['W'] = customColor['w'] = WHITE;
  customColor['X'] = customColor['x'] = WHITE;
  customColor['Y'] = customColor['y'] = WHITE;
  customColor['Z'] = customColor['z'] = NOCOLOR;
}


static void colorSchemeCGP(BelvuContext *bc)
{
  clearResidueColors(bc);
  
  color['C'] = color['c'] = CYAN;
  color['G'] = color['g'] = RED;
  color['P'] = color['p'] = GREEN;
}


static void colorSchemeCGPH(BelvuContext *bc)
{
  clearResidueColors(bc);
  
  color['C'] = color['c'] = CYAN;
  color['G'] = color['g'] = RED;
  color['P'] = color['p'] = GREEN;
  color['H'] = color['p'] = YELLOW;
}


static void colorSchemeEmpty(BelvuContext *bc)
{
  clearResidueColors(bc);
}


static void colorSchemeErik(BelvuContext *bc)
{
  /* Erik's favorite colours:
   
   C        - MIDBLUE
   GP       - CYAN
   HKR      - GREEN
   AFILMVWY - YELLOW
   BDENQSTZ - LIGHTRED
   */
  
  color['A'] = color['a'] = YELLOW;
  color['B'] = color['b'] = NOCOLOR;
  color['C'] = color['c'] = MIDBLUE;
  color['D'] = color['d'] = LIGHTRED;
  color['E'] = color['e'] = LIGHTRED;
  color['F'] = color['f'] = YELLOW;
  color['G'] = color['g'] = CYAN;
  color['H'] = color['h'] = GREEN;
  color['I'] = color['i'] = YELLOW;
  color['K'] = color['k'] = GREEN;
  color['L'] = color['l'] = YELLOW;
  color['M'] = color['m'] = YELLOW;
  color['N'] = color['n'] = LIGHTRED;
  color['P'] = color['p'] = CYAN;
  color['Q'] = color['q'] = LIGHTRED;
  color['R'] = color['r'] = GREEN;
  color['S'] = color['s'] = LIGHTRED;
  color['T'] = color['t'] = LIGHTRED;
  color['V'] = color['v'] = YELLOW;
  color['W'] = color['w'] = YELLOW;
  color['Y'] = color['y'] = YELLOW;
  color['Z'] = color['z'] = NOCOLOR;
}


static void colorSchemeGibson(BelvuContext *bc)
{
  /* Colour scheme by Gibson et. al (1994) TIBS 19:349-353
   
   Listed in Figure 1:
   
   
   Gibson      AA        Here
   
   orange      G         ORANGE (16-colours: LIGHTRED)
   yellow      P         YELLOW
   blue        ACFILMVW  MIDBLUE
   light blue  Y         LIGHTBLUE (16-colours: CYAN)
   green       NQST      GREEN
   purple      DE        PURPLE (16-colours: MAGENTA)
   red         RK        RED
   pink        H         LIGHTRED (16-colours: DARKRED)
   */
  
  color['A'] = color['a'] = MIDBLUE;
  color['B'] = color['b'] = NOCOLOR;
  color['C'] = color['c'] = MIDBLUE;
  color['D'] = color['d'] = PURPLE;
  color['E'] = color['e'] = PURPLE;
  color['F'] = color['f'] = MIDBLUE;
  color['G'] = color['g'] = ORANGE;
  color['H'] = color['h'] = LIGHTRED;
  color['I'] = color['i'] = MIDBLUE;
  color['K'] = color['k'] = RED;
  color['L'] = color['l'] = MIDBLUE;
  color['M'] = color['m'] = MIDBLUE;
  color['N'] = color['n'] = GREEN;
  color['P'] = color['p'] = YELLOW;
  color['Q'] = color['q'] = GREEN;
  color['R'] = color['r'] = RED;
  color['S'] = color['s'] = GREEN;
  color['T'] = color['t'] = GREEN;
  color['V'] = color['v'] = MIDBLUE;
  color['W'] = color['w'] = MIDBLUE;
  color['Y'] = color['y'] = LIGHTBLUE;
  color['Z'] = color['z'] = NOCOLOR;
}


/* Save the current residue colors as the 'custom' color scheme */
void saveCustomColors(BelvuContext *bc)
{
  bc->haveCustomColors = TRUE;
  
  customColor['A'] = customColor['a'] = color['a'];
  customColor['B'] = customColor['b'] = color['b'];
  customColor['C'] = customColor['c'] = color['c'];
  customColor['D'] = customColor['d'] = color['d'];
  customColor['E'] = customColor['e'] = color['e'];
  customColor['F'] = customColor['f'] = color['f'];
  customColor['G'] = customColor['g'] = color['g'];
  customColor['H'] = customColor['h'] = color['h'];
  customColor['I'] = customColor['i'] = color['i'];
  customColor['K'] = customColor['k'] = color['k'];
  customColor['L'] = customColor['l'] = color['l'];
  customColor['M'] = customColor['m'] = color['m'];
  customColor['N'] = customColor['n'] = color['n'];
  customColor['P'] = customColor['p'] = color['p'];
  customColor['Q'] = customColor['q'] = color['q'];
  customColor['R'] = customColor['r'] = color['r'];
  customColor['S'] = customColor['s'] = color['s'];
  customColor['T'] = customColor['t'] = color['t'];
  customColor['V'] = customColor['v'] = color['v'];
  customColor['W'] = customColor['w'] = color['w'];
  customColor['Y'] = customColor['y'] = color['y'];
  customColor['Z'] = customColor['z'] = color['z'];
  
}


/* Set the current colors the the last-saved 'custom' color scheme */
static void colorSchemeCustom(BelvuContext *bc)
{
  color['A'] = color['a'] = customColor['a'];
  color['B'] = color['b'] = customColor['b'];
  color['C'] = color['c'] = customColor['c'];
  color['D'] = color['d'] = customColor['d'];
  color['E'] = color['e'] = customColor['e'];
  color['F'] = color['f'] = customColor['f'];
  color['G'] = color['g'] = customColor['g'];
  color['H'] = color['h'] = customColor['h'];
  color['I'] = color['i'] = customColor['i'];
  color['K'] = color['k'] = customColor['k'];
  color['L'] = color['l'] = customColor['l'];
  color['M'] = color['m'] = customColor['m'];
  color['N'] = color['n'] = customColor['n'];
  color['P'] = color['p'] = customColor['p'];
  color['Q'] = color['q'] = customColor['q'];
  color['R'] = color['r'] = customColor['r'];
  color['S'] = color['s'] = customColor['s'];
  color['T'] = color['t'] = customColor['t'];
  color['V'] = color['v'] = customColor['v'];
  color['W'] = color['w'] = customColor['w'];
  color['Y'] = color['y'] = customColor['y'];
  color['Z'] = color['z'] = customColor['z'];
}


/* This is called when the color scheme has been set to 'by residue'. It updates
 * the colors according to the active color scheme. */
void setResidueSchemeColors(BelvuContext *bc)
{
  switch (bc->residueScheme)
    {
      case BELVU_SCHEME_ERIK:     colorSchemeErik(bc);	    break;
      case BELVU_SCHEME_GIBSON:   colorSchemeGibson(bc);    break;
      case BELVU_SCHEME_CGP:      colorSchemeCGP(bc);	    break;
      case BELVU_SCHEME_CGPH:     colorSchemeCGPH(bc);	    break;
      case BELVU_SCHEME_NONE:     colorSchemeEmpty(bc);	    break;
      case BELVU_SCHEME_CUSTOM:   colorSchemeCustom(bc);    break;
      
      default: 
	g_warning("Program error: unrecognised color scheme '%d'.\n", bc->residueScheme); 
	break;
    }
}


/* This should be called when the color scheme has changed. It updates the
 * colors according to the active scheme. */
void updateSchemeColors(BelvuContext *bc)
{
  /* Set the color scheme if coloring by conservation or if applying a 
   * threshold when coloring by residue */
  if (colorByConservation(bc) || colorByResId(bc))
    setConsSchemeColors(bc);
}


/* Return 1 if c1 has priority over c2, 0 otherwise */
static int colorPriority(BelvuContext *bc, int c1, int c2) 
{
  if (c2 == WHITE) 
    return 1;
  
  if (c2 == *getConsColor(bc, CONS_LEVEL_MAX, FALSE))
    return 0;
  
  if (c2 == *getConsColor(bc, CONS_LEVEL_LOW, FALSE)) 
    {
      if (c1 ==*getConsColor(bc, CONS_LEVEL_LOW, FALSE)) 
        return 0;
      else 
        return 1;
    }
  
  if (c2 == *getConsColor(bc, CONS_LEVEL_MID, FALSE)) 
    {
      if (c1 == *getConsColor(bc, CONS_LEVEL_MAX, FALSE)) 
        return 1;
      else
        return 0;
    }
  
  g_critical("Program error: invalid background colour '%s' when calculating color priority.\n", colorNames[c2]);
  
  return 0 ;
}


/* This is called when the color scheme type has been changed to 'by conservation'. It updates
 * the colors according to the active color scheme. */
void setConsSchemeColors(BelvuContext *bc)
{
  int i, j, k, l, colornr, simCount, n, nseqeff;
  double id, maxid;
  
  if (!bc->conservCount) 
    initConservMtx(bc);
  
  nseqeff = countResidueFreqs(bc);
  
  for (i = 0; i < bc->maxLen; ++i)
    {
      for (k = 1; k < 21; ++k)
	{
	  bc->colorMap[k][i] = WHITE;
	}
    }
  
  for (i = 0; i < bc->maxLen; ++i) 
    {
      maxid = -100.0;
    
      for (k = 1; k < 21; k++) 
	{
	  if (colorBySimilarity(bc)) 
	    {
              /* Convert counts to similarity counts */
              simCount = 0;
              for (j = 1; j < 21; j++) 
                {
                  if (j == k) 
                    {
                      if (1)
                        {
                          simCount += (bc->conservCount[j][i]-1) * bc->conservCount[k][i] * BLOSUM62[j-1][k-1];
                        }
                      else
                        {
                          /* Alternative, less good way */
                          simCount += 
                            (int)floor(bc->conservCount[j][i]/2.0)*
                            (int)ceil(bc->conservCount[k][i]/2.0)*
                            BLOSUM62[j-1][k-1];
                        }
                    }
                  else
                    {
                      simCount += bc->conservCount[j][i] * bc->conservCount[k][i] * BLOSUM62[j-1][k-1];
                    }
                }
	    
              if (bc->ignoreGapsOn) 
                n = bc->conservResidues[i];
              else 
                n = nseqeff;
              
              if (n < 2)
                {
                  id = 0.0;
                }
              else 
                {
                  if (1)
                    id = (double)simCount/(n*(n-1));
                  else
                    /* Alternative, less good way */
                    id = (double)simCount/(n/2.0 * n/2.0);
                }
              
              /* printf("%d, %c:  simCount= %d, id= %.2f\n", i, b2a[k], simCount, id); */
	    
              if (id > bc->lowSimCutoff) 
                {
                  if (id > bc->maxSimCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MAX, FALSE);
                  else if (id > bc->midSimCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MID, FALSE);
                  else
                    colornr = *getConsColor(bc, CONS_LEVEL_LOW, FALSE);
                  
                  if (colorPriority(bc, colornr, bc->colorMap[k][i]))
                    bc->colorMap[k][i] = colornr;
	      
                  /* Colour all similar residues too */
                  for (l = 1; l < 21; l++) 
                    {
                      if (BLOSUM62[k-1][l-1] > 0 && colorPriority(bc, colornr, bc->colorMap[l][i])) 
                        {
                          /*printf("%d: %c -> %c\n", i, b2a[k], b2a[l]);*/
                          bc->colorMap[l][i] = colornr;
                        }
                    }
                }
	    }
	  else 
	    {
              if (bc->ignoreGapsOn && bc->conservResidues[i] != 1)
                id = (double)bc->conservCount[k][i]/bc->conservResidues[i];
              else
                id = (double)bc->conservCount[k][i]/nseqeff;
              
              if (colorByResId(bc)) 
                {
                  if (id*100.0 > bc->colorByResIdCutoff)
                    bc->colorMap[k][i] = color[(unsigned char)(b2a[k])];
                  else
                    bc->colorMap[k][i] = WHITE;
                }
              else if (id > bc->lowIdCutoff) 
                {
                  if (id > bc->maxIdCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MAX, FALSE);
                  else if (id > bc->midIdCutoff) 
                    colornr = *getConsColor(bc, CONS_LEVEL_MID, FALSE);
                  else
                    colornr = *getConsColor(bc, CONS_LEVEL_LOW, FALSE);
                  
                  bc->colorMap[k][i] = colornr;
                  
                  if (bc->consScheme == BELVU_SCHEME_ID_BLOSUM) 
                    {
                      /* Colour all similar residues too */
                      for (l = 1; l < 21; l++) 
                        {
                          if (BLOSUM62[k-1][l-1] > 0 && colorPriority(bc, colornr, bc->colorMap[l][i])) 
                            {
                              /*printf("%d: %c -> %c\n", i, b2a[k], b2a[l]);*/
                              bc->colorMap[l][i] = colornr;
                            }
                        }
                    }
                }
	    }
	
	  if (id > maxid) 
	    {
	      maxid = id;
	    }
	}
      
      bc->conservation[i] = maxid;
    }
}


void initMarkupColors(void)
{
  markupColor['0'] = DARKBLUE;
  markupColor['1'] = BLUE;
  markupColor['2'] = MIDBLUE;
  markupColor['3'] = LIGHTBLUE;
  markupColor['4'] = VIOLET;
  markupColor['5'] = PALEBLUE;
  markupColor['6'] = PALECYAN;
  markupColor['7'] = CYAN;
  markupColor['8'] = CYAN;
  markupColor['9'] = CYAN;
  markupColor['A'] = markupColor['a'] = WHITE;
  markupColor['B'] = markupColor['b'] = RED;
  markupColor['C'] = markupColor['c'] = PALEYELLOW;
  markupColor['D'] = markupColor['d'] = WHITE;
  markupColor['E'] = markupColor['e'] = RED;
  markupColor['F'] = markupColor['f'] = WHITE;
  markupColor['G'] = markupColor['g'] = DARKGREEN;
  markupColor['H'] = markupColor['h'] = DARKGREEN;
  markupColor['I'] = markupColor['i'] = DARKGREEN;
  markupColor['J'] = markupColor['j'] = WHITE;
  markupColor['K'] = markupColor['k'] = WHITE;
  markupColor['L'] = markupColor['l'] = WHITE;
  markupColor['M'] = markupColor['m'] = WHITE;
  markupColor['N'] = markupColor['n'] = WHITE;
  markupColor['O'] = markupColor['o'] = WHITE;
  markupColor['P'] = markupColor['p'] = WHITE;
  markupColor['Q'] = markupColor['q'] = WHITE;
  markupColor['R'] = markupColor['r'] = WHITE;
  markupColor['S'] = markupColor['s'] = YELLOW;
  markupColor['T'] = markupColor['t'] = YELLOW;
  markupColor['V'] = markupColor['v'] = WHITE;
  markupColor['U'] = markupColor['u'] = WHITE;
  markupColor['W'] = markupColor['w'] = WHITE;
  markupColor['X'] = markupColor['x'] = WHITE;
  markupColor['Y'] = markupColor['y'] = WHITE;
  markupColor['Z'] = markupColor['z'] = NOCOLOR;
}


/* Save the current color-by-residue color scheme */
void saveResidueColorScheme(BelvuContext *bc, FILE *fil)
{
  if (!fil) 
    return;
  
  int i = 1;
  for (i = 1; i < 21; i++)
    {
      fprintf(fil, "%c %s\n", b2a[i], colorNames[color[(unsigned char)(b2a[i])]]);
    }
  
  fclose(fil);
}


void readResidueColorScheme(BelvuContext *bc, FILE *fil, int *colorarr)
{
  if (!fil) 
    return;

  char *cp=NULL, line[MAXLINE+1], setColor[MAXLINE+1];
  unsigned char c ;
  int i=0, colornr=0;
  
  while (!feof(fil)) 
    {
      if (!fgets (line, MAXLINE, fil)) 
	break;
    
      /* remove newline */
      if ((cp = strchr(line, '\n'))) 
        *cp = 0 ;
    
      /* Parse color of organism in tree 
         Format:  #=OS BLUE D. melanogaster*/
      if (!strncmp(line, "#=OS ", 5)) 
        {
          cp = line+5;
          sscanf(cp, "%s", setColor);
          
          for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
            {
              if (!strcasecmp(colorNames[i], setColor)) 
                colornr = i;
            }
          
          if (colornr == -1) 
            {
              printf("Unrecognized color: %s, using black instead.\n", setColor);
              colornr = BLACK;
            }
          
          while(*cp == ' ') cp++;
          while(*cp != ' ') cp++;
          while(*cp == ' ') cp++;
          
          /* Find organism and set its colour */
          ALN aln;
          initAln(&aln);
          aln.organism = cp;
          
          int ip = 0;
          if (!arrayFind(bc->organismArr, &aln, &ip, (void*)organism_order))
            g_critical("Cannot find organism \"%s\", specified in color code file. Hope that's ok", aln.organism);
          else
            g_array_index(bc->organismArr, ALN, ip).color = colornr;
        }
      
      /* Ignore comments */
      if (*line == '#') 
        continue;
      
      /* Parse character colours */
      if (sscanf(line, "%c%s", &c, setColor) == 2) 
        {
          c = toupper(c);
          for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
            {
              if (!strcasecmp(colorNames[i], setColor)) 
                colornr = i;
            }
          
          if (colornr == -1) 
            {
              printf("Unrecognized color: %s\n", setColor);
              colornr = 0;
            }
          
          colorarr[(unsigned char)(c)] = colornr;
          
          if (c > 64 && c <= 96)
            colorarr[(unsigned char)(c+32)] = colorarr[(unsigned char)(c)];
          else if (c > 96 && c <= 128)
            colorarr[(unsigned char)(c-32)] = colorarr[(unsigned char)(c)];
        }
    }
  
  fclose(fil);
  
  /* Store the custom colors */
  saveCustomColors(bc);
}


/* Utility function to get the color display name for a color number */
const char* getColorNumName(const int colorNum)
{
  g_assert(colorNum < NUM_TRUECOLORS);
  return colorNames[colorNum];
}

/* Utility function to get the file format name for a format id */
const char* getFileFormatString(const int formatId)
{
  g_assert(formatId < BELVU_NUM_FILE_FORMATS);
  return fileFormatNames[formatId];
}


int* getConsPrintColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground)
{
  int *result = &bc->lowbgPrintColor;

  switch (consLevel)
    {
      case CONS_LEVEL_MAX:
        result = foreground ? &bc->maxfgPrintColor : &bc->maxbgPrintColor;
        break;
        
      case CONS_LEVEL_MID:
        result = foreground ? &bc->midfgPrintColor : &bc->midbgPrintColor;
        break;
        
      case CONS_LEVEL_LOW:
        result = foreground ? &bc->lowfgPrintColor : &bc->lowbgPrintColor;
        break;

      default:
        g_critical("Program error: invalid conservation level '%d' when getting conservation colour.\n", consLevel);
        break;
    };

  return result;
}


/* Utility to get foreground/background the colour number for the given
 * conservation level. Returns a pointer to the value in the context. */
int* getConsColor(BelvuContext *bc, const BelvuConsLevel consLevel, const gboolean foreground)
{
  if (bc->printColorsOn)
    return getConsPrintColor(bc, consLevel, foreground);

  int *result = &bc->lowbgColor;

  switch (consLevel)
    {
      case CONS_LEVEL_MAX:
        result = foreground ? &bc->maxfgColor : &bc->maxbgColor;
        break;
        
      case CONS_LEVEL_MID:
        result = foreground ? &bc->midfgColor : &bc->midbgColor;
        break;
        
      case CONS_LEVEL_LOW:
        result = foreground ? &bc->lowfgColor : &bc->lowbgColor;
        break;

      default:
        g_critical("Program error: invalid conservation level '%d' when getting conservation colour.\n", consLevel);
        break;
    };

  return result;
}


/* Set whether the currently-selected sequence should be excluded
 * from the conservation colours calculations */
void setExcludeFromConsCalc(BelvuContext *bc, const gboolean exclude)
{
  if (!bc->selectedAln) 
    {
      g_critical("Please select a sequence first.\n");
      return;
  }
  
  if ((exclude && bc->selectedAln->nocolor) ||
     (!exclude && !bc->selectedAln->nocolor))
    {
      /* Nothing to do */
      return;
    }
  
  /* Store orignal state in the nocolor field:
   1 = normal sequence
   2 = markup line
   
   This is needed to restore markup lines (nocolor lines are always markups)
   */
  
  if (exclude) 
    {
      if (bc->selectedAln->markup)
        bc->selectedAln->nocolor = 2;
      else
        bc->selectedAln->nocolor = 1;
      
      bc->selectedAln->markup = 1;
    }
  else 
    {
      if (bc->selectedAln->nocolor == 1) 
        bc->selectedAln->markup = 0;
      
      bc->selectedAln->nocolor = 0;
    }
}

/***********************************************************
 *		          Arrays			   *
 ***********************************************************/

/* Finds Entry s from Array  a
 * sorted in ascending order of order()
 * If found, returns TRUE and sets *ip
 * if not, returns FALSE and sets *ip one step left
 */
gboolean arrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer))
{
  int ord;
  int i = 0 , j = a->len, k;

  if (!j || (ord = orderFunc(s, &g_array_index(a, ALN, 0))) < 0)
    { 
      if (ip)
	*ip = -1; 
      return FALSE;
    }   /* not found */

  if (ord == 0)
    { 
      if (ip)
	*ip = 0;
      return TRUE;
    }

  if ((ord = orderFunc(s, &g_array_index(a, ALN, --j))) > 0)
    {
      if (ip)
	*ip = j; 
      return FALSE;
    }
  
  if (ord == 0)
    { 
      if (ip)
	*ip = j;
      return TRUE;
    }

  while (TRUE)
    { 
      k = i + ((j-i) >> 1) ; /* midpoint */

      if ((ord = orderFunc(s, &g_array_index(a, ALN, k))) == 0)
	{ 
          if (ip)
	    *ip = k; 
	  return TRUE;
	}
      
      if (ord > 0) 
	(i = k);
      else
	(j = k) ;
      
      if (i == (j-1) )
        break;
    }
  
  if (ip)
    *ip = i ;
  
  return FALSE;
}


/* Linear search for an exact matching sequence name and coordinates,
 typically to find back highlighted row after sorting 
 */
gboolean alignFind(GArray *alignArr, ALN *obj, int *idx)
{
  for (*idx = 0; *idx < alignArr->len; (*idx)++) 
    {
      ALN *currentAln = &g_array_index(alignArr, ALN, *idx);
      
      if (alphaorder(currentAln, obj) == 0) 
        return TRUE;
    }

  return FALSE;
}


/* Set the ->nr field in Align array in ascending order */
void arrayOrder(GArray *alignArr)
{
  int i = 0;
  
  for (i = 0; i < alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(alignArr, ALN, i);
      alnp->nr = i + 1;
    }
}


/* Set the ->nr field in Align array in ascending order with increments of 10
   Necessary if one wants to insert items.
 */
static void arrayOrder10(GArray *alignArr)
{
  int i;

  for (i = 0; i < alignArr->len; ++i) 
    g_array_index(alignArr, ALN, i).nr = (i+1)*10;
}



/* Separate markuplines to another array before resorting
 */
void separateMarkupLines(BelvuContext *bc)
{
  bc->markupAlignArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN), 100);
  
  ALN aln;
  initAln(&aln);
  
  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  arrayOrder(bc->alignArr);
  
  int i = 0;
  for (i = 0; i < bc->alignArr->len; ) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (alnp->markup) 
        {
          /* printf ("Moving line %d, %s/%d-%d, nseq=%d\n", i, alnp->name, alnp->start, alnp->end, nseq);*/
          g_array_append_val(bc->markupAlignArr, *alnp);
          g_array_sort(bc->markupAlignArr, alphaorder);
          g_array_remove_index(bc->alignArr, i);
        }
      else
        {
          ++i;
        }
    }
  
  arrayOrder(bc->alignArr);
  
  if (bc->selectedAln) 
    {
    int idx = 0;
    
    if (!alignFind(bc->alignArr, &aln, &idx))
      bc->selectedAln = 0;
    else
      bc->selectedAln = &g_array_index(bc->alignArr, ALN, idx);
    }
}


/* Reinsert markuplines after mother sequence or if orphan at bottom */
void reInsertMarkupLines(BelvuContext *bc)
{
  int i, j;
  char tmpname[MAXNAMESIZE+1], *cp;
  
  for (i = bc->markupAlignArr->len - 1; i >=0 ; --i)
    {
    ALN *alnp = &g_array_index(bc->markupAlignArr, ALN, i);
    strcpy(tmpname, alnp->name);
    
    if ((cp = strchr(tmpname, ' ')))
      *cp = 0;
    
    arrayOrder10(bc->alignArr);
    
    for (j = 0; j < bc->alignArr->len; ++j)
      {
      if (!strcmp(tmpname, g_array_index(bc->alignArr, ALN, j).name))
	break;
      }
    
    alnp->nr = (j+1)*10+5;
    
    g_array_append_val(bc->alignArr, *alnp);
    g_array_sort(bc->alignArr, nrorder); /* to do: can we move this out of the loop ? */
    }
  
  arrayOrder(bc->alignArr);
}


int strcmp_(gconstpointer xIn, gconstpointer yIn)
{
  const char *x = (const char*)xIn;
  const char *y = (const char*)yIn;
  
  int retval = strcmp(x, y);
  return retval;
}


GArray *copyAlignArray(GArray *inputArr)
{
  GArray *result = g_array_sized_new(FALSE, FALSE, sizeof(ALN), inputArr->len);
  memcpy(result->data, inputArr->data, inputArr->len * sizeof(ALN));
  
  int i = 0;
  for ( ; i < inputArr->len; ++i)
    {
      ALN *inputAln = &g_array_index(inputArr, ALN, i);
      
      if (alnGetSeq(inputAln))
        g_array_index(result, ALN, i).sequenceStr = g_string_new(alnGetSeq(inputAln));
      else
        g_array_index(result, ALN, i).sequenceStr = NULL;
    }
  
  return result;
}


void columnCopy(GArray *alignArrDest, int destIdx, GArray *alignArrSrc, int srcIdx)
{
  int i;
  
  for (i = 0; i < alignArrSrc->len; ++i)
    {
      ALN *srcAln = &g_array_index(alignArrSrc, ALN, i);
      ALN *destAln = &g_array_index(alignArrDest, ALN, i);

      char *srcSeq = alnGetSeq(srcAln);
      char *destSeq = alnGetSeq(destAln);
      
      if (srcSeq && destSeq && destIdx < alnGetSeqLen(destAln) && srcIdx < alnGetSeqLen(srcAln))
        destSeq[destIdx] = srcSeq[srcIdx];
    }
}



/***********************************************************
 *		          Context			   *
 ***********************************************************/

/* Create the context, which contains all program-wide variables */
BelvuContext* createBelvuContext()
{
  BelvuContext *bc = g_malloc(sizeof *bc);
  
  bc->belvuWindow = NULL;
  bc->spawnedWindows = NULL;
  bc->belvuTree = NULL;
  bc->belvuAlignment = NULL;
  bc->consPlot = NULL;
  bc->orgsWindow = NULL;
  
  bc->defaultCursor = NULL; /* get from gdkwindow once it is shown */
  bc->removeSeqsCursor = gdk_cursor_new(GDK_PIRATE);
  bc->busyCursor = gdk_cursor_new(GDK_WATCH);

  bc->defaultColors = NULL;
  
  bc->alignArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN), 100);  /* was called 'Align' */
  bc->organismArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN), 100);
  bc->markupAlignArr = NULL;
  bc->bootstrapGroups = NULL;
  
  bc->selectedAln = NULL;
  bc->highlightedAlns = NULL;
  
  bc->treeHead = NULL;
  bc->treeBestBalancedNode = NULL;
  
  bc->treeReadDistancesPipe = NULL;
  
  bc->IN_FORMAT = MUL;
  bc->maxScoreLen = 0;
  bc->alignYStart = 0;
  bc->treebootstraps = 0; 
  bc->maxLen = 0;
  bc->maxTreeWidth = 0;
  bc->maxNameLen = 0;   
  bc->maxStartLen = 0; 
  bc->maxEndLen = 0; 
  bc->maxScoreLen = 0; 
  bc->selectedCol = 0;
  bc->highlightedCol = 0;
  
  bc->maxfgColor = DEFAULT_MAX_FG_COLOR;
  bc->midfgColor = DEFAULT_MID_FG_COLOR,
  bc->lowfgColor = DEFAULT_LOW_FG_COLOR;
  bc->maxbgColor = DEFAULT_MAX_BG_COLOR;
  bc->midbgColor = DEFAULT_MID_BG_COLOR;
  bc->lowbgColor = DEFAULT_LOW_BG_COLOR;
  bc->maxfgPrintColor = DEFAULT_MAX_FG_PRINT_COLOR;
  bc->midfgPrintColor = DEFAULT_MID_FG_PRINT_COLOR,
  bc->lowfgPrintColor = DEFAULT_LOW_FG_PRINT_COLOR;
  bc->maxbgPrintColor = DEFAULT_MAX_BG_PRINT_COLOR;
  bc->midbgPrintColor = DEFAULT_MID_BG_PRINT_COLOR;
  bc->lowbgPrintColor = DEFAULT_LOW_BG_PRINT_COLOR;
  
  bc->schemeType = BELVU_SCHEME_TYPE_RESIDUE;
  bc->residueScheme = BELVU_SCHEME_ERIK;
  bc->consScheme = BELVU_SCHEME_BLOSUM;
  
  bc->treeMethod = NJ;
  bc->treeDistCorr = SCOREDIST;
  bc->treePickMode = NODESWAP;
  
  bc->sortType = BELVU_UNSORTED;
  
  bc->annotationList = NULL;
  
  bc->treeBestBalance = 0.0;
  bc->treeBestBalance_subtrees = 0.0;
  bc->tree_y = 0.3;
  bc->lowIdCutoff = DEFAULT_LOW_ID_CUTOFF;
  bc->midIdCutoff = DEFAULT_MID_ID_CUTOFF;
  bc->maxIdCutoff = DEFAULT_MAX_ID_CUTOFF;
  bc->lowSimCutoff = DEFAULT_LOW_SIM_CUTOFF;
  bc->midSimCutoff = DEFAULT_MID_SIM_CUTOFF;
  bc->maxSimCutoff = DEFAULT_MAX_SIM_CUTOFF;
  bc->colorByResIdCutoff = 20.0;
  bc->mksubfamilies_cutoff = 0.0;
  bc->treeScale = 0.3;
  bc->treeLineWidth = 0.3;
  
  bc->gapChar = '.';
  bc->saveSeparator = '/';
  strcpy(bc->treeDistString, SCOREDISTstr);
  strcpy(bc->treeMethodString, NJstr);
  bc->Title[0] = '\0';
  bc->saveFormat = BELVU_FILE_MUL;
  bc->fileName = NULL;
  bc->dirName = NULL;
  bc->organismLabel[0] = 'O';
  bc->organismLabel[0] = 'S';   
  bc->organismLabel[0] = '\0'; 
  
  bc->conservCount = NULL;
  bc->colorMap = NULL;
  bc->conservResidues = NULL;
  bc->conservation = NULL;
  
  bc->treeCoordsOn = TRUE;
  bc->treeReadDistancesOn = FALSE;
  bc->treePrintDistances = FALSE;
  bc->penalize_gaps = FALSE;
  bc->stripCoordTokensOn = TRUE;
  bc->saveCoordsOn = TRUE;
  bc->displayScores = FALSE;
  bc->outputBootstrapTrees = FALSE;
  bc->treebootstrapsDisplay = FALSE;
  bc->treeColorsOn = TRUE;
  bc->treeShowOrganism = TRUE;
  bc->treeShowBranchlen = FALSE;
  bc->matchFooter = FALSE;
  bc->saved = TRUE;
  bc->ignoreGapsOn = FALSE;
  bc->colorByResIdOn = FALSE;
  bc->rmEmptyColumnsOn = TRUE;
  bc->lowercaseOn = FALSE;
  bc->removingSeqs = FALSE;
  bc->displayColors = TRUE;
  bc->haveCustomColors = FALSE;
  bc->printColorsOn = FALSE;
  bc->highlightOrthologs = FALSE;
  bc->useWWWFetch = FALSE;
  
  /* Null out all the entries in the dialogs list */
  int dialogId = 0;
  for ( ; dialogId < BELDIALOG_NUM_DIALOGS; ++dialogId)
    {
      bc->dialogList[dialogId] = NULL;
    }
  
  return bc;
}


/* Destroy the context */
void destroyBelvuContext(BelvuContext **bc)
{
  if (bc && *bc)
    {
    if ((*bc)->alignArr)
      g_array_unref((*bc)->alignArr);
    
    if ((*bc)->organismArr)
      g_array_unref((*bc)->organismArr);
    
    if ((*bc)->markupAlignArr)
      g_array_unref((*bc)->markupAlignArr);
    
    if ((*bc)->bootstrapGroups)
      g_array_unref((*bc)->bootstrapGroups);
    
    g_free(*bc);
    *bc = NULL;
    }
}


/***********************************************************
 *                           Utilities                     *
 ***********************************************************/

/* Utility to return the sequence data in the given alignment; returns null if
 * not set. */
char* alnGetSeq(ALN *aln)
{
  return (aln && aln->sequenceStr ? aln->sequenceStr->str : NULL);
}


/* Utility to return the length sequence data in the given alignment; returns
 * 0 if not set. */
int alnGetSeqLen(ALN *aln)
{
  return (aln && aln->sequenceStr ? aln->sequenceStr->len : 0);
}


/* This function just returns the value in the b2a array at the given index */
char b2aIndex(const int idx)
{
  return b2a[idx];
}


void drawText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const char *text, int *textWidth, int *textHeight)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, text);
  
  if (drawable)
    gdk_draw_layout(drawable, gc, x, y, layout);
  
  /* Return the width and height of the layout, if requested */
  pango_layout_get_size(layout, textWidth, textHeight);
  
  if (textWidth)
    *textWidth /= PANGO_SCALE;
  
  if (textHeight)
    *textHeight /= PANGO_SCALE;
  
  g_object_unref(layout);
}


/* Utility to draw the given integer as text. The text is right-aligned, so 
 * the input x coord must be the RIGHTMOST EDGE of where you want the text. */
void drawIntAsText(GtkWidget *widget, 
                   GdkDrawable *drawable, 
                   GdkGC *gc, 
                   const int x, 
                   const int y, 
                   const int value)
{
  char *tmpStr = blxprintf("%d", value);
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, tmpStr);
  g_free(tmpStr);

  int textWidth;
  pango_layout_get_pixel_size(layout, &textWidth, NULL);

  gdk_draw_layout(drawable, gc, x - textWidth, y, layout);
  g_object_unref(layout);
}


/* Utility to draw the given double as text. The text is right-aligned, so 
 * the input x coord must be the RIGHTMOST EDGE of where you want the text. */
void drawDoubleAsText(GtkWidget *widget, 
                      GdkDrawable *drawable, 
                      GdkGC *gc, 
                      const int x, 
                      const int y, 
                      const double value)
{
  char *tmpStr = blxprintf("%.1f", value);
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, tmpStr);
  g_free(tmpStr);
  
  int textWidth;
  pango_layout_get_pixel_size(layout, &textWidth, NULL);
  
  gdk_draw_layout(drawable, gc, x - textWidth, y, layout);
  g_object_unref(layout);
}


/* Calculate percent identity of two strings */
double identity(char *s1, char *s2, const gboolean penalize_gaps)
{
    int n, id;

    for (n = id = 0; *s1 && *s2; s1++, s2++) 
      {
	if (isGap(*s1) && isGap(*s2)) 
          continue;
        
	if (isGap(*s1) || isGap(*s2))
          {
	    if (!penalize_gaps) 
              continue;
          }
        
	n++;
	if (toupper(*s1) == toupper(*s2)) 
          id++;
      }

    if (n)
      return (double)id/n*100;
    else
      return 0.0;
}


/* Calculate score of two sequences */
static double score(char *s1, char *s2, const gboolean penalize_gaps)
{
    double sc=0.0;

    for (;*s1 && *s2; s1++, s2++) 
      {
	if (isGap(*s1) && isGap(*s2)) 
          {
	    continue;
          }
	else if (isGap(*s1) || isGap(*s2)) 
          {
	    if (penalize_gaps) sc -= 0.6;
          }
	else
          {
            int val1 = a2b[(unsigned char)(*s1)];
            int val2 = a2b[(unsigned char)(*s2)];
            
            if (val1 > 0 && val2 > 0)
              sc += (double) BLOSUM62[val1 - 1][val2 - 1];
          }
      }
    
    return sc;
}


gboolean isGap(char c) 
{
  if (c == '.' || c == '-' ||
      c == '[' || c == ']' /* Collapse-control chars */ ) 
    return TRUE;
  else
    return FALSE;
}


static gboolean isAlign(char c)
{
  if (isalpha(c) || isGap(c) || c == '*')
    return TRUE;
  else 
    return FALSE;
}


int GCGchecksum(BelvuContext *bc, char *seq)
{
  int  
  check = 0, 
  count = 0, 
  i=0;
  
  for (i = 0; i < bc->maxLen; ++i) 
    {
    ++count;
    check += count * toupper((int) seq[i]);
    
    if (count == 57) 
      count = 0;
    }
  
  return (check % 10000);
}


int GCGgrandchecksum(BelvuContext *bc)
{
  int 
  i=0,
  grand_checksum=0;
  
  for(i=0; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      grand_checksum += GCGchecksum(bc, alnGetSeq(alnp));
    }
  
  return (grand_checksum % 10000);
}


/* Read labels of highlighted sequence and spread them */
void readLabels(BelvuContext *bc, FILE *fil)
{
  char *labelseq = g_malloc(bc->maxLen + 1); /* The raw sequence of labels */
  char *label = g_malloc(bc->maxLen + 1);    /* The mapped sequence of labels, 1-maxlen */
  char line[MAXLENGTH+1];
  
  /* read file */
  char *cq = labelseq;
  while (!feof (fil)) 
    {
      if (!fgets (line, MAXLENGTH, fil)) 
        break;
      
      char *cp = line;
      while (*cp) 
        {
          if (isalpha(*cp)) 
            *cq++ = *cp;
          cp++;
        }
    }

  fclose(fil);
  
  /* Warn if seq too long, return if too short */
  int seqlen = bc->selectedAln->end - bc->selectedAln->start + 1;
  if (strlen(labelseq) > seqlen)
    {
      g_critical("The sequence of labels is longer (%d) than the sequence (%d).\nHope that's ok", 
		 (int)strlen(labelseq), seqlen);
    }
  else if (strlen(labelseq) < seqlen) 
    {
      g_critical("The sequence of labels is shorter (%d) than the sequence (%d).\nAborting", 
		 (int)strlen(labelseq), seqlen);
      return;
    }
  
  /* map labels to alignment */
  int col = 0;
  int seqpos = 0;
  
  for (col = 0, seqpos = 0; col < bc->maxLen; col++) 
    {
      label[col] = labelseq[seqpos];
      const char *selectedSeq = alnGetSeq(bc->selectedAln);
      
      if (selectedSeq && isalpha(selectedSeq[col]) && labelseq[seqpos+1]) 
        seqpos++;
    }
  
  int row = 0;
  for (row = 0; row < bc->alignArr->len; ++row) 
    {
      ALN *alnrow = &g_array_index(bc->alignArr, ALN, row);
      
      for (col = 0; col < bc->maxLen; col++)
        {
          char *rowSeq = alnGetSeq(alnrow);
          
          if (rowSeq && isalpha(rowSeq[col]))
            rowSeq[col] = label[col];
        }
    }
  
  g_free(labelseq);
}

/***********************************************************
 *		          			   *
 ***********************************************************/

static void makeSegList(BelvuContext *bc, SEG **SegList, char *line)
{
  int n, i;
  char *p, *linecopy;
  SEG *seg, *prevseg=0;

  /* Count coordinates - has to be multiple of 4 */
  linecopy = g_malloc(strlen(line)+1);
  strcpy(linecopy, line);
  n = 0;
  
  if (atoi(strtok(linecopy, " "))) 
    n++;;
  
  while ( (p = strtok(0, " ")) && atoi(p) ) 
    n++;
  
  g_free(linecopy);
    
  if (!n || n % 4) 
    g_error("Segments not multiple of 4 ints (%s)", line);

  for (i = 0; i < n/4; i++) 
    {
      seg = (SEG *)g_malloc(sizeof(SEG));
      if (prevseg) 
        prevseg->next = seg;
      else
        *SegList = seg;
      
      prevseg = seg;

      seg->qstart = atoi(i ? strtok(0, " ") : strtok(line, " "));
      seg->qend = atoi(strtok(0, " "));
      seg->start = atoi(strtok(0, " "));
      seg->end = atoi(strtok(0, " "));
      seg->next = NULL;
    }

  for (seg = *SegList; seg; seg = seg->next) 
    {
      DEBUG_OUT("%d %d %d %d\n", seg->qstart, seg->qend, seg->start, seg->end);

      if (seg == *SegList && seg->qstart != 1)
        g_error("Bad qstart: Must start on 1");

      if (seg->qstart < 1 || seg->qstart > bc->maxLen)
        g_error("Bad qstart: %d.  Range: 1-%d", seg->qstart, bc->maxLen);
      if (seg->qend < 1 || seg->qend > bc->maxLen)
        g_error("Bad qend: %d.  Range: 1-%d", seg->qend, bc->maxLen);
      if (seg->start < 1 || seg->start > bc->maxLen)
        g_error("Bad start: %d.  Range: 1-%d", seg->start, bc->maxLen);
      if (seg->end < 1 || seg->end > bc->maxLen)
        g_error("Bad end: %d.  Range: 1-%d", seg->end, bc->maxLen);

      if (seg->qstart > seg->qend)
        g_error("qstart > qend  (%d > %d)", seg->qstart, seg->qend);
      if (seg->start > seg->end)
        g_error("start > end  (%d > %d)", seg->start, seg->end);
    }
}


static int countInserts(SEG *seg)
{
  int   
    Align_gap, Query_gap, gap, insert_counter=0;

  if (!seg)
    return insert_counter ;

  while (seg->next) 
    {
      Align_gap = seg->next->start - seg->end - 1;
      if (Align_gap < 0) 
        g_error("Negative Align_gap: %d (%d-%d)", Align_gap, seg->start, seg->next->end );
		
      Query_gap = seg->next->qstart - seg->qend - 1;
      if (Query_gap < 0) 
        g_error("Negative Query_gap: %d (%d-%d)", Query_gap, seg->qstart, seg->next->qend );

      gap = Query_gap - Align_gap;
      if (gap > 0) 
        {
          insert_counter += gap;
	}
      
      seg = seg->next;
    }
  
  return insert_counter;
}


/* Inserts n gap columns after position p, in sequence coordinate, 1...maxLen 
 */
static void insertColumns(BelvuContext *bc, int p, int n)
{
  int
    i, j;

  ALN 
    *alni;
  char
    *dest, *src, *seq;

  printf("Inserting %d columns after column %d\n", n, p);

  bc->maxLen += n;

  for (i = 0; i < bc->alignArr->len; ++i) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
	
      seq = g_malloc(bc->maxLen + 1);

      dest = seq;
      src = alnGetSeq(alni);

      for (j = 0; j < p;  j++) 
        {
          *dest++ = *src++;
	}
      for (; j < p+n;  j++) 
        {
          *dest++ = '.';
	}
      for (; j < bc->maxLen;  j++) 
        {
          *dest++ = *src++;
	}

      g_string_free(alni->sequenceStr, TRUE);
      alni->sequenceStr = g_string_new(seq);
    }

  bc->saved = FALSE;
}


/* readMatch displays a matching sequences to the alignment so that
   the match can be displayed in relation to the whole family.

   Considerations:

   - Must work both from ACEDB and command line (and from Web server).

   Display:

   - The score should be visible (turn on score column.)

   - Should matches be displayed on separate lines or in the middle of the alignment?
     On top would be nice, but would require a lot of extra programming (status bar etc).
     Instead add to alignment (SEED usually fits on one screen) 
     Draw names with red background.  
     Add on top.


   Steps needed on command line:

   hmmb -W tmp HMM align
   selex2alignmap tmp > alignmap

   hmm?s -F HMM query | ?smapback map=alignmap > query.mapmatch  ( -> acedb)
   belvu -m query.mapmatch align

   The two last lines should be scripted in hmm?sBelvu


*/
void readMatch(BelvuContext *bc, FILE *fil)
{
  /* Format:
       
  Using segments is more difficult than residues, but is
  necessary to display in ACEDB Pepmap.

  seg_match/3-30                   (matching segment name)
  WLPLHTLinsertAACGEFYLVDSLLKH     (only matching part, no pads!)
  1 7 3 9  14 28 10 24

  seq_match2...
  */
  int	orig_maxLen, tmp,
    Align_gap, Query_gap, gap, inserts, seglen,
    done_one=0, len,  was_saved=bc->saved, insert_counter;
  char *cp, *seq, *rawseq = NULL, *seqp;
  SEG	*SegList = NULL, *seg;
  char line[MAXLENGTH+1];
  
  while (!done_one)
    {
      if (!fgets(line, MAXLENGTH, fil))
	break;

      /*printf("%s\n", line); continue;*/

      if (*line != '\n' && *line != '#')
	{
	  if (done_one)
	    {
	      g_critical("Sorry, can only do one match");
	      continue;
	    }

          ALN aln;
          initAln(&aln);

	  /* Name */
	  if (!strtok(line, "/"))
	    g_error("Bad format: %s", line);

	  strncpy(aln.name, line, MAXNAMESIZE);
	  aln.name[MAXNAMESIZE] = 0;
          
	  if ((cp = strchr(aln.name, ' ')))
	    *cp = 0;

	  if (bc->maxNameLen  < strlen(aln.name))
	    bc->maxNameLen = strlen(aln.name);

	  /* Start & End */
	  if (!(cp = strtok(0, "-"))) 
            g_error("Bad start: %s", cp);

	  aln.start = atoi(cp);
          
	  if ((len=strlen(cp)) > bc->maxStartLen)
	    bc->maxStartLen = len;
          
          char *tmpStr = blxprintf("%d", aln.start);

	  if (bc->maxStartLen < (tmp = strlen(tmpStr)))
	    bc->maxStartLen = tmp;
                                   
          g_free(tmpStr);
          tmpStr = NULL;

	  if (!(cp = strtok(0, " ")))
	    g_error("Bad end: %s", cp);
          
	  aln.end = atoi(cp);
          
	  if ((len=strlen(cp)) > bc->maxEndLen)
	    bc->maxEndLen = len;
          
          tmpStr = blxprintf("%d", aln.end);
          
	  if (bc->maxEndLen < (tmp = strlen(tmpStr)))
	    bc->maxEndLen = tmp;
    
          g_free(tmpStr);
          tmpStr = NULL;
          
          int ip = 0;
	  if (alignFind(bc->alignArr, &aln, &ip))
	    {
	      g_critical("Sorry, already have sequence %s/%d-%d\n", 
                         aln.name, aln.start, aln.end);
	      return;
	    }

	  /* Score */
	  if (!(cp = strtok(0, "\n")))
	    g_error("Bad score: %s", cp);

	  aln.score = atof(cp);
          
          tmpStr = blxprintf("%.1f", aln.score);
          
	  if ((len=strlen(tmpStr)) > bc->maxScoreLen)
	    bc->maxScoreLen = len;
          
          g_free(tmpStr);
          tmpStr = NULL;

	  bc->displayScores = TRUE;

	  /* Sequence */
	  if (!fgets(line, MAXLENGTH, fil))
	    break;

	  if ((cp = strchr(line, '\n')))
	    *cp = 0 ;
          
	  rawseq = g_malloc(strlen(line)+1);
	  strcpy(rawseq, line);

	  /* Matching segments */
          
	  if (!fgets(line, MAXLENGTH, fil))
	    break;

	  if ((cp = strchr(line, '\n')))
            *cp = 0;

	  makeSegList(bc, &SegList, line);

	  inserts = countInserts(SegList);
	  seq = g_malloc(bc->maxLen + inserts + 1);
	  memset(seq, '.', bc->maxLen + inserts);
	    
	  orig_maxLen = bc->maxLen;

	  /* Add first segment */
	  seg = SegList;

	  seqp = seq + seg->start-1;
	  seglen = seg->end - seg->start+1;
	  strncpy(seqp, rawseq, seglen);
	  seqp += seglen; 

	  /* Add subsequent segments */
	  insert_counter = 0;
	  while (seg->next)
	    {
	      Align_gap = seg->next->start - seg->end - 1;
	      Query_gap = seg->next->qstart - seg->qend - 1;


	      if (Query_gap > 0)
		{ 
		  strncpy(seqp, rawseq + seg->qend, Query_gap);
		  seqp += Query_gap; 
		}

	      gap = Query_gap - Align_gap;
	      if (gap > 0)
		{
		  insertColumns(bc, seg->end+insert_counter, gap);
		  insert_counter += gap;
		}
	      else if (gap < 0)
		{
		  seqp += -gap;
		  DEBUG_OUT("inserting %d gaps at %d\n", -gap, (int)(seqp-seq));
		}

	      seg = seg->next;

	      /* Add sequence segment */
	      seglen = seg->end - seg->start+1;
	      strncpy(seqp, rawseq + seg->qstart-1, seglen);
	      seqp += seglen;
	    }

	  /* Add final pads */
	  memset(seqp, '.', orig_maxLen - seg->end);

	  aln.sequenceStr = g_string_new(seq);
          g_free(seq);
          
	  aln.color = RED;
	  aln.nr = 0;

          g_array_append_val(bc->alignArr, aln);
          g_array_sort(bc->alignArr, nrorder); /* to do: can we do the sort outside the loop? */
          
	  done_one = 1;

	} /* End of this matchSeq */
    }

  g_free(rawseq);

  arrayOrder(bc->alignArr);
  bc->saved = was_saved;

  return ;
}


void checkAlignment(BelvuContext *bc)
{
  int i, g, cres, nres, tmp;

  bc->maxNameLen = bc->maxStartLen = bc->maxEndLen = 0;

  for (i = 0; i < bc->alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);

      if (!alnp->markup) 
        {
          char *alnSeq = alnGetSeq(alnp);
          
          /* Count residues */
          for (cres = g = 0; g < bc->maxLen; g++) 
            {
              if (alnSeq && isAlign(alnSeq[g]) && !isGap(alnSeq[g])) 
                cres++;
	    }
	    
          if (!alnp->start) 
            {
              /* No coords provided - reconstruct them */
              alnp->start = 1;
              alnp->end = cres;
	    }
          else 
            {
              /* Check if provided coordinates are correct */
              nres = abs(alnp->end - alnp->start) + 1;
              
              if ( nres != cres)
                {
                  fprintf(stderr, "Found wrong number of residues in %s/%d-%d: %d instead of %d\n",
                          alnp->name, alnp->start, alnp->end,
                          cres,
                          nres);
                }
	    }
	}

      /* Find max string length of name, start, and end */
      if (bc->maxNameLen  < strlen(alnp->name))
        bc->maxNameLen = strlen(alnp->name);

      char *tmpStr = blxprintf("%d", alnp->start);
      
      if (bc->maxStartLen < (tmp = strlen(tmpStr))) 
        bc->maxStartLen = tmp;

      g_free(tmpStr);
      tmpStr = NULL;

      tmpStr= blxprintf("%d", alnp->end);

      if (bc->maxEndLen   < (tmp = strlen(tmpStr))) 
        bc->maxEndLen = tmp;
    }
}    


static void initConservMtx(BelvuContext *bc)
{
  int i;

  bc->conservCount = (int **)g_malloc(21*sizeof(int *));
  bc->colorMap = (int **)g_malloc(21*sizeof(int *));
    
  for (i = 0; i < 21; ++i) 
    {
      bc->conservCount[i] = (int *)g_malloc(bc->maxLen*sizeof(int));
      bc->colorMap[i] = (int *)g_malloc(bc->maxLen*sizeof(int));
    }

  bc->conservResidues = (int*)g_malloc(bc->maxLen*sizeof(int));
  bc->conservation = (double*)g_malloc(bc->maxLen*sizeof(double));
}	


/* Calculate conservation in each column */
static int countResidueFreqs(BelvuContext *bc)
{
  int   i, j, nseqeff=0;

  if (!bc->conservCount) 
    initConservMtx(bc);

    /* Must reset since this routine may be called many times */
  for (i = 0; i < bc->maxLen; ++i)
    {
      for (j = 0; j < 21; j++)
        {
          bc->conservCount[j][i] = 0;
          bc->colorMap[j][i] = 0;
        }
      
      bc->conservResidues[i] = 0;
    }
    
  for (j = 0; j < bc->alignArr->len; ++j)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, j);
      
      if (alnp->markup)
        continue;

      nseqeff++;

      for (i = 0; i < bc->maxLen; ++i)
        {
          char *alnSeq = alnGetSeq(alnp);
          
          bc->conservCount[a2b[(unsigned char)(alnSeq[i])]][i]++;

          if (isalpha(alnSeq[i]) || alnSeq[i] == '*')
            bc->conservResidues[i]++;
        }
    }

  return nseqeff;
}


/* Calculate overhang between two aligned sequences 
   Return nr of residues at the ends unique to s1 (overhanging s2).
*/
static int alnOverhang(char *s1, char *s2)
{
  int overhang = 0,
    s1started = 0,
    s2started = 0;
  char *s1save = s1, *s2save = s2;

  for (; *s1 && *s2; s1++, s2++) 
    {
      if (!isGap(*s1)) s1started = 1;
      if (!isGap(*s2)) s2started = 1;
      if (s1started && !s2started) overhang++;
    }	
  
  s1--; s2--;
  s1started = s2started = 0;
  for (; s1>=s1save && s2>=s2save; s1--, s2--) 
    {
      if (!isGap(*s1)) s1started = 1;
      if (!isGap(*s2)) s2started = 1;
      if (s1started && !s2started) overhang++;
    }	

  /* printf ("%s\n%s\nOverhang=%d\n", s1save, s2save, overhang);*/

  return overhang;
}


/***********************************************************
 *		Removing alignments/columns		   *
 ***********************************************************/

void rmFinaliseColumnRemoval(BelvuContext *bc)
{
  rmGappySeqs(bc, 100.0);
  rmFinalise(bc);
}

/* Removes all columns between from and to, inclusively.
 * Note that from and to are sequence coordinates, 1...maxLen !!!!!  */
void rmColumn(BelvuContext *bc, const int from, const int to)
{
  int
    i=0, j=0, len = to - from + 1;
  ALN 
    *alni=NULL;
  
  fprintf(stderr, "Removing Columns %d-%d.\n", from, to);

  for (i = 0; i < bc->alignArr->len; i++) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
      char *alnSeq = alnGetSeq(alni);

      /* If N or C terminal trim, change the coordinates */
      if (from == 1)
        for (j = from; j <= to;  j++) 
          {
            /* Only count real residues */
            
            if (!isGap(alnSeq[j-1]))
              (alni->start < alni->end ? alni->start++ : alni->start--); 
          }

      if (to == bc->maxLen)
        for (j = from; j <= to;  j++) 
          {
            /* Only count real residues */
            if (!isGap(alnSeq[j-1]))
              (alni->start < alni->end ? alni->end-- : alni->end++); 
          }
	
      /* Remove the columns */
      for (j = 0; j < bc->maxLen-to+1 /* Including terminator 0 */;  j++) 
        {
          alnSeq[from+j-1] = alnSeq[to+j];
	}
      j--;
      
      if (alnSeq[from+j-1] || alnSeq[to+j])
        printf("Still a bug in rmColumn(): End=%c, Oldend=%c\n", 
               alnSeq[from+j-1],
               alnSeq[to+j]);
    }

  bc->maxLen -= len;

  bc->saved = FALSE;
}


/* Remove columns whose conservation is below the given cutoff */
void rmColumnCutoff(BelvuContext *bc, const double from, const double to)
{
  int 
    i, j, max, removed=0, oldmaxLen=bc->maxLen;
  static double cons ;
  
  
  for (i = bc->maxLen-1; i >= 0; i--) 
    {
      if (colorBySimilarity(bc)) 
        {
          cons = bc->conservation[i];
        }
      else 
        {
          max = 0;
          for (j = 1; j < 21; j++)
            {
              if (bc->conservCount[j][i] > max) 
                {
                  max = bc->conservCount[j][i];
                }
            }	
          cons = (double)max / bc->alignArr->len;
        }
      
      if (cons > from && cons <= to) 
        {
          printf("removing %d, cons= %.2f\n", i+1, cons);
          rmColumn(bc, i+1, i+1);
          if (++removed == oldmaxLen) 
            {
              g_critical("You have removed all columns.  Prepare to exit Belvu");
              exit(EXIT_SUCCESS);
            }
        }
    }
  
  bc->saved = FALSE;
  rmFinaliseColumnRemoval(bc);
}


void rmEmptyColumns(BelvuContext *bc, double cutoff)
{
  int
    i=0, j=0, gaps=0, totseq=0, removed=0, oldmaxLen=bc->maxLen;
  ALN 
    *alni=NULL;
  char
    c='\0';

  for (totseq = i = 0; i < bc->alignArr->len; ++i) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
      
      if (!alni->markup) 
        totseq++;
    }

  for (j = bc->maxLen-1; j >= 0; j--) 
    {
      for (gaps = i = 0; i < bc->alignArr->len; i++) 
        {
          alni = &g_array_index(bc->alignArr, ALN, i);
          char *alnSeq = alnGetSeq(alni);
          
          if (!alni->markup) 
            {
              c = alnSeq[j];
              
              if (isGap(c) || c == ' ') 
                gaps++;
	    }
	}
      
      if ((double)gaps/totseq >= cutoff - MACHINE_RES) 
        {
          rmColumn(bc, j+1, j+1);
          if (++removed == oldmaxLen)
            {
              g_critical("You have removed all columns.  Prepare to exit Belvu");
              exit(0);
            }
	}
    }
}


/* Certain actions need to be enabled/disabled depending on certain properties */
void greyOutInvalidActionsForGroup(BelvuContext *bc, GtkActionGroup *action_group)
{
  if (!action_group)
    return;
  
  enableMenuAction(action_group, "Output", bc->displayScores);
  enableMenuAction(action_group, "rmScore", bc->displayScores);
  enableMenuAction(action_group, "scoreSort", bc->displayScores);
  
  enableMenuAction(action_group, "toggleColorByResId", !colorByConservation(bc));
  enableMenuAction(action_group, "colorByResId", !colorByConservation(bc));
  enableMenuAction(action_group, "saveColorScheme", !colorByConservation(bc));
  enableMenuAction(action_group, "loadColorScheme", !colorByConservation(bc));
  
  enableMenuAction(action_group, "colorSchemeCustom", bc->haveCustomColors);
  
  /* ignoreGaps used to be greyed out in old belvu when in color-by-residue
   * mode; it does affect the display when in that mode if using a %ID
   * threshold, though, so it should probably be enabled (or alternatively the 
   * logic should be adjusted so that it does not affect the residue colors, if 
   * that was the intent...) */
  /* enableMenuAction(action_group, "ignoreGaps", colorByConservation(bc)); */
  enableMenuAction(action_group, "printColors", colorByConservation(bc));
  
  enableMenuAction(action_group, "excludeHighlighted", bc->selectedAln != NULL);
  
  enableMenuAction(action_group, "RecalcTree", bc->treeHead == NULL);
}


/* Certain actions need to be enabled/disabled depending on certain properties */
void greyOutInvalidActions(BelvuContext *bc)
{
  /* Need to do the action groups for the main window and the tree (wrapped
   * windows currently don't have anything that gets greyed out but if they
   * do in future they will need updating here too) */
  
  GtkActionGroup *actionGroup = belvuWindowGetActionGroup(bc->belvuWindow);
  greyOutInvalidActionsForGroup(bc, actionGroup);

  actionGroup = belvuTreeGetActionGroup(bc->belvuTree);
  greyOutInvalidActionsForGroup(bc, actionGroup);
}


/* Set the head node of the main tree (null to reset) */
void setTreeHead(BelvuContext *bc, TreeNode *headNode)
{
  bc->treeHead = headNode;
  
  /* Whether a tree exists affects whether some menu items are greyed out */
  greyOutInvalidActions(bc);
}


/* This function should be called after removing sequences or columns */
static void rmFinalise(BelvuContext *bc) 
{
  /*    ruler[maxLen] = 0;*/
  checkAlignment(bc);
  setConsSchemeColors(bc);
  
  /* Removing seqs/cols invalidates the tree, so set the tree head to NULL. */
  setTreeHead(bc, NULL);
  
  /* Removing seqs/cols invalidates the conservation plot, so recalculate it */
  belvuConsPlotRecalcAll(bc->consPlot);
}


/* Get rid of seqs that start or end with a gap.
 */
void rmPartialSeqs(BelvuContext *bc)
{
  int i=0, n=0;
  ALN *alni = NULL;

  for (i = 0; i < bc->alignArr->len; ) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
      char *alnSeq = alnGetSeq(alni);
      
      if (isGap(alnSeq[0]) || isGap(alnSeq[bc->maxLen-1])) 
        {
          /* Remove entry */
          n++;

          if (bc->selectedAln == alni) 
            bc->selectedAln = 0;
          
          g_array_remove_index(bc->alignArr, i);
          bc->saved = 0;
	}
      else 
        {
          ++i;
        }
    }

  g_message("%d partial sequences removed.  %d seqs left.\n\n", n, bc->alignArr->len);

  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Get rid of seqs that are too gappy.
 */
void rmGappySeqs(BelvuContext *bc, const double cutoff)
{
  int i=0, j=0, n=0, gaps;
  ALN *alni = NULL;

  for (i = 0; i < bc->alignArr->len; ) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
      char *alnSeq = alnGetSeq(alni);

      for (gaps = j = 0; j < bc->maxLen; j++)
        if (isGap(alnSeq[j])) 
          gaps++;

      if ((double)gaps/bc->maxLen >= cutoff/100.0) 
        {
          /* Remove entry */
          n++;
          
          if (bc->selectedAln == alni) 
            bc->selectedAln = 0;
          
          g_array_remove_index(bc->alignArr, i);
          bc->saved = 0;
	}
      else 
	{
	  i++;
	}
    }

    g_message("%d gappy sequences removed.  %d seqs left.\n\n", n, bc->alignArr->len);

    arrayOrder(bc->alignArr);
}


/* Remove empty (gappy) columns if the 'remove empty columns' option
 * is enabled. */
void rmFinaliseGapRemoval(BelvuContext *bc)
{
  if (bc->rmEmptyColumnsOn) 
    rmEmptyColumns(bc, 1.0);
  
  rmFinalise(bc);
}


/* Get rid of seqs that are more than x% identical with others. 
 * Keep the  first one.
 */
void mkNonRedundant(BelvuContext *bc, const double cutoff)
{
  int i=0,j=0, n=0, overhang;
  ALN *alni=NULL, *alnj=NULL;
  double id = 0.0;

  for (i = 0; i < bc->alignArr->len; i++) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
      char *alniSeq = alnGetSeq(alni);

      for (j = 0; j < bc->alignArr->len; j++) 
        {
          if (i == j) 
            continue;

          alnj = &g_array_index(bc->alignArr, ALN, j);
          char *alnjSeq = alnGetSeq(alnj);
	    
          overhang = alnOverhang(alnjSeq, alniSeq);
          id = identity(alniSeq, alnjSeq, bc->penalize_gaps);

          if (!overhang && id > cutoff)
	    {
              g_message("%s/%d-%d and %s/%d-%d are %.1f%% identical. "
			"The first includes the latter which was removed.\n",
			alni->name, alni->start, alni->end,
			alnj->name, alnj->start, alnj->end,
			id);
		
              /* Remove entry j */
              n++;
		
              if (bc->selectedAln == alnj) 
                bc->selectedAln = NULL;
              
              g_array_remove_index(bc->alignArr, j);
              bc->saved = FALSE;
              
              if (j < i) 
                {
                  i--;
                  alni = &g_array_index(bc->alignArr, ALN, i);
		}
              j--;
	    }
	}
    }

  g_message("%d sequences removed at the %.0f%% level.  %d seqs left.\n\n", n, cutoff, bc->alignArr->len);
  
  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Get rid of seqs that are less than x% identical with any of the others. 
 */
void rmOutliers(BelvuContext *bc, const double cutoff)
{
  int i=0,j=0, n=0;
  ALN *alni=NULL, *alnj=NULL;
  
  for (i = 0; i < bc->alignArr->len - 1; ) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);
    
      double maxid = 0;
      for (maxid=0, j = 0; j < bc->alignArr->len; ++j) 
	{
	  if (i == j) 
	    continue;
	
	  alnj = &g_array_index(bc->alignArr, ALN, j);
	  double id = identity(alnGetSeq(alni), alnGetSeq(alnj), bc->penalize_gaps);
	
	  if (id > maxid) 
	    maxid = id;
	}
    
      if (maxid < cutoff) 
	{
	  g_message("%s/%d-%d was max %.1f%% identical to any other sequence and was removed.\n",
		    alni->name, alni->start, alni->end, maxid);
	
	  /* Remove entry */
	  n++;
	
	  if (bc->selectedAln == alni) 
	    bc->selectedAln = NULL;
	
	  g_array_remove_index(bc->alignArr, i);
	  bc->saved = FALSE;
	}
      else 
	{
	  i++;
	}
    }
  
  g_message("%d sequences removed at the %.0f%% level.  %d seqs left.\n\n", n, cutoff, bc->alignArr->len);
  
  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Remove sequences that have a score below the given value */
void rmScore(BelvuContext *bc, const double cutoff)
{
  scoreSort(bc);
  
  /* Save bc->selectedAln */
  ALN aln;
  initAln(&aln);
  
  if (bc->selectedAln)
    alncpy(&aln, bc->selectedAln);
  
  int numRemoved = 0;
  int i = 0;
  
  for (i = 0; i < bc->alignArr->len; ) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
    
      if (alnp->score < cutoff) 
	{ 
	  ++numRemoved;
	  
	  g_message("Removing %s/%d-%d (score %.1f)\n",
		    alnp->name, alnp->start, alnp->end, alnp->score);
	
	  if (bc->selectedAln == alnp) 
	    bc->selectedAln = NULL;
	  
	  g_array_remove_index(bc->alignArr, i);
	  bc->saved = FALSE;
	}
      else
	{
	  i++;
	}
    }
  
  arrayOrder(bc->alignArr);
  
  g_message("%d sequences with score < %.1f removed.  %d seqs left.\n\n", numRemoved, cutoff, bc->alignArr->len);
  
  /* Find bc->selectedAln in new array */
  if (bc->selectedAln) 
    { 
      int ip = 0;
      ALN aln;
      initAln(&aln);
      
      if (!arrayFind(bc->alignArr, &aln, &ip, (void*)scoreorder)) 
	{
	  bc->selectedAln = NULL;
	}
      else
	{
	  bc->selectedAln = &g_array_index(bc->alignArr, ALN, ip);
	}
  }
  
  bc->alignYStart = 0;
  
  rmFinaliseGapRemoval(bc);
}


/***********************************************************
 *		      Read/write files 			   *
 ***********************************************************/

/* Save the given tree to the given file in New Hampshire format */
void saveTreeNH(TreeNode *headNode, TreeNode *node, FILE *file)
{
  if (!node) 
    return;
  
  if (node->left && node->right) 
    {
      fprintf(file, "(\n");
      saveTreeNH(headNode, node->left, file);
      fprintf(file, ",\n");
      saveTreeNH(headNode, node->right, file);
      fprintf(file, ")\n");
      
      if (node != headNode)	/* Not exactly sure why this is necessary, but njplot crashes otherwise */
        fprintf(file, "%.0f", node->boot+0.5);
    }
  else
    {
      fprintf(file, "%s", node->name);
    }
  
  if (node != headNode)	/* Not exactly sure why this is necessary, but njplot crashes otherwise */
    fprintf(file, ":%.3f", node->branchlen/100.0);
}


void writeMSF(BelvuContext *bc, FILE *pipe) /* c = separator between name and coordinates */
{
  int i=0, j=0;
  int maxfullnamelen = bc->maxNameLen + 1 + bc->maxStartLen + 1 + bc->maxEndLen;
  int paragraph=0, alnstart=0, alnend=0, alnlen=0, linelen=50, blocklen=10;
  
  if (bc->saveCoordsOn)
    maxfullnamelen = bc->maxNameLen + 1 + bc->maxStartLen + 1 + bc->maxEndLen;
  else
    maxfullnamelen = bc->maxNameLen;
  
  /* Title */
  fprintf(pipe, "PileUp\n\n %s  MSF: %d  Type: X  Check: %d  ..\n\n",
          bc->Title, bc->maxLen, GCGgrandchecksum(bc));

  /* Names and checksums */
  for(i=0; i < bc->alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (alnp->markup) 
        continue;

      char *tmpStr = bc->saveCoordsOn 
        ? blxprintf("%s%c%d-%d", alnp->name, bc->saveSeparator, alnp->start, alnp->end)
        : blxprintf("%s", alnp->name);

      fprintf(pipe, "  Name: %-*s  Len:  %5d  Check:  %5d  Weight: %.4f\n",
              maxfullnamelen,
              tmpStr,
              bc->maxLen,
              GCGchecksum(bc, alnGetSeq(alnp)),
              1.0);

      g_free(tmpStr);
    }
  
  fprintf(pipe, "\n//\n\n");

  /* Alignment */
  while (paragraph * linelen < bc->maxLen)
    {
      for (i = 0; i < bc->alignArr->len; ++i)
	{
          alnstart = paragraph * linelen; 
          alnlen = ( (paragraph+1)*linelen < bc->maxLen ? linelen : bc->maxLen - alnstart );
          alnend = alnstart + alnlen;

          ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
          
          if (alnp->markup) 
            continue;

          char *tmpStr = bc->saveCoordsOn 
            ? blxprintf("%s%c%d-%d", alnp->name, bc->saveSeparator, alnp->start, alnp->end)
            : blxprintf("%s", alnp->name);
          
          fprintf(pipe, "%-*s  ", maxfullnamelen, tmpStr);
          g_free(tmpStr);
          
          for (j = alnstart; j < alnend; ) 
            {
              char *alnpSeq = alnGetSeq(alnp);
              fprintf(pipe, "%c", alnpSeq[j]);
              j++;
              
              if ( !((j-alnstart) % blocklen) ) 
                fprintf(pipe, " ");
	    }
          
          fprintf(pipe, "\n");
	}
      
      fprintf(pipe, "\n");
      paragraph++;
    }

  fclose(pipe);
  fflush(pipe);
  
  bc->saved = TRUE;
}


static void readMSF(BelvuContext *bc, FILE *pipe)
{
  char seq[1001], *cp=NULL, *cq=NULL;
  char line[MAXLENGTH + 1];
  line[0] = 0;
  seq[0] = 0;
  
  ALN aln;
  initAln(&aln);
  
  fprintf(stderr, "\nDetected MSF format\n");

  /* Read sequence names */
  while (!feof (pipe))
    { 
      if (!fgets (line, MAXLENGTH, pipe)) 
        break;

      if (!strncmp(line, "//", 2)) 
        {
          break;
        }
      else if (strstr(line, "Name:") && (cp = strstr(line, "Len:")) && strstr(line, "Check:")) 
	{
          int len;

          sscanf(cp+4, "%d", &len);

          if (bc->maxLen < len) 
            bc->maxLen = len;

          cp = strstr(line, "Name:")+6;
          
          while (*cp && *cp == ' ') 
            cp++;

          parseMulLine(bc, cp, &aln);

          aln.sequenceStr = g_string_new("");
          aln.nr = bc->alignArr->len + 1;

          g_array_append_val(bc->alignArr, aln);
          g_array_sort(bc->alignArr, alphaorder);
	}
    }

  /* Read sequence alignment */
  while (!feof (pipe))
    { 
      if (!fgets (line, MAXLENGTH, pipe)) 
        break;

      cp = line;
      cq = seq;
      while (*cp && *cp == ' ') cp++;
      while (*cp && *cp != ' ') cp++; /* Spin to sequence */
      
      while (*cp && cq-seq < 1000) 
        {
          if (isAlign(*cp)) *cq++ = *cp;
          cp++;
	}
	
      *cq = 0;

      if (*seq) 
        {
          cp = line;
          while (*cp && *cp == ' ') 
            cp++;

          parseMulLine(bc, cp, &aln);

          int ip = 0;
          if (arrayFind(bc->alignArr, &aln, &ip, (void*)alphaorder)) 
            {
              ALN *alnp = &g_array_index(bc->alignArr, ALN, ip);
              g_string_append(alnp->sequenceStr, seq);
	    }
          else
            {
              g_critical("Cannot find back %s %d %d seq=%s.\n", 
                         aln.name, aln.start, aln.end, seq);
	    }
	}
    }
  
  bc->saveFormat = BELVU_FILE_MSF;
}


/* Parse annotation lines from #=GS lines in the Mul file format */
static void parseMulAnnotationLine(BelvuContext *bc, const char *seqLine)
{
  const char *cp = seqLine;

  char
    *namep,		/* Seqname */
    name[MAXNAMESIZE+1],/* Seqname */
    *labelp,		/* Label (OS) */
    *valuep;		/* Organism (c. elegans) */

  
  /* Ignore anything that's not a 'GS' line */
  if (strncmp(cp, "#=GS", 4) != 0)
    return;
      
  /* Copy the portion after  the '#=GS' into 'name' */
  strcpy(name, cp+4);
  namep = name;
  
  /* Ignore leading whitespace */
  while (*namep == ' ') 
    namep++;
      
  /* Find the end of the name */
  labelp = namep;
  while (*labelp != ' ') 
    labelp++;
  
  *labelp = 0;		/* Terminate namep */
  ++labelp;		/* Walk across terminator */
  
  while (*labelp == ' ') /* Find the start of the label */
    ++labelp;
      
  valuep = labelp;
  while (*valuep != ' ') valuep++; /* Find the end of the label */
  while (*valuep == ' ') valuep++; /* Find the start of the value */
  
  if (strncasecmp(labelp, "LO", 2) == 0)
    {
      int i, colornr;
      
      /* Check the value to see if it matches one of our colour names */
      for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
        {
          if (strcasecmp(colorNames[i], valuep) == 0) 
            {
              colornr = i;
              break;
            }
        }
      
      if (colornr == -1)
        {
          printf("Unrecognized color: %s\n", valuep);
          colornr = 0;
        }
          
      /* Create an ALN struct from the name */
      ALN aln;
      initAln(&aln);
      str2aln(bc, namep, &aln);

      /* Find the corresponding sequence */
      if (!arrayFind(bc->alignArr, &aln, &i, (void*)alphaorder))
        {
          g_critical("Cannot find '%s' [%d %d] in alignment.\n", aln.name, aln.start, aln.end);
          return;
        }
      
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      alnp->color = colornr;
    }
  else if (strncasecmp(labelp, bc->organismLabel, 2) == 0)
    {
      /* Add organism to table of organisms.  This is necessary to make all 
         sequences of the same organism point to the same place and to make a 
         non-redundant list of present organisms */
      
      /* Find string in permanent stack */
      if (!(valuep = strstr(cp, valuep)))
        {
          g_critical("Failed to parse organism properly.\n");
          return;
        }
      
      /* Create an ALN struct from this organism */
      ALN aln;
      initAln(&aln);
      
      aln.organism = valuep; /* Find organism string in permanent stack */
      g_array_append_val(bc->organismArr, aln);
      g_array_sort(bc->organismArr, organism_order);
          
      if (strchr(cp, '/') && strchr(cp, '-'))
        {
          str2aln(bc, namep, &aln);
              
          /* Find the corresponding sequence */
          int ip = 0;
          
          if (!arrayFind(bc->alignArr, &aln, &ip, (void*)alphaorder)) 
            {
              g_critical("Cannot find '%s' [%d %d] in alignment.\n", aln.name, aln.start, aln.end);
              return;
            }
             
          ALN *alnp = &g_array_index(bc->alignArr, ALN, ip);
              
          /* Store pointer to unique organism in ALN struct */
          aln.organism = valuep;
          ip = 0;
          
          arrayFind(bc->organismArr, &aln, &ip, (void*)organism_order);
          alnp->organism = g_array_index(bc->organismArr, ALN, ip).organism;
        }
      else 
        {
          /* Organism specified for entire sequences.
             Find all ALN instances of this sequences.
          */
          int i = 0;
          for (i = 0; i < bc->alignArr->len; ++i) 
            {
              ALN *alnp = &g_array_index(bc->alignArr, ALN, i);

              if (strcmp(alnp->name, namep) == 0) 
                {
                  /* Store pointer to unique organism in ALN struct */
                  aln.organism = valuep;
                  int ip = 0;

                  arrayFind(bc->organismArr, &aln, &ip, (void*)organism_order);
                  alnp->organism = g_array_index(bc->organismArr, ALN, ip).organism;
                }
            }
        }
    }
}


/* Add sequence data to an alignment from a line of text of the format:
 * SEQ_NAME     SEQUENCE_DATA,
 * where alnstart gives the position of the start of SEQUENCE_DATA. */
static void appendSequenceDataToAln(BelvuContext *bc, char *line, const int alnstart)
{
  /* Find the alignment name */
  ALN aln;
  initAln(&aln);
  parseMulLine(bc, line, &aln);

  /* See if this alignment is in the alignments array */
  int ip = 0;
  if (arrayFind(bc->alignArr, &aln, &ip, (void*)alphaorder))
    {
      /* Append this bit of sequence to the existing alignment */
      ALN *alnp = &g_array_index(bc->alignArr, ALN, ip);
      g_string_append(alnp->sequenceStr, line + alnstart);
      
      /* Recalculate the max alignment length */
      if (alnp->sequenceStr->len > bc->maxLen)
        bc->maxLen = alnp->sequenceStr->len;
    }
  else
    {
      /* Create a new alignment and add this bit of sequence */
      aln.sequenceStr = g_string_new(line + alnstart);
      aln.nr = bc->alignArr->len + 1;
      
      g_array_append_val(bc->alignArr, aln);
      g_array_sort(bc->alignArr, alphaorder);
      
      /* Recalculate the max alignment length */
      if (aln.sequenceStr->len > bc->maxLen)
        bc->maxLen = aln.sequenceStr->len;
    }
}


/* ReadMul 
 * parses an alignment file in mul (stockholm) or selex format and puts it in the Align array
 *
 * Assumes header contains ' ' or '#'
 *
 * Alignment can have one of the following formats:
 *  CSW_DROME  VTHIKIQNNGDFFDLYGGEKFATLP
 *  CSW_DROME  51   75   VTHIKIQNNGDFFDLYGGEKFATLP
 *  CSW_DROME  51   75   VTHIKIQNNGDFFDLYGGEKFATLP P29349
 *  KFES_MOUSE/458-539    .........WYHGAIPW.....AEVAELLT........HTGDFLVRESQG
 *
 */
static void readMul(BelvuContext *bc, FILE *pipe)
{
  char line[MAXLENGTH+1];
  line[0] = 0;

  /* Read raw alignment into stack
   *******************************/
  
  int alnstart = MAXLENGTH;
  GSList *alnList = NULL;
  
  while (!feof (pipe))
    { 
      /* EOF checking to make acedb calling work */
      if (!fgets (line, MAXLENGTH, pipe) || (unsigned char)*line == (unsigned char)EOF ) 
        break;
      
      if (!strncmp(line, "PileUp", 6)) 
        {
          readMSF(bc, pipe);
          return;
        }

      /* Remove any trailing newline */
      char *cp = strchr(line, '\n');
      if (cp) 
        *cp = 0;
      
      int lineLen = strlen(line);

      if (lineLen > 0 && *line != '#' && strcmp(line, "//") != 0)
	{
          /* Sequence line */
          char *cq = strchr(line, ' ');
          
          if (!cq) 
            g_error("Error reading selex file; no spacer between name and sequence in the following line:\n%s", line);
          
          /* Find which column the alignment starts in */
          int i = 0;
          for (i = 0; line[i] && line[i] != ' '; ++i);
          for (; line[i] && !isAlign(line[i]); ++i);
          
          /* Remember the leftmost start position of any alignment. We'll assume
           * all alignments start in the same column. */
          if (i < alnstart) 
            alnstart = i;

          /* Remove optional accession number at end of alignment */
          /* FOR PRODOM STYLE ALIGNMENTS W/ ACCESSION TO THE RIGHT - MAYBE MAKE OPTIONAL
           * This way it's incompatible with alignments with ' ' gapcharacters */
          for (; line[i] && isAlign(line[i]); i++);
          line[i] = 0;

          /* Store the line for processing later (once alnstart has been calculated) */
          alnList = g_slist_append(alnList, g_strdup(line));
	}
      else if (!strncmp(line, "#=GF ", 5) || 
               !strncmp(line, "#=GS ", 5)) 
        {
	  /* Store all annotation lines in a list. Prepend the items because that
	   * is more efficient, and then reverse the list at the end */
	  bc->annotationList = g_slist_prepend(bc->annotationList, g_strdup(line));
          parseMulAnnotationLine(bc, line);
        }
      else if (!strncmp(line, "#=GC ", 5) || 
               !strncmp(line, "#=GR ", 5) || 
               !strncmp(line, "#=RF ", 5)) 
        {
          /* These are markup lines that are shown in the alignment list */
          alnList = g_slist_append(alnList, g_strdup(line));
        }
      else if (!strncmp(line, "# matchFooter", 13)) 
        {
          /* Match Footer  */
          bc->matchFooter = TRUE;
          break;
        }
    }
  
  /* Reverse the annotation list, because we prepended items instead of appending them */
  bc->annotationList = g_slist_reverse(bc->annotationList);
  
  /* Loop through all of the alignment lines and extract the sequence string */
  GSList *alnItem = alnList;
  for ( ; alnItem; alnItem = alnItem->next)
    {
      char *alnStr = (char*)(alnItem->data);
      appendSequenceDataToAln(bc, alnStr, alnstart);
      
      /* Free the line string, now we're done with it */
      g_free(alnStr);
      alnItem->data = NULL;
    }

  g_slist_free(alnList);
  alnList = NULL;
  
/* For debugging * /
   for (i = 0; i < nseq; i++) {
   alnp = arrp(Align, i, ALN);
   printf("\n%-10s %4d %4d %s %d %s\n", 
   alnp->name, alnp->start, alnp->end, alnp->seq, 
   alnp->len, alnp->organism);
   }
   for (i=0; i < arrayMax(organismArr); i++) 
   printf("%s\n", arrp(organismArr, i, ALN)->organism);
   */
  
  if (bc->alignArr->len == 0 || bc->maxLen == 0) 
    g_error("Unable to read sequence data");
  
  bc->saveFormat = BELVU_FILE_MUL;
}


/* 
 Format the name/start-end string 
 
 For convenience, used by writeMul. Returns a newly allocated string
 which must be freed with g_free
 */
static char *writeMulName(BelvuContext *bc, ALN *aln) 
{  
  char *name = NULL; /* result */
  
  char *cp, *namep, GRname[MAXNAMESIZE*2+2], GRfeat[MAXNAMESIZE+1];
  
  if (aln->markup == GC) 
    {
      name = blxprintf("#=GC %s", aln->name);
      return name;
    }
  
  /* NOTE: GR lines have the feature concatenated in aln->name - must separate */
  if (aln->markup == GR) 
    {
      namep = GRname;
      strncpy(namep, aln->name, 50);
      cp = strchr(namep, ' ');
      strncpy(GRfeat, cp+1, 50);
      *cp = 0;
    }
  else
    {
      namep = aln->name;
    }
  
  if (!bc->saveCoordsOn) 
    {
      name = blxprintf("%s", namep);
    }
  else 
    {
      name = blxprintf("%s%c%d-%d", namep, bc->saveSeparator, aln->start, aln->end);
    }
  
  if (aln->markup == GR)
    name = blxprintf("#=GR %s %s", name, GRfeat);
  
  return name;
}


void writeMul(BelvuContext *bc, FILE *fil)
{
  /* Write Annotation */
  GSList *annotationItem = bc->annotationList;
  
  for ( ; annotationItem; annotationItem = annotationItem->next)
    {
      const char *line = (const char*)(annotationItem->data);
      fprintf(fil, "%s\n", line);
    }

  /* Find max width of name column */
  int i = 0;
  int maxWidth = 0;
  
  for ( ; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);

      char *mulName = writeMulName(bc, alnp);
      int curWidth = strlen(mulName);
      g_free(mulName);
    
      if ( curWidth > maxWidth) 
	maxWidth = curWidth;
    }
    
  /* Write alignment */
  for (i = 0; i < bc->alignArr->len; i++)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      char *mulName = writeMulName(bc, alnp);
      const char *alnpSeq = alnGetSeq(alnp);
    
      fprintf(fil, "%-*s %s\n", maxWidth, mulName, (alnpSeq ? alnpSeq : ""));
    
      g_free(mulName);
    }
  
  fprintf(fil, "//\n");
  
  fclose(fil);
  fflush(fil);

  bc->saved = TRUE;
}


/* readFile 
 * Determines the format of the input file and calls the appropriate 
 * parser. Valid file formats are fasta, MSF, mul (stockholm) or selex.
 */
void readFile(BelvuContext *bc, FILE *pipe)
{
  char   ch = '\0';
  char line[MAXLENGTH+1];
  line[0] = 0;
  
  /* Parse header to check for MSF or Fasta format */
  while (!feof (pipe))
    { 
      ch = fgetc(pipe);
      
      if (!isspace(ch)) 
        {
          if (ch == '>') 
            {
              ungetc(ch, pipe);
              return readFastaAln(bc, pipe);
            }
          else
            {
              break;
            }
        }
      else if (ch == '\n')
        {
          if (!fgets (line, MAXLENGTH, pipe)) 
            break;
        }
      
      if (strstr(line, "MSF:") && strstr(line, "Type:") && strstr(line, "Check:")) 
        {
          return readMSF(bc, pipe);
        }
    }
  
  if (!feof(pipe))
    ungetc(ch, pipe);

  return readMul(bc, pipe);
}


void writeFasta(BelvuContext *bc, FILE *pipe)
{
  int i=0, n=0;
  char *cp = NULL;

  for (i = 0; i < bc->alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (alnp->markup) 
        continue;
	
      if (bc->saveCoordsOn)
        fprintf(pipe, ">%s%c%d-%d\n", alnp->name, bc->saveSeparator, alnp->start, alnp->end);
      else
        fprintf(pipe, ">%s\n", alnp->name);
	
      for (n=0, cp = alnGetSeq(alnp); *cp; cp++)
        {
	  if (bc->saveFormat == BELVU_FILE_ALIGNED_FASTA)
            {
              fputc(*cp, pipe);
              n++;
	    }
          else
            {
              if (!isGap(*cp)) 
                {
                  fputc(*cp, pipe);
                  n++;
		}
	    }
	    
          if (n && !(n % 80) ) 
            {
              fputc('\n', pipe);
              n = 0;
	    }
	}
	
      if (n)
        fputc('\n', pipe);
    }

  fclose(pipe);
  fflush(pipe);
  
  bc->saved = TRUE;
}


static int getMatchStates(BelvuContext *bc)
{
  int i=0, j=0, n=0, retval=0;

  for (j = 0; j < bc->maxLen; ++j) 
    {
      n = 0;
      for (i = 0; i < bc->alignArr->len; ++i) 
        {
          ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
          char *alnpSeq = alnGetSeq(alnp);
          
          if (isalpha(alnpSeq[j])) 
            n++;
	}
      
      if (n > bc->alignArr->len / 2) 
        retval++;
    }

  return retval;
}


/* Output a.a. probabilities in HMMER format.

   Disabled counting of match-state residues only (script HMM-random does that already)

   Disabled pseudocounts - not sure what the best way is.
*/
void outputProbs(BelvuContext *bc, FILE *fil)
{
  /* From swissprot 33:

     Ala (A) 7.54   Gln (Q) 4.02   Leu (L) 9.31   Ser (S) 7.19
     Arg (R) 5.15   Glu (E) 6.31   Lys (K) 5.94   Thr (T) 5.76
     Asn (N) 4.54   Gly (G) 6.86   Met (M) 2.36   Trp (W) 1.26
     Asp (D) 5.29   His (H) 2.23   Phe (F) 4.06   Tyr (Y) 3.21
     Cys (C) 1.70   Ile (I) 5.72   Pro (P) 4.91   Val (V) 6.52
  */
  
  double f_sean[21] = {
    0.0, 
    .08713, .03347, .04687, .04953, .03977, .08861, .03362, .03689, .08048, .08536, 
    .01475, .04043, .05068, .03826, .04090, .06958, .05854, .06472, .01049, .02992
  };

  double p = 0.0;
  int i = 0, c[21], n=0, nmat = 0;
  char *cp = NULL;

  for (i = 1; i <= 20; i++) 
    c[i] = 0;

  for (i = 0; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      cp = alnGetSeq(alnp);

      for (; *cp; cp++)
        {
          if (a2b_sean[(unsigned char)(*cp)])
            {
              c[a2b_sean[(unsigned char)(*cp)]]++;
              n++;
            }
        }
    }

  if (!n)
    g_error("No residues found");

  if (0) 
    {
      nmat = getMatchStates(bc);	/* Approximation, HMM may differ slightly */

      if (!nmat) 
        g_error("No match state columns found");
      
      printf("Amino\n");
      
      for (i = 1; i <= 20; i++) 
        {
          /* One way:  p = ((double)c[i]/nmat + 20*f_sean[i]) / (double)(n/nmat+20);*/
          /* Other way: */
          p = ((double)c[i] + 20*f_sean[i]) / (double)(n+20);

          if (p < 0.000001) 
            p = 0.000001; /* Must not be zero */

          printf("%f ", p);
          /* printf("   %d %f %d \n", c[i], f[i], n);*/
        }
    }
  else
    {
      printf("Amino\n");
      for (i = 1; i <= 20; ++i) 
        {
          p = (double)c[i] / (double)n;
          if (p < 0.000001) p = 0.000001; /* Must not be zero */
          printf("%f ", p);
	}
    }
    
  printf("\n");
  
  fclose(fil);
  fflush(fil);

  bc->saved = TRUE;
}



/* Calculate the identity of each sequence against each other, and print the
 * results to stdout. */
void listIdentity(BelvuContext *bc)
{
  /* This may take some time, so feed back to the user what is happening. */
  g_message("Outputting identities...\n");
  setBusyCursor(bc, TRUE);

  int i=0,j=0,n=0 ;
  double totsc=0, maxsc=0, minsc=1000000,
         totid=0.0, maxid=0.0, minid=100.0;
  
  for (i = n = 0; i < bc->alignArr->len - 1; ++i) 
    {
      /* Update the display periodically, otherwise the busy cursor might not
       * get updated (e.g. if we were outside the main window when we set the 
       * cursor and the user later enters the window). */
      while (gtk_events_pending())
        gtk_main_iteration();

      ALN *alni = &g_array_index(bc->alignArr, ALN, i);
      
      if (alni->markup > 0) /* ignore markup lines */
        continue;

      for (j = i+1; j < bc->alignArr->len; ++j) 
        {
          ALN *alnj = &g_array_index(bc->alignArr, ALN, j);

          if (alnj->markup > 0) /* ignore markup lines */
            continue;

          double id = identity(alnGetSeq(alni), alnGetSeq(alnj), bc->penalize_gaps);
          totid += id;

          if (id > maxid) 
            maxid = id;
          
          if (id < minid) 
            minid = id;
      
          double sc = score(alnGetSeq(alni), alnGetSeq(alnj), bc->penalize_gaps);
          totsc += sc;
          
          if (sc > maxsc) 
            maxsc = sc;
          
          if (sc < minsc) 
            minsc = sc;
          
          printf("%s/%d-%d and %s/%d-%d are %.1f%% identical, score=%f\n",
                 alni->name, alni->start, alni->end,
                 alnj->name, alnj->start, alnj->end,
                 id, sc);
          
          ++n;
        }
      printf("\n");
    }

  printf("Maximum %%id was: %.1f\n", maxid);
  printf("Minimum %%id was: %.1f\n", minid);
  printf("Mean    %%id was: %.1f\n", totid/n);
  printf("Maximum score was: %.1f\n", maxsc);
  printf("Minimum score was: %.1f\n", minsc);
  printf("Mean    score was: %.1f\n", (double)totsc/n);

  setBusyCursor(bc, FALSE);
  g_message("Finished outputting identities.\n");
}


static void onDestroySpawnedWindow(GtkWidget *window, gpointer data)
{
  /* We must remove the window from the list of spawned windows */
  BelvuContext *bc = (BelvuContext*)data;
  bc->spawnedWindows = g_slist_remove(bc->spawnedWindows, window);
}


void fetchAln(BelvuContext *bc, ALN *alnp)
{
  GError *error = NULL;

  if (bc->useWWWFetch)
    {
      /* Get the URL from the enviroment var (or use the hard-coded default
       * value) */
      char  *env, url[1025]="";
      
      if ((env = getenv(FETCH_URL_ENV_VAR)) )
	strcpy(url, env);
      else
	strcpy(url, DEFAULT_FETCH_URL);

      /* Add the sequence name to the url. The URL should be a format string
       * containing '%s' for the name. */
      char *cp = strchr(url, '%');
      
      if (cp)
        ++cp;
      
      if (!cp || *cp != 's')
        {
          g_critical("Invalid URL string %s.\nThe URL must contain the search string \"%%s\", which will be replaced with the sequence name.", url);
        }
      else
        {
          char *link = blxprintf(url, alnp->name);
          
          g_message("Opening URL: %s\n", link);
          seqtoolsLaunchWebBrowser(link, &error);
          
          g_free(link);
        }
    }
  else
    {
      char  *env, fetchProg[1025]="";
      
      if ((env = getenv(FETCH_PROG_ENV_VAR)) )
	strcpy(fetchProg, env);
      else
	strcpy(fetchProg, DEFAULT_FETCH_PROG);

      char *cmd = blxprintf("%s '%s' &", fetchProg, alnp->name);
      
      GtkWidget *pfetchWin = externalCommand(cmd, BELVU_TITLE, bc->belvuAlignment, &error);
      
      if (pfetchWin)
        {
          /* Add the window to our list of spawned windows */
          bc->spawnedWindows = g_slist_prepend(bc->spawnedWindows, pfetchWin);
          g_signal_connect(G_OBJECT(pfetchWin), "destroy", G_CALLBACK(onDestroySpawnedWindow), bc);
        }
      
      g_free(cmd);
    }
  
  reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);

}



/* Utility to return true if the given alignment is highlighted (i.e. has the
 * same name as the selected alignment) */
gboolean alignmentHighlighted(BelvuContext *bc, ALN *alnp)
{
  /* Return true if this alignment has the same name as the selected alignment */
  return (bc->selectedAln && stringsEqual(alnp->name, bc->selectedAln->name, TRUE));
}


/***********************************************************
 *                        Cursors                          *
 ***********************************************************/

/* Utility to set/unset the busy cursor on all open toplevel windows */
void setBusyCursor(BelvuContext *bc, const gboolean busy)
{
  GdkCursor *cursor = bc->defaultCursor;

  if (busy)
    cursor = bc->busyCursor;
  else if (bc->removingSeqs)
    cursor = bc->removeSeqsCursor;

  if (bc->belvuWindow)
    gdk_window_set_cursor(bc->belvuWindow->window, cursor);

  if (bc->belvuTree)
    gdk_window_set_cursor(bc->belvuTree->window, cursor);

  if (bc->consPlot)
    gdk_window_set_cursor(bc->consPlot->window, cursor);

  if (bc->orgsWindow)
    gdk_window_set_cursor(bc->orgsWindow->window, cursor);

  /* Force cursor to change immediately */
  while (gtk_events_pending())
    gtk_main_iteration();
}

