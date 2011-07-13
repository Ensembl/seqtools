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


/*  Description of color_by_similarity algorithm in setConservColors():

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

#include <stdarg.h>
/*#include <stdlib.h> / * Needed for RAND_MAX but clashes with other stuff */
#include <sys/types.h>
#include <sys/time.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <math.h>



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




static void                    fillParents(TreeNode *parent, TreeNode *node);










#ifdef OLD_BELVU_CODE
#include <wh/regular.h>
#include <wh/aceio.h>
#include <wh/graph.h>
#include <wh/gex.h>
#include <wh/key.h>
#define BELVU						    /* Fix horrible include mess in dotter_.h */
#include <wh/dotter_.h>
#include <wh/menu.h>


#define MAXLINE   512
#define MAXLENGTH 100000
#define DBL_CLK_DELAY 2 /* seconds */
#define HSCRWID    1.0  /* Width of Horizontal scroll bar */
#define SCRBACKCOLOR PALEGRAY


#define myGraphDestroyDeclare(FUNCNAME)        \
  static void FUNCNAME(void) ;

#define myGraphDestroy(FUNCNAME, FUNCGRAPH)    \
static void FUNCNAME(void)                     \
{                                              \
  messAssert(graphActive() == FUNCGRAPH) ;     \
                                               \
  graphDestroy() ;                             \
                                               \
  FUNCGRAPH = GRAPH_NULL ;                     \
                                               \
  return ;                                     \
}








static void     externalCommand();

/* Global variables *******************/

static int    **conservCount=0,	   /* Matrix of conservation values  - 21 x maxLen */
              **colorMap=0,	   /* Matrix of conservation colours - 21 x maxLen */
               *conservResidues=0, /* Array of number of residues present in each column */

                colorScheme = COLORSIM,  /* Current colour scheme mode (used for checkmarks) */
                color_by_conserv = 1,    /* 1 => by id or similarity; 0 => by residue  */
                id_blosum=1,	         /* Paint similar residues too */
                color_by_similarity = 1, /* 0 => by id */


                maxLen=0,	/* Columns in alignment */
                nseq=0,		/* Rows in alignment */
                verbose=0,
                debug=0,
                xmosaic = 0,
                ip,		/* Int pointer used in array operations */
                init_sort=0, 
                init_rmPartial=0, 
                rmPickLeftOn,
                pickedCol = 0,	/* Last picked column */
                colorRectangles=1,
                maxNameLen,	/* Max string length of any sequence name */
                maxStartLen=0,	/* Max string length of any sequence start */
                maxEndLen=0,	/* Max string length of any sequence end */
                maxScoreLen=0,	/* Max string length of any sequence score */
                AlignXstart=0,	/* x Start position of alignment in alignment window */
                AlignYstart=0,	/* y Start position of alignment in alignment window */
                Alignwid,	/* Witdh of shown alignment */
                Alignheig,	/* Height of shown alignment */
                statsBox,
                HscrollSliderBox,
                VscrollSliderBox,
                scrLeftBox,
                scrRightBox,
                scrUpBox,
                scrDownBox,
                scrLeftBigBox,
                scrRightBigBox,
                scrUpBigBox,
                scrDownBigBox,
                colorButtonBox,
                sortButtonBox,
                editButtonBox,
                maxBox,
                midBox,
                lowBox,
                gridOn = 0,
                IN_FORMAT = MUL,
                saved = 1,
                displayScores = 0,
                Vcrosshair,
                Highlight_matches,
                colorByResIdOn = 0, /* colour by residue type above identity cutoff */
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
                saveFormat[50],
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
static void saveMul(void);
static void saveFasta(void);
static void saveMSF(void);
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
static void treeBootstrapStats(treeNode *tree);
static void wrapAlignmentWindow(void);
static void mkNonRedundant(double cutoff);
static void mkNonRedundantPrompt(void);
static void rmOutliers(void);
static void rmPartial(void);
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


static void menuCheck(MENU menu, int mode, int thismode, char *str)
{
    if (mode == thismode)
	menuSetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
    else
	menuUnsetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
}
static void setMenuCheckmarks(void)
{
    menuCheck(colorMenu, colorScheme, COLORBYRESIDUE, colorByResidueStr);
    menuCheck(colorMenu, colorScheme, COLORSCHEMESTANDARD, colorSchemeStandardStr);
    menuCheck(colorMenu, colorScheme, COLORSCHEMEGIBSON, colorSchemeGibsonStr);
    menuCheck(colorMenu, colorScheme, COLORSCHEMECYS, colorSchemeCysStr);
    menuCheck(colorMenu, colorScheme, COLORSCHEMEEMPTY, colorSchemeEmptyStr);
    menuCheck(colorMenu, colorScheme, COLORSIM, colorSimStr);
    menuCheck(colorMenu, colorScheme, COLORID, colorIdStr);
    menuCheck(colorMenu, colorScheme, COLORIDSIM, colorIdSimStr);
    menuCheck(colorMenu, 1, colorRectangles, colorRectStr);
    menuCheck(colorMenu, 1, printColorsOn, printColorsStr);
    menuCheck(colorMenu, 1, ignoreGapsOn, ignoreGapsStr); 
    menuCheck(editMenu, 1, bc->rmEmptyColumnsOn, rmEmptyColumnsStr);

    if (!graphExists(showColorGraph))
	graphNewBoxMenu(colorButtonBox, colorMenu);
    else
	graphNewBoxMenu(colorButtonBox, colorEditingMenu);

    graphNewBoxMenu(sortButtonBox, sortMenu);

    graphNewBoxMenu(editButtonBox, editMenu);

    menuCheck(belvuMenu, 1, xmosaic, xmosaicFetchStr);
    graphNewMenu(belvuMenu);
}


static void Help(void)
{
    graphMessage (messprintf("\
Belvu - a multiple alignment viewer.\n\
\n\
Version %s\n\
Copyright (C) Erik Sonnhammer, 1995\n\
\n\
\n\
ALIGNMENT:\n\
  %d sequences, %d residues wide.\n\
\n\
LEFT MOUSE BUTTON:\n\
  Click on residue to show coordinate.\n\
  Double click on sequence to fetch.\n\
  Scrollbars: drag and page.\n\
  Toolbar buttons:\n\
     File:   Display this text.\n\
     Edit:   Remove many sequences.\n\
     Colour: Edit colour scheme window.\n\
     Sort:   Sort alphabetically.\n\
\n\
MIDDLE MOUSE BUTTON:\n\
  In alignment: centre on crosshair.\n\
  Scrollbars: drag and jump.\n\
\n\
RIGHT MOUSE BUTTON:\n\
  Toolbar menus (file menu elsewhere).\n\
\n\
KEYBOARD:\n\
  Arrow keys: Move screenfulls.\n\
  Home: Go to Top.\n\
  End: Go to Bottom.\n\
  Insert: Go to Start of line.\n\
  Delete: Go to End of line.\n\
For more details, see http://www.sanger.ac.uk/~esr/Belvu.html\n\
", belvuVersion, nseq, maxLen));
}


static void fetchAln(ALN *alnp)
{
  if (xmosaic)
    {
      static char *browser = NULL ;

      if (!browser)
	{
	  printf("Looking for WWW browsers ...\n") ;

	  if (!findCommand("netscape", &browser) &&
	      !findCommand("Netscape", &browser) &&
	      !findCommand("Mosaic", &browser) &&
	      !findCommand("mosaic", &browser) &&
	      !findCommand("xmosaic", &browser))
	    {
	      messout("Couldn't find any WWW browser.  Looked for "
		      "netscape, Netscape, Mosaic, xmosaic & mosaic");
	      return ;
	    }
	}

      printf("Using WWW browser %s\n", browser);
      fflush(stdout);
      system(messprintf("%s http://www.sanger.ac.uk/cgi-bin/seq-query?%s&", 
			browser,
			alnp->fetch));
      /* OLD: system(messprintf("xfetch '%s' &", alnp->fetch));*/
    }
  else
    {
      char  *env, fetchProg[1025]="";

      if ((env = getenv("BELVU_FETCH")) )
	strcpy(fetchProg, env);
      else
	strcpy(fetchProg, "efetch");
	    
      externalCommand(messprintf("%s '%s' &", fetchProg, alnp->fetch));
    }

  return ;
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


static void HscrollUp(double x, double y) 
{
    double X = x - SliderOffset;

    AlignXstart = (X - HscrollStart)/(HscrollEnd - HscrollStart) * maxLen;
    unregister();
    belvuRedraw();

    return ;
}

static void HscrollDrag(double x, double y) 
{
    double X = x - SliderOffset;

    if (X < HscrollStart) X = HscrollStart;
    if (X + HsliderLen > HscrollEnd) X = HscrollEnd - HsliderLen;

    graphBoxShift(HscrollSliderBox, X, VscrollEnd+HSCRWID);
    graphRegister(LEFT_UP, HscrollUp);

    return ;
}


static void VscrollUp(double x, double y) 
{
    double Y = y - SliderOffset;

    AlignYstart = (Y - VscrollStart)/(VscrollEnd - VscrollStart) * nseq;
    unregister();
    belvuRedraw();

    return ;
}

static void VscrollDrag(double x, double y) 
{
    double Y = y - SliderOffset;

    if (Y < VscrollStart) Y = VscrollStart;
    if (Y + VsliderLen > VscrollEnd) Y = VscrollEnd - VsliderLen;

    graphBoxShift(VscrollSliderBox, HscrollEnd+VSCRWID, Y);
    graphRegister(LEFT_UP, VscrollUp);

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

static void updateStatus(ALN *alnp, int col)
{
    int i, gaps=0, star;

    if (!alnp)
      return;

    strcpy(stats, " ");

    if (col)
      strcat(stats, messprintf("Column %d: ", col));

    strcat(stats, messprintf("%s/%d-%d", alnp->name, alnp->start, alnp->end));
    
    if (col)
      {
	strcat(stats, messprintf("  %c = ", alnp->seq[col-1]));
	
	for (star = i = 0; i < col; i++)
	  {
	    if (isGap(alnp->seq[i])) gaps++;
	    else if (alnp->seq[i] == '*') star = 1;
	  }
	
	if (star)
	  {
	    strcat(stats, "(unknown position due to insertion)");
	  }
	else
	  {
	    strcat(stats, messprintf("%d", col-1 + alnp->start - gaps));
	  }
      }

    strcat(stats, messprintf(" (%d match", Highlight_matches));
    if (Highlight_matches > 1)
      strcat(stats, "es");
    strcat(stats, ")");

    graphBoxDraw(statsBox, BLACK, LIGHTBLUE);

    return ;
}

static void xorVcrosshair(int mode, double x)
{
    static double xold;

    if (mode == 1)
	graphXorLine(xold, VscrollStart-HSCRWID, xold, VscrollEnd+HSCRWID);

    graphXorLine(x, VscrollStart-HSCRWID, x, VscrollEnd+HSCRWID);

    xold = x;
}

static void middleDrag(double x, double y) 
{
    pickedCol = x2col(&x);

    if (pickedCol == Vcrosshair) return;

    updateStatus(bc->highlightedAln, pickedCol);

    Vcrosshair = pickedCol;

    xorVcrosshair(1, x);
}
static void middleUp(double x, double y) 
{
    pickedCol = x2col(&x);

    AlignXstart = pickedCol - Alignwid/2;
    unregister();
    belvuRedraw();
}

static void leftDown(double x, double y) 
{
    double oldlinew = graphLinewidth(0.05);
    
    graphRectangle(0.5, floor(y), wrapLinelen+0.5, floor(y)+1);
    graphLinewidth(oldlinew);

    graphRedraw();

}


static void middleDown(double x, double y) 
{
    if (x > HscrollStart && x < HscrollEnd &&
	y > VscrollEnd+HSCRWID && y < VscrollEnd+2*HSCRWID)
      {
	graphRegister(MIDDLE_DRAG, HscrollDrag);
	graphRegister(MIDDLE_UP, HscrollUp);
	SliderOffset = HsliderLen/2;
	HscrollDrag(x, y);
      }
    else if (y > VscrollStart && y < VscrollEnd &&
	x > HscrollEnd+VSCRWID && x < HscrollEnd+2*VSCRWID)
      {
	graphRegister(MIDDLE_DRAG, VscrollDrag);
	graphRegister(MIDDLE_UP, VscrollUp);
	SliderOffset = VsliderLen/2;
	VscrollDrag(x, y);
      }
    else if (x > HscrollStart-VSCRWID && x < HscrollEnd+VSCRWID &&
	     y > VscrollStart-HSCRWID && y < VscrollEnd+HSCRWID) 
      {
	/* In alignment */

	graphRegister(MIDDLE_DRAG, middleDrag);
	graphRegister(MIDDLE_UP, middleUp);

	pickedCol = x2col(&x);
	updateStatus(bc->highlightedAln, pickedCol);

	Vcrosshair = pickedCol;

	xorVcrosshair(0, x);
    }
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


/* Highlight the box of bc->highlightedAln
 */
static void highlight(int ON)
{
    int i, box;

    if (!bc->highlightedAln) return;

    box = bc->highlightedAln->nr - AlignYstart;

    if (box > 0 && box <= Alignheig) {

	/* The alignment */
	if (ON) 
	    graphBoxDraw(box, WHITE, BLACK);
	else
	    graphBoxDraw(box, BLACK, WHITE);
    }

    /* All names * /
    for (i = Highlight_matches = 0; i < Alignheig; i++)
	if (!strcmp(arrp(Align, AlignYstart+i, ALN)->name, bc->highlightedAln->name)) {
	    if (ON) {
		graphBoxDraw(Alignheig+i+1, WHITE, BLACK);
		Highlight_matches++;
	    }
	    else
		graphBoxDraw(Alignheig+i+1, BLACK, WHITE);
	}
*/

    /* All names - also count all matches */
    for (i = Highlight_matches = 0; i < nseq; i++)
	if (!strcmp(arrp(Align, i, ALN)->name, bc->highlightedAln->name)) {
	    Highlight_matches++;

	    if (i >= AlignYstart && i < AlignYstart+Alignheig) {
		if (ON) 
		    graphBoxDraw(i-AlignYstart+1+Alignheig, WHITE, BLACK);
		else
		    graphBoxDraw(i-AlignYstart+1+Alignheig, BLACK, WHITE);
	    }
	}
}

/* General purpose routine to convert a string to ALN struct.
   Note: only fields Name, Start, End are filled!
 */
static void str2aln(char *src, ALN *alnp) {

    char *tmp = g_malloc(strlen(src)+1) ;

    strcpy(tmp, src);
    stripCoordTokens(tmp);
    if (sscanf(tmp, "%s%d%d", alnp->name, &alnp->start, &alnp->end) != 3)
      {
	messout("Name to field conversion failed for %s (%s)", src, tmp);
	return;
      }

    if (strlen(alnp->name) > MAXNAMESIZE)
      fatal("buffer overrun in %s !!!!!!!!!\n", "str2aln") ;

    g_free(tmp);
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
	bc->highlightedAln = treenodePicked->aln;
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

/* Click on name -> fetch 
 * Click on sequence -> update status bar 
 *
 * NOTE: the first Alignheig boxes are the alignment ones, then come Alignheig name ones.
 */
static void boxPick (int box, double x, double y, int modifier_unused)
{
  static int lastbox=0;
  static Graph  ga, lastga ;
  static time_t   lasttime;

    pickedCol = (int)x+1+AlignXstart;

    ga = graphActive();

    /* printf("Picked box %d\n", box); */

    if (box > Alignheig*2)
      {
	/* Click outside alignment - check for scrollbar controls */
	if (box == scrLeftBox) scrLeft();
	else if (box == scrRightBox) scrRight();
	else if (box == scrUpBox) scrUp();
	else if (box == scrDownBox) scrDown();
	else if (box == scrLeftBigBox) scrLeftBig();
	else if (box == scrRightBigBox) scrRightBig();
	else if (box == scrUpBigBox) scrUpBig();
	else if (box == scrDownBigBox) scrDownBig();
	else if (box == HscrollSliderBox)
	  {
	    graphRegister(LEFT_DRAG, HscrollDrag);
	    SliderOffset = x;
	  }
	else if (box == VscrollSliderBox)
	  {
	    graphRegister(LEFT_DRAG, VscrollDrag);
	    SliderOffset = y;
	  }

	return;
      }


    /* Alignment clicked */

    /* Turn last box off */
    if (lastga && lastbox) {
	graphActivate(lastga);
    }
    highlight(0);

    /* Turn all selections off */
    if (!box) {
	bc->highlightedAln = 0;
	return;
    }

    aln.nr = (box > Alignheig ? box-Alignheig : box);
    aln.nr += AlignYstart;
    if (!arrayFind(Align, &aln, &ip, (void*)nrorder)) {
	messout("boxPick: Cannot find row %d in alignment array", aln.nr);
	return;
    }
    bc->highlightedAln = alnp = arrp(Align, ip, ALN);
    

    /* Double click */
    if (box == lastbox && (time(0) - lasttime < DBL_CLK_DELAY)) 
      {
	lasttime = time(0) - DBL_CLK_DELAY;  /* To avoid triple clicks */

	/* 'Remove many sequences' mode */
	if (rmPickLeftOn)
	  {
	    rmPicked();
	    return;
	  }
	
	/* if (box > Alignheig) */
	fetchAln(alnp);
      }
    else
      lasttime = time(0);

    /* Single click */

    /* Turn this box on */
    graphActivate(ga);
    highlight(1);
    lastbox = box;
    lastga = ga;

    if (box <= Alignheig)
      {
	/* In alignment - update status bar */
	updateStatus(alnp, pickedCol);
      }
    else
      {
	updateStatus(alnp, 0);
	graphPostBuffer(alnp->name) ;
      }

    return ;
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



static void centerHighlighted(void)
{
    if (bc->highlightedAln) 
	AlignYstart = bc->highlightedAln->nr - Alignheig/2;
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

  /* arrows */
  scrLeftBox = graphBoxStart();
  graphTriangle(HscrollStart-VSCRWID, ny-0.5-0.5*HSCRWID,
		HscrollStart, ny-0.5-HSCRWID,
		HscrollStart, ny-0.5);
  graphBoxEnd();
  graphBoxDraw(scrLeftBox, BLACK, SCRBACKCOLOR);
  scrRightBox = graphBoxStart();
  graphTriangle(HscrollEnd+VSCRWID, ny-0.5-0.5*HSCRWID,
		HscrollEnd, ny-0.5-HSCRWID,
		HscrollEnd, ny-0.5);
  graphBoxEnd();
  graphBoxDraw(scrRightBox, BLACK, SCRBACKCOLOR);
    
  /* big-scroll boxes */
  scrLeftBigBox = graphBoxStart();
  graphRectangle(HscrollStart, ny-0.5-HSCRWID,
		 x1, ny-0.5);
  graphBoxEnd();
  graphBoxDraw(scrLeftBigBox, BLACK, SCRBACKCOLOR);
  scrRightBigBox = graphBoxStart();
  graphRectangle(x2, ny-0.5-HSCRWID,
		 HscrollEnd, ny-0.5);
  graphBoxEnd();
  graphBoxDraw(scrRightBigBox, BLACK, SCRBACKCOLOR);

  /* slider */
  HscrollSliderBox = graphBoxStart();
  graphScrollBar(x1, ny-0.5-HSCRWID, x2, ny-0.5);
  graphBoxEnd();


  /* Vertical scrollbar */

  /* arrows */
  scrUpBox = graphBoxStart();
  graphTriangle(nx-0.5-0.5*VSCRWID, VscrollStart-HSCRWID,
		nx-0.5-VSCRWID, VscrollStart,
		nx-0.5, VscrollStart);
  graphBoxEnd();
  graphBoxDraw(scrUpBox, BLACK, SCRBACKCOLOR);
  scrDownBox = graphBoxStart();
  graphTriangle(nx-0.5-0.5*VSCRWID, VscrollEnd+HSCRWID,
		nx-0.5-VSCRWID, VscrollEnd,
		nx-0.5, VscrollEnd);
  graphBoxEnd();
  graphBoxDraw(scrDownBox, BLACK, SCRBACKCOLOR);

  /* big-scroll boxes */
  scrUpBigBox = graphBoxStart();
  graphRectangle(nx-0.5-VSCRWID, VscrollStart,
		 nx-0.5, y1);
  graphBoxEnd();
  graphBoxDraw(scrUpBigBox, BLACK, SCRBACKCOLOR);
  scrDownBigBox = graphBoxStart();
  graphRectangle(nx-0.5-VSCRWID, y2,
		 nx-0.5, VscrollEnd);
  graphBoxEnd();
  graphBoxDraw(scrDownBigBox, BLACK, SCRBACKCOLOR);

  /* slider */
  VscrollSliderBox = graphBoxStart();
  graphScrollBar(nx-0.5-VSCRWID, y1, nx-0.5, y2);
  graphBoxEnd();

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




/* 
   Format the name/start-end string 

   For convenience, used by writeMul
*/
char *writeMulName(ALN *aln) {

    static char name[MAXNAMESIZE+50];
    char *cp, *namep, GRname[MAXNAMESIZE*2+2], GRfeat[MAXNAMESIZE+1];

    if (aln->markup == GC) 
	return messprintf("#=GC %s", aln->name);

    /* NOTE: GR lines have the feature concatenated in aln->name - must separate */
    if (aln->markup == GR) {
	namep = GRname;
	strncpy(namep, aln->name, 50);
	cp = strchr(namep, ' ');
	strncpy(GRfeat, cp+1, 50);
	*cp = 0;
    }
    else
	namep = aln->name;

    if (!saveCoordsOn) {
	strcpy(name, messprintf("%s", namep));
    }
    else {
	strcpy(name, messprintf("%s%c%d-%d", namep, saveSeparator, aln->start, aln->end));
    }

    if (aln->markup == GR)
	return messprintf("#=GR %s %s", name, GRfeat);
    else
	return name;
}


static void saveAlignment(void)
{
    if (!strcmp(saveFormat, MSFStr))
	saveMSF();
    else if (!strcmp(saveFormat, FastaAlnStr))
	saveFasta();
    else if (!strcmp(saveFormat, FastaStr))
	saveFasta();
    else
	saveMul();
}


static void graphButtonCheck(char* text, void (*func)(void), double x, double *y, int On)
{
    char buttont[1024];

    strcpy(buttont, messprintf("%s %s", On ? "*" : " ", text));

    graphButton(buttont, func, x, *y);
    *y += 2;
}

static void saveMSFSelect(void) {
    strcpy(saveFormat, MSFStr);
    saveRedraw();
}
static void saveMulSelect(void){
    strcpy(saveFormat, MulStr);
    saveRedraw();
}
static void saveFastaAlnSelect(void){
    strcpy(saveFormat, FastaAlnStr);
    saveRedraw();
}
static void saveFastaSelect(void){ /* unaligned */
    strcpy(saveFormat, FastaStr);
    saveRedraw();
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
    setTreeScaleCorr();
    treeSettings();
}
static void treeKIMURAselect(void){
    strcpy(treeDistString, KIMURAstr);
    treeDistCorr = KIMURA;
    setTreeScaleCorr();
    treeSettings();
}
static void treeSTORMSONNselect(void){
    strcpy(treeDistString, STORMSONNstr);
    treeDistCorr = STORMSONN;
    setTreeScaleCorr();
    treeSettings();
}
static void treeScoredistselect(void){
    strcpy(treeDistString, Scorediststr);
    treeDistCorr = Scoredist;
    setTreeScaleCorr();
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





/* Convert plot coordinates to screen coordinates */
static double xTran(double x)
{
    return(5 + x*conservationXScale);
}
static double yTran(double y)
{
    return(3 + conservationRange - y*conservationYScale);
}


static void conservationSetWindow(void)
{
    ACEIN cons_in;

    if (!(cons_in = messPrompt("Window size:", 
			       messprintf("%d", conservationWindow), 
			       "iz", 0)))
	return;
    
    aceInInt(cons_in, &conservationWindow);
    aceInDestroy(cons_in);
    cons_in = NULL ;

    conservationPlot();
}
static void conservationSetLinewidth(void)
{
    ACEIN cons_in;
    
    if (!(cons_in = messPrompt("Line width:", 
		     messprintf("%.2f", conservationLinewidth), 
			       "fz", 0)))
	return;
    
    aceInDouble(cons_in, &conservationLinewidth);
    aceInDestroy(cons_in);
    cons_in = NULL ;

    conservationPlot();
}
static void conservationSetScale(void)
{
    if (!(ace_in = messPrompt("X and Y scales:", 
		     messprintf("%.1f %.1f", conservationXScale, conservationYScale), 
			      "ffz", 0)))
	return;
    
    aceInDouble(ace_in, &conservationXScale);
    aceInDouble(ace_in, &conservationYScale);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    conservationPlot();
}


static void conservationPlot()
{
    int i;
    double 
	oldlinew, mincons, maxcons, avgcons,
	XscaleBase, YscaleBase,
	*smooth, sum, f;

    /* init window */
    if (!graphActivate(conservationGraph))
      {
	conservationGraph = graphCreate (TEXT_FULL_SCROLL, "Belvu conservation profile", 0, 0, 1, 0.4);
	graphTextFormat(FIXED_WIDTH);
      }

    oldlinew = graphLinewidth(conservationLinewidth);
    graphClear();
    graphMenu(conservationMenu);

    /* smooth  conservation by applying a window */
    smooth = g_malloc(maxLen*sizeof(double));
    sum = 0.0;
    for (i = 0; i < conservationWindow; i++) 
	sum += conservation[i];
    smooth[conservationWindow/2] = sum/conservationWindow;
    for (     ; i < maxLen; i++) {
	sum -= conservation[i-conservationWindow];
	sum += conservation[i];
	smooth[i-conservationWindow/2] = sum/conservationWindow;
    }


    /* Find max and min and avg conservation */
    maxcons = -1;
    mincons = 10000;
    avgcons = 0;
    for (i = 0; i < maxLen; i++) {
	if (smooth[i] > maxcons) maxcons = smooth[i];
	if (smooth[i] < mincons) mincons = smooth[i];
	avgcons += conservation[i];
    }
    avgcons /= maxLen*1.0;
    conservationRange = (maxcons - mincons)*conservationYScale;
    graphTextBounds(maxLen*conservationXScale+30, conservationRange*conservationYScale+10);

    graphText(messprintf("Window = %d", conservationWindow), 10, 1);

    /* Draw x scale */
    YscaleBase = -1;
    graphLine(xTran(0), yTran(YscaleBase), xTran(maxLen), yTran(YscaleBase));
    for (i=0; i < maxLen; i += 10) {
	graphLine(xTran(i), yTran(YscaleBase), 
		  xTran(i), yTran(YscaleBase)+0.5);
	if (i+5 < maxLen-1)
	    graphLine(xTran(i+5), yTran(YscaleBase), 
		      xTran(i+5), yTran(YscaleBase)+0.25);
	graphText(messprintf("%d", i), xTran(i)-0.5, yTran(YscaleBase)+1);
    }
    
    /* Draw y scale */
    XscaleBase = -1;
    graphLine(xTran(XscaleBase), yTran(0), xTran(XscaleBase), yTran(maxcons));
    for (f=mincons; f < maxcons; f += 1) {
	graphLine(xTran(XscaleBase), yTran(f), 
		  xTran(XscaleBase-0.5), yTran(f));
	if (f+0.5 < maxcons) graphLine(xTran(XscaleBase), yTran(f+0.5), 
				       xTran(XscaleBase-0.25), yTran(f+0.5));
	graphText(messprintf("%.0f", f), xTran(XscaleBase)-3, yTran(f)-0.5);
    }
    
    /* Draw average line */
    graphColor(RED);
    graphLine(xTran(0), yTran(avgcons), xTran(maxLen), yTran(avgcons));
    graphText("Average conservation", xTran(maxLen), yTran(avgcons));
    graphColor(BLACK);

    /* Plot conservation */
    for (i=1; i <maxLen; i++) {
	graphLine(xTran(i), yTran(smooth[i-1]), 
		  xTran(i+1), yTran(smooth[i]));
    }
    
    graphLinewidth(oldlinew);
    graphRedraw();
    g_free(smooth);

}


static void alphaSort(void)
{
    if (bc->highlightedAln)
	alncpy(&aln, bc->highlightedAln);

    arraySort(Align, (void*)alphaorder);

    if (bc->highlightedAln) {
	if (!alignFind(Align, &aln, &ip)) {
	    messout("Cannot find back highlighted seq after sort - probably a bug");
	    bc->highlightedAln = 0;
	}
	else
	    bc->highlightedAln = arrp(Align, ip, ALN);
    }

    arrayOrder();
    centerHighlighted();
    belvuRedraw();
}


static void organismSort(void)
{
    if (bc->highlightedAln)
	alncpy(&aln, bc->highlightedAln);

    arraySort(Align, (void*)organismorder);

    if (bc->highlightedAln) {
	if (!alignFind(Align, &aln, &ip)) {
	    messout("Cannot find back highlighted seq after sort - probably a bug");
	    bc->highlightedAln = 0;
	}
	else
	    bc->highlightedAln = arrp(Align, ip, ALN);
    }

    arrayOrder();
    centerHighlighted();
    belvuRedraw();
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


static void treeSwapNode(treeNode *node)
{
    void *tmp;

    tmp = node->left;
    node->left = node->right;
    node->right = tmp;
}



/* Save tree in New Hampshire format */
static void treePrintNH_init(void) {
    FILE *file;
    treeStruct *treestruct;
    
    if (!(file = filqueryopen(dirName, fileName, "","w", "Write New Hampshire format tree to file:"))) 
      return;

    if (!graphAssFind(assVoid(1), &treestruct))
      {
	messout("Could not find tree for this graph");
	return;
      }

    treePrintNH(treestruct, treestruct->head, file);
    fprintf(file, ";\n");
    fclose(file);

    return ;
}


static void treePrintNode(treeNode *node) {
    if (node->name) printf("%s ", node->name);
}


static void treeFindOrthologsRecur(treeStruct *treestruct, treeNode *node, char *org1, char *org2) {

    if (!node || !node->left || !node->right) return;

    if (debug) {
	printf("\n 1 (%s, seq=%s):  ", node->left->organism, node->left->name);
	printf("\n 2 (%s, seq=%s)\n: ", node->right->organism, node->right->name);
    }

    /* The open-minded way */
    if (node->left->organism && node->right->organism &&
	node->left->organism != node->right->organism) {
	
	graphBoxDraw(node->box, WHITE, BLACK);
	printf("\nSpecies 1 (%s):  ", node->left->organism);
	treeTraverse(node->left, treePrintNode);
	printf("\nSpecies 2 (%s): ", node->right->organism);
	treeTraverse(node->right, treePrintNode);
	printf("\n");
    }
    
    /* The narrow-minded way * /
    if ((node->left->organism == org1 && node->right->organism == org2) ||
	(node->left->organism == org2 && node->right->organism == org1)) {
	
	printf("Found orthologs");
    }*/

    else {
	treeFindOrthologsRecur(treestruct, node->left, org1, org2);
	treeFindOrthologsRecur(treestruct, node->right, org1, org2);
    }    
}

static void treeFindOrthologs(void) {
    treeStruct *treestruct;
    char *org1 = NULL, *org2 = NULL ;

    
    if (!graphAssFind(assVoid(1), &treestruct))
      {
	messout("Could not find tree for this graph");
	return;
      }

    /* Dialog to define org1 and org2 */

    treeFindOrthologsRecur(treestruct, treestruct->head, org1, org2);
}


static void printMtx(double **mtx) {
    int i, j;

    printf ("\n");
    for (i = 0; i < nseq; i++) {
	for (j = 0; j < nseq; j++)
	    printf("%6.2f ", mtx[i][j]);
	printf ("\n");
    }
}



/* Find back the ALN pointers lost in treeSortBatch
*/
static void treeFindAln(treeNode *node) {

    if (!node->name) return;

    str2aln(node->name, &aln);

    if (!alignFind(Align, &aln, &ip)) {
	messout("Cannot find back seq after sort - probably a bug");
    }
    else
	node->aln = arrp(Align, ip, ALN);
}

static void treeSort(void)
{
    ALN aln;

    if (bc->highlightedAln)
	alncpy(&aln, bc->highlightedAln);
    
    treeSortBatch();

    /* After treeSortBatch the treeNodes point to the wrong ALN,
       because arraySort only copies contents and not entire
       structs. This can be fixed by either a linear string search or
       by remaking the tree... */
    treeTraverse(treeHead, treeFindAln);

    if (bc->highlightedAln) {
	if (!alignFind(Align, &aln, &ip)) {
	    messout("Cannot find back highlighted seq after sort - probably a bug");
	    bc->highlightedAln = 0;
	}
	else
	    bc->highlightedAln = arrp(Align, ip, ALN);
	centerHighlighted();
    }
    
    treeDraw(treeHead);

    belvuRedraw();
}


static void simSort(void) {
    highlightScoreSort('P');
}

static void idSort(void) {
    highlightScoreSort('I');
}



static void scoreSort(void)
{
    if (bc->highlightedAln)
	alncpy(&aln, bc->highlightedAln);

    arraySort(Align, (void*)scoreorder);

    if (bc->highlightedAln) {
	if (!alignFind(Align, &aln, &ip)) {
	    messout("Cannot find back highlighted seq after sort - probably a bug");
	    bc->highlightedAln = 0;
	}
	else
	    bc->highlightedAln = arrp(Align, ip, ALN);
    }

    arrayOrder();
    AlignYstart = 0;
    belvuRedraw();
}


static void printScore(void)
{
    if (!bc->highlightedAln) {
	messout("Select a line first");
	return;
    }

    printf("%.1f %s/%d-%d\n", 
	   bc->highlightedAln->score,
	   bc->highlightedAln->name,
	   bc->highlightedAln->start,
	   bc->highlightedAln->end);
    fflush(stdout);
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
	pw, ph,	 /* pixel width and height */
	output_probs = 0,
	init_tree = 0,
	only_tree = 0,
        show_ann = 0;
    double   
        makeNRinit = 0.0,
        init_rmEmptyColumns = 0.0,
        init_rmGappySeqs = 0.0,
	cw;      /* character width of initial alignment */

/*    extern int printOnePage;*/

    int          optc;
    extern int   optind;
    extern char *optarg;
    char        *optstring="aBb:CcGgil:L:m:n:O:o:PpQ:q:RrS:s:T:t:uX:z:";

    static char *cc_date = 
#if defined(__DATE__)
    __DATE__
#else
    ""
#endif
    ;

    char *usage;
    static char usageText[] = "\
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
 -X <#>      Print UPGMA-based subfamilies at cutoff #.\n\
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
 -b <#>      Apply boostrap analysis with # bootstrap samples\n\
 -B          Print out bootstrap trees and exit\n\
               (Negative value -> display bootstrap trees on screen)\n\
 -O <label>  Read organism info after this label (default OS)\n\
 -t <title>  Set window title.\n\
 -g          Draw grid line (for debugging).\n\
 -u          Start up with uncoloured alignment (faster).\n\
\n\
 Some X options:\n\
 -acefont <font>   Main font.\n\
 -font    <font>   Menu font.\n\
\n\
 Note: X options only work after \"setenv POSIXLY_CORRECT\"\n\
\n\
 setenv BELVU_FETCH to desired sequence fetching program.\n\
\n\
 For documentation, see:\n\
 http://sonnhammer.sbc.su.se/Belvu.html\n\
\n\
 by Erik.Sonnhammer@sbc.su.se\n\
 Version ";


    usage = g_malloc(strlen(usageText) + strlen(belvuVersion) + strlen(cc_date) + 20);
    sprintf(usage, "%s%s, compiled %s\n", usageText, belvuVersion, cc_date);

    /* Set up tree defaults */
    strcpy(treeMethodString, NJstr);
    treeMethod = NJ;
    strcpy(treeDistString, Scorediststr);
    treeDistCorr = Scoredist;
    strcpy(treePickString, SWAPstr);
    treePickMode = NODESWAP;
    setTreeScaleCorr();

    while ((optc = getopt(argc, argv, optstring)) != -1)
	switch (optc) 
	{
	case 'a': show_ann = 1;                  break;
	case 'B':  outputBootstrapTrees = 1; break;
	case 'b': treebootstraps = atoi(optarg); break;
	case 'C': saveCoordsOn = 0;              break;
	case 'c': verbose = 1;                   break;
	case 'G': penalize_gaps = 1;             break;
	case 'g': gridOn = 1;                    break;
	case 'l': 
	    colorCodesFile = g_malloc(strlen(optarg)+1);
	    strcpy(colorCodesFile, optarg);      break;
	case 'i': 
	    ignoreGapsOn = 1;                           break;
	case 'L': 
	    markupColorCodesFile = g_malloc(strlen(optarg)+1);
	    strcpy(markupColorCodesFile, optarg);break;
	case 'm': 
	    readMatchFile = g_malloc(strlen(optarg)+1);
	    strcpy(readMatchFile, optarg); 
	                                         break;
	case 'n':
	    makeNRinit = atof(optarg);           break;
	case 'O': 
	    strncpy(OrganismLabel, optarg, 2);    break;
	case 'o': 
	    output_format = g_malloc(strlen(optarg)+1);
	    strcpy(output_format, optarg);       break;
	case 'P': init_rmPartial = 1;            break;
	case 'Q': init_rmEmptyColumns = atof(optarg); break;
	case 'q': init_rmGappySeqs = atof(optarg); break;
	case 'p': output_probs = 1;              break;
	case 'R': stripCoordTokensOn = saveCoordsOn = 0; break;
	case 'r': IN_FORMAT = RAW;               break;
	case 'S': 
	    switch (*optarg)
	    {
	    case 'a': init_sort = SORT_ALPHA;    break;
	    case 'o': init_sort = SORT_ORGANISM; break;
	    case 's': init_sort = SORT_SCORE;    break;
	    case 'S': init_sort = SORT_SIM;      break;
	    case 'i': init_sort = SORT_ID;       break;
	    case 'n': 
		treeMethod = NJ;
		init_sort = SORT_TREE;           break;
	    case 'u': 
		treeMethod = UPGMA; 
		init_sort = SORT_TREE;           break;
	    default : fatal("Illegal sorting order: %s", optarg);
	    }                                    break;
	case 's': 
	    scoreFile = g_malloc(strlen(optarg)+1);
	    strcpy(scoreFile, optarg);           break;
	case 'T': 
	    for (optargc = optarg; *optargc; optargc++) {
		switch (*optargc)
		{
		case 'n': 
		    strcpy(treeMethodString, NJstr);
		    treeMethod = NJ;          break;
		case 'u': 
		    strcpy(treeMethodString, UPGMAstr);
		    treeMethod = UPGMA;       break;
		case 'c':
		    treeColorsOn = 0;         break;
		case 'd':
		    treeShowBranchlen = 1;    break;
		case 'I':
		    only_tree=1;
		case 'i':
		    init_tree = 1;            break;
		case 'j':
		    strcpy(treeDistString, JUKESCANTORstr);
		    treeDistCorr = JUKESCANTOR;
		    setTreeScaleCorr();         break;
		case 'k':
		    strcpy(treeDistString, KIMURAstr);
		    treeDistCorr = KIMURA;
		    setTreeScaleCorr();         break;
		case 'o':
		    treeCoordsOn = 0;         break;
		case 'p':
		    treePrintDistances = 1;  
		    init_tree = 1;            break;
		case 'R':
		    treeReadDistancesON = 1;  break;
		case 's':
		    strcpy(treeDistString, STORMSONNstr);
		    treeDistCorr = STORMSONN;
		    setTreeScaleCorr();         break;
		case 'b':
		    strcpy(treeDistString, SCOREDISTstr);
		    treeDistCorr = SCOREDIST;
		    setTreeScaleCorr();         break;
		case 'r':
		    strcpy(treeDistString, UNCORRstr);
		    treeDistCorr = UNCORR;
		    treeScale = 1.0;          break;
		default : fatal("Illegal sorting order: %s", optarg);
		}
	    }                                 break;
	case 't': 
	    strncpy(Title, optarg, 255);         break;
	case 'u': colorRectangles = 0;           break;
	case 'X': mksubfamilies_cutoff = atof(optarg);  break;
	case 'z': saveSeparator = *optarg;       break;
	default : fatal("Illegal option");
	}

    if (argc-optind < 1) {
	fprintf(stderr, "%s\n", usage); 
	exit(1);
    }

    if (!strcmp(argv[optind], "-")) {
	pipe = stdin;
	if (!*Title) strcpy(Title, "stdin");
    }
    else {
	if (!(pipe = fopen(argv[optind], "r")))
	    fatal("Cannot open file %s", argv[optind]);
	if (!*Title) strncpy(Title, argv[optind], 255);
	strncpy(fileName, argv[optind], FIL_BUFFER_SIZE);
    }

/*    printOnePage = TRUE;*/

 #if defined(LINUX) 
    Align = arrayCreate(100000, ALN);
    organismArr = arrayCreate(10000, ALN);
    /* Note: this excessive initialization size is a temporary fix to stop
    Belvu from segfaulting under linux.  The segfaulting happens during 
    arrayInsert expansions but only under linux. As I can't find any bugs
    in Belvu it seems the array package may be bugged under linux. For
    now, avoiding arrayInsert expansions by initializing it to over twice the largest
    family in Pfam 040401 buys me some time to find the problem. */    
#else
    Align = arrayCreate(100, ALN);
    organismArr = arrayCreate(100, ALN);
#endif

    if (treeReadDistancesON) 
      {
        /* Should this really be either or?  Problem: cannot read organism info when reading tree */
	treeReadDistancesPipe = pipe;
	treeReadDistancesNames();
	
	init_tree = 1;
	only_tree = 1;
	treeCoordsOn = 0;
      }
    else
      {
        readMul(pipe);
      }

    if (!arrayMax(organismArr))
      suffix2organism();

    setOrganismColors();

    if (scoreFile) readScores(scoreFile);
    
    init_sort_do();

    if (!matchFooter && readMatchFile) {
	if (!(file = fopen(readMatchFile, "r"))) 
	    fatal("Cannot open file %s", readMatchFile);
	readMatch(file);
	fclose(file);
    }
    else if (matchFooter) 
    {	 
        readMatch(pipe);
	fclose(pipe);
    }

    if (!treeReadDistancesON) {
	checkAlignment();
	setConservColors();
    }

    if (verbose)
      {
	/* Print conservation statistics */
	int i, j, max, consensus = 0 ;
	double totcons = 0.0;

	printf("\nColumn Consensus        Identity       Conservation\n");
	printf  ("------ ---------  -------------------  ------------\n");

	for (i = 0; i < maxLen; i++)
	  {
	    max = 0;

	    for (j = 1; j < 21; j++)
	      {
		if (conservCount[j][i] > max)
		  {
		    max = conservCount[j][i];
		    consensus = j;
		  }
	      }

	    printf("%4d       %c      %4d/%-4d = %5.1f %%  %4.1f\n", 
		   i+1, b2a[consensus], max, nseq, (double)max/nseq*100, conservation[i]);
	    totcons += conservation[i];
	  }

	printf ("\nAverage conservation = %.1f\n", totcons/(maxLen*1.0));

	exit(0);
      }

    initResidueColors();
    if (colorCodesFile) {
	if (!(file = fopen(colorCodesFile, "r"))) 
	    fatal("Cannot open file %s", colorCodesFile);
	readColorCodes(file, color);
 	colorScheme = COLORBYRESIDUE;
	bc->color_by_conserv = bc->colorByResIdOn = 0;
   }
    initMarkupColors();
    if (markupColorCodesFile) {
	if (!(file = fopen(markupColorCodesFile, "r"))) 
	    fatal("Cannot open file %s", markupColorCodesFile);
	readColorCodes(file, markupColor);
    }
    
    if (makeNRinit)
	mkNonRedundant(makeNRinit);
    
    if (init_rmPartial)
	rmPartial();

    if (init_rmEmptyColumns)
	rmEmptyColumns(init_rmEmptyColumns/100.0);

    if (init_rmGappySeqs) {
	rmGappySeqs(init_rmGappySeqs);
	rmFinaliseGapRemoval();
    }

    if (output_format) {
	if (!strcasecmp(output_format, "Stockholm") ||
	    !strcasecmp(output_format, "Mul") ||
	    !strcasecmp(output_format, "Selex"))
	    writeMul(stdout);
	else if (!strcasecmp(output_format, "MSF"))
	    writeMSF(stdout);
	else if (!strcasecmp(output_format, "FastaAlign")) {
	    strcpy(saveFormat, FastaAlnStr);
	    writeFasta(stdout);
	}
	else if (!strcasecmp(output_format, "Fasta")) {
	    strcpy(saveFormat, FastaStr);
	    writeFasta(stdout);
	}
	else if (!strcasecmp(output_format, "tree")) {
	    treeStruct *treestruct = g_malloc(sizeof(treeStruct));
	    separateMarkupLines();
	    treestruct->head = treeMake(1);
	    treePrintNH(treestruct, treestruct->head, stdout);
	    printf(";\n");
	}
	else 
	    fatal("Illegal output format: %s", output_format);
	exit(0);
    }
    
    if (outputBootstrapTrees && treebootstraps > 0) {
	treeBootstrap();
	exit(0);
    } 
   
    if (output_probs) {
	outputProbs(stdout);
	exit(0);
    }

    if (bc->mksubfamilies_cutoff) {
	mksubfamilies(bc->mksubfamilies_cutoff);
	exit(0);
    
    }
    fprintf(stderr, "\n%d sequences, max len = %d\n", nseq, maxLen);

    /* Try to get 8x13 font for menus, if not set on command line */
    for ( i=0; i < argc; i++)
	if (!strcmp(argv[i], "-font")) break;
    if (i == argc) {
	argvAdd(&argc, &argv, "-font");
	argvAdd(&argc, &argv, "8x13");
    }

    graphInit(&argc, argv);
    gexInit(&argc, argv);

    graphScreenSize(&screenWidth, &screenHeight, &fontwidth, &fontheight, &pw, &ph);
    Aspect = (pw/fontwidth)/(ph/fontheight);
    VSCRWID = HSCRWID/Aspect;

    if (show_ann) showAnnotation();

    /* Calculate screen width of alignment */
    cw = 1 + maxNameLen+1;
    if (maxStartLen) cw += maxStartLen+1;
    if (maxEndLen)   cw += maxEndLen+1;
    if (maxScoreLen) cw += maxScoreLen+1;
    cw += maxLen + ceil(VSCRWID) + 2;

     colorMenu = menuInitialise ("color", (MENUSPEC*)colorMENU);
    colorEditingMenu = menuInitialise ("color", (MENUSPEC*)colorEditingMENU);
    sortMenu = menuInitialise ("sort", (MENUSPEC*)sortMENU);
    editMenu = menuInitialise ("edit", (MENUSPEC*)editMENU);
    showColorMenu =  menuInitialise ("", (MENUSPEC*)showColorMENU);
    saveMenu =  menuInitialise ("", (MENUSPEC*)saveMENU);
    belvuMenu = menuInitialise ("belvu", (MENUSPEC*)mainMenu);
    treeGUIMenu = menuInitialise ("", (MENUSPEC*)treeGUIMENU);
    treeDistMenu = menuInitialise ("", (MENUSPEC*)treeDistMENU);
    treePickMenu = menuInitialise ("", (MENUSPEC*)treePickMENU);
    if (!displayScores) {
	menuSetFlags(menuItem(sortMenu, "Sort by score"), MENUFLAG_DISABLED);
	menuSetFlags(menuItem(editMenu, "Remove sequences below given score"), MENUFLAG_DISABLED);
	menuSetFlags(menuItem(belvuMenu, "Print score and coords of line"), MENUFLAG_DISABLED);
	menuSetFlags(menuItem(belvuMenu, "Output score and coords of line"), MENUFLAG_DISABLED);
    }
    menuSetFlags(menuItem(colorMenu, thresholdStr), MENUFLAG_DISABLED);

   if (init_tree)
      {
	treeDisplay();

	if (only_tree)
	  {

	    /*ACEOUT out = aceOutCreateToFile("t", "w", 0);
	    graphGIF(treeGraph, out, 0);*/

	    graphLoop(FALSE);
	    graphFinish();
	    
	    exit(0);
	  }
      }

    if (outputBootstrapTrees && treebootstraps < 0)
      {	/* Display [treebootstraps] bootstrap trees */
	bc->treebootstrapsDisplay = TRUE;
	
	treebootstraps = -treebootstraps;

	treeBootstrap();

	/* graphLoop(FALSE);
	graphFinish();
	exit(0); */
      }



    /* Create the main belvu graph display of aligned sequences. */

    belvuGraph = graphCreate (TEXT_FIT, messprintf("Belvu:  %s", Title), 0, 0, 
			      cw*screenWidth/fontwidth, 0.44*screenHeight);


    graphRegister(PICK, boxPick);
    graphRegister(MIDDLE_DOWN, middleDown);
    graphRegister(RESIZE, belvuRedraw);
    graphRegister(KEYBOARD, keyboard);
    graphRegister(DESTROY, belvuDestroy) ;

    if (!colorCodesFile) colorSim() ;

    belvuRedraw() ;

    graphLoop(FALSE) ;

    graphFinish() ;

    return(0) ;
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
    
    setConservColors();

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

static void settingsPick(int box, double x_unused, double y_unused, int modifier_unused)
{
    if (box == maxBox) graphTextScrollEntry(maxText,0,0,0,0,0);
    else if (box == midBox) graphTextScrollEntry(midText,0,0,0,0,0);
    else if (box == lowBox) graphTextScrollEntry(lowText,0,0,0,0,0);
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

	    if (bc->id_blosum) paintMessage(&i);

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

	setConservColors();
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

	if (bc->colorByResIdOn) setConservColors();
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
      setConservColors();

      belvuRedraw();
    }

  return ;
}


static void colorRect(void)
{
    colorRectangles = (!colorRectangles);
    belvuRedraw();
}


static void readColorCodesMenu(void)
{
    FILE *fil;

    if (!(fil = filqueryopen(dirName, fileName, "","r", "Read file:"))) return;

    readColorCodes(fil, color);
    belvuRedraw();
}


static void saveColorCodes(void)
{
    FILE *fil;
    int i;

    if (!(fil = filqueryopen(dirName, fileName, "","w", "Save to file:"))) return;

    for (i = 1; i < 21; i++)
      {
	fprintf(fil, "%c %s\n", b2a[i], colorNames[color[(unsigned char)(b2a[i])]]);
      }

    fclose(fil);
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


static void colorSchemeCys(void)
{
    colorScheme = COLORSCHEMECYS;

    initResidueColors();
    color['C'] = color['c'] = CYAN;
    color['G'] = color['g'] = RED;
    color['P'] = color['p'] = GREEN;

    colorRes();
}


static void colorSchemeEmpty(void)
{
    colorScheme = COLORSCHEMEEMPTY;

    initResidueColors();
    colorRes();
}


static void colorSchemeStandard(void)
{
    colorScheme = COLORSCHEMESTANDARD;

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

    colorRes();
}


static void colorSchemeGibson(void)
{
    colorScheme = COLORSCHEMEGIBSON;

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
		    if (bc->color_by_conserv)
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


static void rmFinaliseColumnRemoval(void)
{
    rmGappySeqs(100.0);
    rmFinalise();
}



/* Get rid of seqs that are less than x% identical with any of the others. 
 */
static void rmOutliers(void)
{
    int i,j, n=0;
    ALN *alni, *alnj;
    static double cutoff = 20.0, id, maxid;

    if (!(ace_in = messPrompt ("Remove sequences less identical to others than (%) :", 
		      messprintf("%.0f", cutoff), 
			       "fz", 0)))
	return;

    aceInDouble(ace_in, &cutoff);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    for (i = 0; i < nseq-1; ) {
	alni = &g_array_index(bc->alignArr, ALN, i);

	for (maxid=0, j = 0; j < nseq; j++) {
	    if (i == j) continue;
	    alnj = arrp(Align, j, ALN);
	    id = identity(alni->seq, alnj->seq);
	    if (id > maxid) maxid = id;
	}

	if (maxid < cutoff) {

	    printf("%s/%d-%d was max %.1f%% identical to any other sequence and was removed.\n",
		   alni->name, alni->start, alni->end, maxid);
	    
	    /* Remove entry */
	    n++;
	    nseq--;
	    if (bc->highlightedAln == alni) bc->highlightedAln = 0;
	    arrayRemove(Align, alni, (void*)nrorder);
	    saved = 0;
	}
	else i++;
    }

    printf("%d sequences removed at the %.0f%% level.  %d seqs left\n\n", n, cutoff, nseq);

    arrayOrder();
    rmFinaliseGapRemoval();
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

static void readLabels(void)
{
    int row, col, seqpos, seqlen;
    ALN *alnrow;
    char 
      *labelseq,		/* The raw sequence of labels */
      *label,			/* The mapped sequence of labels, 1-maxlen */
      *cp, *cq;
    FILE *fil;

    if (!bc->highlightedAln) {
        messout("pick a sequence first");
	return;
    }

    labelseq = g_malloc(maxLen+1);
    label = g_malloc(maxLen+1);

    if (!(fil = filqueryopen(dirName, fileName, "","r", 
			     messprintf("Read labels of %s from file:", bc->highlightedAln->name))))
      return;

    /* read file */
    cq = labelseq;
    while (!feof (fil)) { 
	if (!fgets (line, MAXLENGTH, fil)) break;
	cp = line;
	while (*cp) {
	    if (isalpha(*cp)) *cq++ = *cp;
	    cp++;
	}
    }
    fclose(fil);
    
    /* Warn if seq too long, return if too short */
    seqlen = bc->highlightedAln->end - bc->highlightedAln->start +1;
    if (strlen(labelseq) > seqlen)
        messout(messprintf("The sequence of labels is longer (%d) than the sequence (%d).  "
			   "Hope that's ok", 
			   strlen(labelseq), seqlen));
    if (strlen(labelseq) < seqlen) {
        messout(messprintf("The sequence of labels is shorter (%d) than the sequence (%d).  "
			   "Aborting", 
			   strlen(labelseq), seqlen));
	return;
    }

    /* map labels to alignment */
    for (col = 0, seqpos = 0; col < maxLen; col++) {
        label[col] = labelseq[seqpos];
	if (isalpha(bc->highlightedAln->seq[col]) && labelseq[seqpos+1]) seqpos++;
    }

    for (row = 0; row < nseq; row++) {
	alnrow = arrp(Align, row, ALN);
	for (col = 0; col < maxLen; col++)
  	    if (isalpha(alnrow->seq[col]))
	        alnrow->seq[col] = label[col];
    }
      
    g_free(labelseq);
    belvuRedraw();
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



static void mkNonRedundantPrompt(void)
{
    static double cutoff = 80.0;

    if (!(ace_in = messPrompt ("Remove sequences more identical than (%):", 
		      messprintf("%.0f", cutoff), 
			       "fz", 0)))
	return;
    
    aceInDouble(ace_in, &cutoff);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    mkNonRedundant(cutoff);
}


static void rmScore(void)
{
    int 
	i;
    static double 
	cutoff = 20.0;

    if (!(ace_in = messPrompt ("Remove sequences scoring less than: ", 
		      messprintf("%.0f", cutoff), 
			       "fz", 0)))
	return;

    aceInDouble(ace_in, &cutoff);
    aceInDestroy(ace_in);
    ace_in = NULL ;

    scoreSort();

    /* Save bc->highlightedAln */
    if (bc->highlightedAln)
	alncpy(&aln, bc->highlightedAln);

    for (i = 0; i < nseq; ) {
	alnp = &g_array_index(bc->alignArr, ALN, i);
	if (alnp->score < cutoff) {
	    fprintf(stderr, "Removing %s/%d-%d (score %.1f\n",
		   alnp->name, alnp->start, alnp->end, alnp->score);
		
	    nseq--;
	    if (bc->highlightedAln == alnp) bc->highlightedAln = 0;
	    arrayRemove(Align, alnp, (void*)nrorder);
	    saved = 0;
	}
	else
	    i++;
    }
    
    arrayOrder();

    /* Find bc->highlightedAln in new array */
    if (bc->highlightedAln) {
	if (!arrayFind(Align, &aln, &ip, (void*)scoreorder)) {
	    bc->highlightedAln = 0;
	}
	else
	    bc->highlightedAln = arrp(Align, ip, ALN);
	centerHighlighted();
    }

    AlignYstart = 0;

    rmFinaliseGapRemoval();
}


static void listIdentity(void)
{
    int 
      i,j,n ;
    ALN 
	*alni, *alnj;
    double 
	sc, totsc=0, maxsc=0, minsc=1000000,
	id, totid=0.0, maxid=0.0, minid=100.0;

    for (i = n = 0; i < nseq-1; i++) {
	alni = &g_array_index(bc->alignArr, ALN, i);
	for (j = i+1; j < nseq; j++, n++) {
	    alnj = arrp(Align, j, ALN);
	    
	    id = identity(alni->seq, alnj->seq);
	    totid += id;
	    if (id > maxid) maxid = id;
	    if (id < minid) minid = id;
	    
	    sc = score(alni->seq, alnj->seq);
	    totsc += sc;
	    if (sc > maxsc) maxsc = sc;
	    if (sc < minsc) minsc = sc;
	    
	    printf("%s/%d-%d and %s/%d-%d are %.1f%% identical, score=%f\n",
		   alni->name, alni->start, alni->end,
		   alnj->name, alnj->start, alnj->end,
		   id, sc);
		
	}
	printf("\n");
    }
    printf("Maximum %%id was: %.1f\n", maxid);
    printf("Minimum %%id was: %.1f\n", minid);
    printf("Mean    %%id was: %.1f\n", totid/n);
    printf("Maximum score was: %.1f\n", maxsc);
    printf("Minimum score was: %.1f\n", minsc);
    printf("Mean    score was: %.1f\n", (double)totsc/n);
}


static void markup (void)
{
    if (!bc->highlightedAln) {
	messout ("Pick a sequence first!");
	return;
    }

    /* Store orignal state in the nocolor field:
       1 = normal sequence
       2 = markup line

       This is needed to restore markup lines (nocolor lines are always markups)
    */

    if (!bc->highlightedAln->nocolor) {
	if (bc->highlightedAln->markup)
	    bc->highlightedAln->nocolor = 2;
	else
	    bc->highlightedAln->nocolor = 1;
	bc->highlightedAln->markup = 1;
    }
    else {
	if (bc->highlightedAln->nocolor == 1) 
	    bc->highlightedAln->markup = 0;
	bc->highlightedAln->nocolor = 0;
    }
    
    belvuRedraw();
}


static void hide (void)
{
    if (!bc->highlightedAln) {
	messout ("Pick a sequence first!");
	return;
    }

    bc->highlightedAln->hide = 1;
    belvuRedraw();
}


static void unhide (void)
{
    int i;

    for (i = 0; i < nseq; i++)
	arrp(Align, i, ALN)->hide = 0;
    
    belvuRedraw();
}


static void rmPicked(void)
{
    if (!bc->highlightedAln) {
	messout ("Pick a sequence first!");
	return;
    }

    printf("Removed %s/%d-%d.  ", bc->highlightedAln->name, bc->highlightedAln->start, bc->highlightedAln->end);

    nseq--;
    arrayRemove(Align, bc->highlightedAln, (void*)nrorder);
    saved = 0;
    arrayOrder();
    bc->highlightedAln = 0;

    printf("%d seqs left\n\n", nseq);

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


static void rmColumnCutoff(void)
{
    int 
	i, j, max, removed=0, oldmaxLen=maxLen;
    static double 
	from = -1,
	to = 0.9,
      cons ;

    if (!bc->color_by_conserv) {
	messout("You must use a conservation coloring scheme");
	return;
    }

    if (!(ace_in = messPrompt ("Remove columns with a (maximum) conservation between x1 and x2 ( i.e. if (c > x1 && c <= x2) ).  Provide x1 and x2:", 
		      messprintf("%.2f %.2f", from, to), 
			       "ffz", 0)))
	return;

    aceInDouble(ace_in, &from);
    aceInDouble(ace_in, &to);
    aceInDestroy(ace_in);
    ace_in = NULL ;


    for (i = maxLen-1; i >= 0; i--) {

	if (color_by_similarity) {
	    cons = conservation[i];
	}
	else {
	    max = 0;
	    for (j = 1; j < 21; j++) {
		if (conservCount[j][i] > max) {
		    max = conservCount[j][i];
		}
	    }	
	    cons = (double)max/nseq;
	}


	if (cons > from && cons <= to) {
	    printf("removing %d, cons= %.2f\n", i+1, cons);
	    rmColumn(i+1, i+1);
	    if (++removed == oldmaxLen) {
		messExit("You have removed all columns.  Prepare to exit Belvu");
	    }
	}
    }

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
	from=1, to=pickedCol;
    
    if (!pickedCol) {
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
	from=pickedCol, to=maxLen;
    
    if (!pickedCol) {
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
    setConservColors();
    belvuRedraw();
}


static void ignoreGaps (void)
{
    ignoreGapsOn = (ignoreGapsOn ? 0 : 1);
    setConservColors();
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




/* Global variables */
static double treeScale;

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


/* Residue colors */
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
static gboolean            arrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer));



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




/* Reset an ALN struct to default values */
void resetALN(ALN *alnp)
{
  *alnp->name = *alnp->fetch = 0;
  alnp->start = alnp->end = alnp->len = alnp->nr = 0;
  alnp->seq = 0;
  alnp->score = 0.0;
  alnp->color = WHITE;
  alnp->markup = 0;
  alnp->organism = 0;
}


void setTreeScaleCorr(const int treeMethod) 
{
  if (treeMethod == UPGMA)
      treeScale = 1.0;
  else if (treeMethod == NJ)
      treeScale = 0.3;
}


void setTreeScale(const double newScale) 
{
  treeScale = newScale;
}


static void readFastaAlnFinalise(BelvuContext *bc, ALN *aln, char *seq)
{
  aln->len = strlen(seq);
  aln->seq = g_strdup(seq);
  
  if (bc->maxLen) 
    {
      if (aln->len != bc->maxLen) 
        g_error("Differing sequence lengths: %d %d", bc->maxLen, aln->len);
    }
  else
    {
      bc->maxLen = aln->len;
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


/* Read in fasta sequences from a file and create a sequence in the given
 * alignments array for each of the fasta sequences. */
static void readFastaAln(BelvuContext *bc, FILE *pipe)
{
  char line[MAXLENGTH+1];
  GString *currentSeq = NULL;
  ALN currentAln;

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

  strcpy(bc->saveFormat, FastaAlnStr);
  
  return ;
}


static void alncpy(ALN *dest, ALN *src)
{
  strncpy(dest->name, src->name, MAXNAMESIZE);
  dest->start = src->start;
  dest->end = src->end;
  dest->len = src->len;
  /* dest->seq = g_malloc(strlen(src->seq));
     strcpy(dest->seq, src->seq); */
  dest->seq = src->seq;
  dest->nr = src->nr;			
  strncpy(dest->fetch, src->fetch, MAXNAMESIZE+10);
  dest->score = src->score;
  dest->color = src->color;
}



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
  
  resetALN(aln);
  
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


/* Finds Entry s from Array  a
 * sorted in ascending order of order()
 * If found, returns TRUE and sets *ip
 * if not, returns FALSE and sets *ip one step left
 */
static gboolean arrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer))
{
  int ord;
  int i = 0 , j = a->len, k;

  if (!j || (ord = orderFunc(s, &g_array_index(a, BootstrapGroup, 0))) < 0)
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

  if ((ord = orderFunc(s, &g_array_index(a, BootstrapGroup, --j))) > 0)
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

      if ((ord = orderFunc(s, &g_array_index(a, BootstrapGroup, k))) == 0)
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
    g_array_index(alignArr, ALN, i).nr = i + 1;
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



/* Calculate percent identity of two strings */
static double identity(char *s1, char *s2, const gboolean penalize_gaps)
{
    int n, id;

    for (n = id = 0; *s1 && *s2; s1++, s2++) {
	if (isGap(*s1) && isGap(*s2)) continue;
	if (isGap(*s1) || isGap(*s2))
	    if (!penalize_gaps) continue;
	n++;
	if (toupper(*s1) == toupper(*s2)) id++;
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

    for (;*s1 && *s2; s1++, s2++) {

	if (isGap(*s1) && isGap(*s2)) 
	    continue;
	else if (isGap(*s1) || isGap(*s2)) {
	    if (penalize_gaps) sc -= 0.6;
	}
	else
	  sc += (double) BLOSUM62[a2b[(unsigned char)(*s1)]-1][a2b[(unsigned char)(*s2)]-1];
    }
    
    return sc;
}


/* Separate markuplines to another array before resorting
*/
void separateMarkupLines(BelvuContext *bc)
{
  bc->markupAlignArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN), 100);
    
  ALN aln;
  if (bc->highlightedAln)
    alncpy(&aln, bc->highlightedAln);

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

  if (bc->highlightedAln) 
    {
      int idx = 0;
      
      if (!alignFind(bc->alignArr, &aln, &idx))
        bc->highlightedAln = 0;
      else
        bc->highlightedAln = &g_array_index(bc->alignArr, ALN, idx);
    }
}


/* Reinsert markuplines after mother sequence or if orphan at bottom */
static void reInsertMarkupLines(BelvuContext *bc)
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


/*
* Expect format:
*
* name1 name2 name3
* 1-1   1-2   1-3
* 2-1   2-2   2-3
* 3-1   3-2   3-3
*/
static void treeReadDistances(BelvuContext *bc, double **pairmtx)
{
  char   *p;
  int    i = 0, j = 0 ;
  double  d;
  char line[MAXLENGTH+1];
  
  while (!feof (bc->treeReadDistancesPipe))
    { 
      if (!fgets (line, MAXLENGTH, bc->treeReadDistancesPipe))
	break;
      
      j = 0;
      
      while ( (p = strtok(j ? 0 : line, " \t")) && !isspace(*p))
	{
	  d = atof(p);
	  DEBUG_OUT("%d  %d  %f\n", i, j, d);
	  pairmtx[i][j] = d;
	  j++;
	}

      if (j != bc->alignArr->len)
	g_error("nseq = %d, but read %d distances in row %d", bc->alignArr->len, j, i);

      i++;
    }

  if (j != bc->alignArr->len)
    g_error("nseq = %d, but read %d distance rows", bc->alignArr->len, i);

  return ;
}


/* Correct an observed distance to an estimated 'true' distance.

   Adapted from Lasse Arvestad's lapd program

   od = observed distance in percent 

*/
static double treeJUKESCANTOR(double od)
{
    double cd;

    od /= 100;

    od = 20.0/19.0*od;
    if  (od > 0.95) od = 0.95; /* Limit to 300 PAM */

    cd = -19.0/20.0 * log(1-od);

    return cd*100;
}



/* Correct an observed distance to an estimated 'true' distance.

   Adapted from Lasse Arvestad's lapd program

   od = observed distance in percent 

*/
static double treeKimura(double od)
{
    double adjusted, cd;

    od /= 100;

    adjusted = od + 0.2 * od*od;

    if (adjusted > 0.95) adjusted = 0.95;   /* Limit to 300 PAM */

    cd = -log(1 - adjusted);

    return cd*100;
}



/* Correct an observed distance to an estimated 'true' distance.

   Based on Christian Storm's 5th order polynomial curve fitting of ROSE simulated data

   od = observed distance in percent 

*/
static double treeSTORMSONN(double od)
{
    double cd;			/* Corrected distance */

    if (od > 91.6) return 1000.0;	/* Otherwise log of negative value below */

    od /= 100;

    cd= -log(1 
             -0.95844*od 
             -0.69957*od*od
             +2.4955*od*od*od
             -4.6353*od*od*od*od
             +2.8076*od*od*od*od*od);

    /* printf(" od=%f  cd=%f\n", od, cd); */
   
    if (cd > 3.0) cd=3.0; /* Limit to 300 PAM */
    
    return cd*100;
}


static double treeSCOREDIST(char *seq1, char *seq2, BelvuContext *bc)
{
    int len,
	sc = 0, 
	s1sc = 0, 
	s2sc = 0;
    double od, cd, 
	maxsc, 
	expect; 
    char *s1, *s2;
    
#define mtx BLOSUM62
    

    s1 = seq1;
    s2 = seq2;


    /* Calc scores */
    for (len=0; *s1; s1++, s2++) 
      {
        if (bc->penalize_gaps) 
          {
	    if (isGap(*s1) || isGap(*s2)) sc -= 0.6;	
          }

	if (!isGap(*s1) && !isGap(*s2)) 
          {
            sc += mtx[a2b[(unsigned char)(*s1)]-1][a2b[(unsigned char)(*s2)]-1];
            s1sc += mtx[a2b[(unsigned char)(*s1)]-1][a2b[(unsigned char)(*s1)]-1];
            s2sc += mtx[a2b[(unsigned char)(*s2)]-1][a2b[(unsigned char)(*s2)]-1];
	    len++;
          }
      }
    
    maxsc = (s1sc + s2sc) / 2.0;
    
    /* Calc expected score */
    expect =  -0.5209 * len;

    od = ((double)sc - expect) / (maxsc - expect);

    if (!len || od < 0.05) od = 0.05;  /* Limit to 300 PAM;  len==0 if no overlap */
    if (od > 1.0) od = 1.0; 
    
    cd = -log(od);

    /* printf ("len=%d  sc=%.2f  maxsc=%.2f  expect=%.2f  maxsc-expect=%.2f  od=%.3f\n", len, sc, maxsc, expect, maxsc-expect, od); */

    cd = cd * 100.0 * 1.337 /* Magic scaling factor optimized for Dayhoff data */ ;

    if (cd > 300) cd = 300;  /* Limit to 300 PAM */
    
/*    if (!bc->penalize_gaps) {
	g_free(s1);
	g_free(s2);
    } */

    return cd;
}


static void treeTraverse(BelvuContext *bc, TreeNode *node, void (*func)(BelvuContext *bc, TreeNode *treeNode)) 
{
  if (!node) 
    return;

  treeTraverse(bc, node->left, func);
  func(bc, node);
  treeTraverse(bc, node->right, func);
}



/* Sum branchlengths, allow going up parents

   If going up parents, arg1 = parent  ;  arg2 = child.

   If you're not interested in going up parents, simply call me with
   the same node in both arguments.
   
 */
static double treeSize3way(TreeNode *node, TreeNode *fromnode) 
{
  int root = 0;
  TreeNode *left, *right;
  double len;
  
  if (!node) 
    return 0.0;

  /* Get the correct branch length */
  if (node->left == fromnode || node->right == fromnode) /* Going up the tree */
    len = fromnode->branchlen;
  else
    len = node->branchlen;

  if (node->left == fromnode) 
    {
      left = node->parent;
      if (!left) 
        root = 1;
    }
  else
    {
      left = node->left;
    }

  if (node->right == fromnode) 
    {
      right = node->parent;
      if (!right) 
        root = 1;
    }
  else 
    {
      right = node->right;
    }

  if (root)
    {
      double retval;

      /* go across tree root */
      if (left) /* Coming from right */
        retval = treeSize3way(left, left);
      else  /* Coming from left */
        retval = treeSize3way(right, right);

      DEBUG_OUT("Returning (%.1f + %.1f) = %.1f\n", fromnode->branchlen, retval, fromnode->branchlen + retval);

      return fromnode->branchlen + retval;
    }
  else 
    {
      double 
        l = treeSize3way(left, node),
        r = treeSize3way(right, node),
        retval = (l + r)/2.0 + len;

      DEBUG_OUT("Returning (%.1f + %.1f)/2 + %.1f = %.1f\n", l, r, len, retval);

      return retval;
    }
}



/* Calculate the difference between left and right trees if the tree
   were to be rerooted at this node.

   What is calculated (bal) is the difference between 'left' and
   'right' subtrees minus the length of the branch itself. If this
   difference is negative, a perfectly balanced tree can be made.  For
   imperfectly balanced trees we want to root at the branch that gives
   the best balance.  However, perfectly balanced trees are all
   'perfect' so here we chose the branch with most equal subtrees.

   Actually it is not "left" and "right" but "down" and "up" subtrees.  */
static void treeCalcBalance(BelvuContext *bc, TreeNode *node) 
{
  double bal, lweight, rweight;

  if (node == bc->treeBestBalancedNode) 
    return;

  DEBUG_OUT("Left/Downstream weight\n");

  lweight = treeSize3way(node, node);

  DEBUG_OUT("Right/Upstream weight\n");

  rweight = treeSize3way(node->parent, node);

  bal = fabsf(lweight - rweight) - node->branchlen;
  
  DEBUG_OUT("Node=%s (branchlen = %.1f).  Weights = %.1f  %.1f. Bal = %.1f\n", 
            node->name, node->branchlen, lweight, rweight, bal);

  if (bal < bc->treeBestBalance) 
    { /* better balance */
      if (bc->treeBestBalance > 0.0 ||
          /* If previous tree was not perfectly balanced, or
             If previous tree was perfectly balanced - choose root with best subtree balance */
          fabsf(lweight - rweight) < bc->treeBestBalance_subtrees)
        {
          DEBUG_OUT("            %s has better balance %.1f < %.1f\n", node->name, bal, bc->treeBestBalance);
          
          bc->treeBestBalancedNode = node;
          bc->treeBestBalance = bal;
          bc->treeBestBalance_subtrees = fabsf(lweight - rweight);
        }
    }
}


static TreeNode *treeParent2leaf(TreeNode *newparent, TreeNode *curr)
{
  if (!curr->parent) 
    { /* i.e. the old root */
      if (curr->left == newparent)
        {
          newparent->branchlen += curr->right->branchlen; /* Add second part of this vertex */
          return curr->right;
        }
      else
        {
          newparent->branchlen += curr->left->branchlen;
          return curr->left;
        }
    } 
  else
    {
      if (curr->left == newparent) 
        {
          /* Change the link to the new parent to the old parent */
          curr->left = treeParent2leaf(curr, curr->parent);
          curr->left->branchlen = curr->branchlen;
        }
      else
        {
          curr->right = treeParent2leaf(curr, curr->parent);
          curr->right->branchlen = curr->branchlen;
        }
    }

  return curr;
}


/* Balance the two sides of a tree branch by the weights of the two subtrees.
   The branch has two sides: llen and rlen.

   Rationale: The correct center point balances the treesizes on both sides.
   
   Method: Calculate half the unbalance, add it to the appropriate side.

   No real theory for this, but it seems to work in easy cases 
*/
static void treeBalanceByWeight(TreeNode *lnode, TreeNode *rnode, double *llen, double *rlen)
{
    double adhocRatio = 0.95;

    double 
	halfbal, 
	branchlen = *rlen+*llen;
    double lweight = treeSize3way(lnode, lnode) /*- lnode->branchlen*/;
    double rweight = treeSize3way(rnode, rnode) /*- rnode->branchlen*/;

    halfbal = fabsf(lweight-rweight) / 2.0;

    if (halfbal < *llen && halfbal < *rlen) {

	if (lweight < rweight) {
	    *llen += halfbal;
	    *rlen -= halfbal;
	}
	else {
	    *llen -= halfbal;
	    *rlen += halfbal;
	}
    }
    else {
	/* The difference is larger than the final branch -
	   give nearly all weight to the shorter one (cosmetic hack) */
	if (lweight < rweight) {
	    *llen = branchlen*adhocRatio;
	    *rlen = branchlen*(1.0 - adhocRatio);
	}
	else {
	    *rlen = branchlen*adhocRatio;
	    *llen = branchlen*(1.0 - adhocRatio);
	}
    }
}



/* Rerooting works roughly this way:

   - A new node is created with one child being the node chosen as new root 
     and the other child the previous parent of it.  The sum of left and right
     branch lengths of the new root should equal the branch length of the chosen node.
   - All previous parent nodes are visited and converted so that:
         1. The previous parent node becomes a child.
	 2. The new parent is the node that calls.
	 3. The branchlength of the previous parent node is assigned to its parent.
   - When the previous root is reached, it is deleted and the other child
     becomes the child of the calling node.

     Note that treeReroot destroys the old tree, since it reuses the nodes.  
     Use treecpy if you still need it for something later.
*/
static TreeNode *treeReroot(BelvuContext *bc, TreeNode *node)
{
  TreeNode *newroot = g_malloc(sizeof(TreeNode));

  newroot->left = node;
  newroot->right = treeParent2leaf(node, node->parent);

  fillParents(newroot, newroot->left);
  fillParents(newroot, newroot->right);

  newroot->left->branchlen = newroot->right->branchlen = node->branchlen / 2.0;
  treeBalanceByWeight(newroot->left, newroot->right, 
                      &newroot->left->branchlen, &newroot->right->branchlen);

  bc->treeHead = newroot;
  return newroot;
}


/* Find the node which has most equal balance, return new tree with this as root.
 */
static TreeNode *treeFindBalance(BelvuContext *bc, TreeNode *tree) 
{
  double lweight = treeSize3way(tree->left, tree->left);
  double rweight = treeSize3way(tree->right, tree->right);

  bc->treeBestBalancedNode = tree;
  bc->treeBestBalance = fabsf(lweight - rweight);

  bc->treeBestBalance_subtrees = 
    fabsf((lweight - tree->left->branchlen) - (rweight - tree->right->branchlen));
  
  DEBUG_OUT("Initial weights = %.1f  %.1f. Bal = %.1f\n", lweight, rweight, bc->treeBestBalance);

  treeTraverse(bc, tree, treeCalcBalance);
  
  if (bc->treeBestBalancedNode == tree)
    return tree;
  else
    return treeReroot(bc, bc->treeBestBalancedNode);
}


static void fillParents(TreeNode *parent, TreeNode *node)
{
  if (!node) 
    return;

  node->parent = parent;
  fillParents(node, node->left);
  fillParents(node, node->right);
}


static char *fillOrganism(TreeNode *node)
{
  char 
    *leftorganism,
    *rightorganism;

  if (node->name) 
    return node->organism;

  leftorganism = fillOrganism(node->left);
  rightorganism = fillOrganism(node->right);

  /* printf("\nFill: (left=%s, right=%s):  ", leftorganism, rightorganism); */

  node->organism = (leftorganism == rightorganism ? leftorganism : 0);

  return node->organism;
}


static int strcmp_(gconstpointer xIn, gconstpointer yIn)
{
  const char *x = (const char*)xIn;
  const char *y = (const char*)yIn;
  
  int retval = strcmp(x, y);
  return retval;
}


static int BSorder(gconstpointer xIn, gconstpointer yIn)
{
  const BootstrapGroup *x = (const BootstrapGroup*)xIn;
  const BootstrapGroup *y = (const BootstrapGroup*)yIn;

  return strcmp(x->s, y->s);
}



/* Combines left and right sequence groups and insert to bootstrap group list */
static GArray* fillBootstrapGroups(BelvuContext *bc, TreeNode *node, TreeNode *rootnode, int maintree) 
{
  int i;
  int ssize=0;
  BootstrapGroup *BS = (BootstrapGroup *)g_malloc(sizeof(BootstrapGroup));

  if (!node->name)
    {
      /* Internal node */
      /* Combine left node sequences into right node array */

      GArray *left = fillBootstrapGroups(bc, node->left, rootnode, maintree);
      GArray *right = fillBootstrapGroups(bc, node->right, rootnode, maintree);

      /* No need to do root */
      if (node == rootnode)
        return NULL ;

      /* Combine left and right groups */
      for (i = 0 ; i < left->len; ++i) 
        {
          char *s = g_array_index(left, char*, i);
          g_array_append_val(right, s);
          g_array_sort(right, strcmp_);
        }

      g_array_unref(left);

      /* Create string with group members */
      for (i = 0 ; i < right->len ; ++i) 
        ssize += (strlen(g_array_index(right, char*, i)) + 1);
      
      BS->s = g_malloc(ssize+1);
      
      for (i = 0 ; i < right->len ; ++i) 
        {
          strcat(BS->s, g_array_index(right, char*, i));
          strcat(BS->s, " ");
        }

      /* printf("   New bootstrap group: %s\n", BS->s); */
      
      if (maintree) 
        {
          /* Associate this node with the group string */
          BS->node =  node;
          
          /* Add group string to array of bootstrap groups */
          g_array_append_val(bc->bootstrapGroups, BS);
          g_array_sort(bc->bootstrapGroups, BSorder);
        }
      else
        {
          /* Find group string and increment counter if exists */
          BootstrapGroup *BS2;
          
          int ip = 0;
          if (arrayFind(bc->bootstrapGroups, BS, &ip, (void *)BSorder)) 
            {
              BS2 = &g_array_index(bc->bootstrapGroups, BootstrapGroup, ip);
              BS2->node->boot++;
              /* printf("Found bootgroup %s\n", BS->s); */
            }
          else
            {
              /* printf("Did not find bootgroup %s\n", BS->s); */
            }
          
          g_free(BS->s);
          g_free(BS);
        }
      
      return right;
    }
  else
    {
      /* Leaf node - return as embryonic array with one entry */
      GArray *leaf = g_array_sized_new(FALSE, FALSE, sizeof(char*), bc->alignArr->len);
      g_array_append_val(leaf, node->name);
      g_array_sort(leaf, strcmp_);
      return leaf;
    }

}


static void normaliseBootstraps(BelvuContext *bc, TreeNode *node) 
{
  node->boot = 100*node->boot/(double)bc->treebootstraps;
}


/* Bootstrap the internal nodes in a tree.
   1. Set up a list (array) of all internal nodes with leaf content and pointer to node
   2. In bootstrap tree, for each internal node, check if its contents exists in list. If so, increment node's bootstrap count
   3. Turn increments to percentages
*/
static void treeBootstrapStats(BelvuContext *bc, TreeNode *tree)
{
  /* Traverse tree, fill array bootstrapGroups */
  bc->bootstrapGroups = g_array_sized_new(FALSE, FALSE, sizeof(BootstrapGroup), bc->alignArr->len);
  fillBootstrapGroups(bc, tree, tree, 1);
  
  treeBootstrap(bc);
    
  /* Turn increments to percentages */
  treeTraverse(bc, tree, normaliseBootstraps);
}


static GArray *copyAlignArray(GArray *inputArr)
{
  GArray *result = g_array_sized_new(FALSE, FALSE, sizeof(ALN), inputArr->len);
  memcpy(result->data, inputArr->data, inputArr->len * sizeof(ALN));

  int i = 0;
  for ( ; i < inputArr->len; ++i)
    {
      g_array_index(result, ALN, i).seq = g_strdup(g_array_index(inputArr, ALN, i).seq);
    }

  return result;
}


static void columnCopy(GArray *alignArrDest, int destIdx, GArray *alignArrSrc, int srcIdx)
{
    int i;

    for (i = 0; i < alignArrSrc->len; ++i)
      {
        g_array_index(alignArrDest, ALN, i).seq[destIdx] = g_array_index(alignArrSrc, ALN, i).seq[srcIdx];
      }
}


/* Draws clickable boxes for horizontal lines.  The routine must be in sync with
   treeDrawNode, but can not be integrated since a box may
   accidentally overwrite some text.  Therefore all boxes must be 
   drawn before any text or lines.
*/
static double treeDrawNodeBox(BelvuContext *bc, Tree *tree, TreeNode *node, double x) 
{
  double y, yl, yr;
  int box;

  if (!node) 
    return 0.0;

  yl = treeDrawNodeBox(bc, tree, node->left, x + node->branchlen*treeScale);
  yr = treeDrawNodeBox(bc, tree, node->right, x + node->branchlen*treeScale);

  if (yl) 
    {
      y = (yl + yr) / 2.0;
    }
  else 
    {
      y = bc->tree_y++;
    }

  /* Make box around horizontal lines */
//  box = graphBoxStart();
//  graphLine(x + node->branchlen*treeScale, y, x, y);
//  oldlinew = graphLinewidth(0.0); graphColor(BG);
//  graphRectangle(x + node->branchlen*treeScale, y-0.5, x, y+0.5);
//  graphColor(BLACK); graphLinewidth(oldlinew);
//  graphBoxEnd();
    
//    graphAssociate(assVoid(100+box), node);

  node->box = box;
  tree->lastNodeBox = box;

  return y;
}


/* The actual tree drawing routine.
   Note: must be in sync with treeDrawNodeBox, which draws clickable
   boxes first.
*/
static double treeDrawNode(BelvuContext *bc, Tree *tree, TreeNode *node, double x) 
{
  double y, yl, yr;
  int leftcolor;

  if (!node) 
    return 0.0;

  yl = treeDrawNode(bc, tree, node->left, x + node->branchlen*treeScale);
  //  leftcolor = graphColor(BLACK);
  yr = treeDrawNode(bc, tree, node->right, x + node->branchlen*treeScale);

  if (yl) 
    {
      /* internal node */
      y = (yl + yr) / 2.0;

      /* connect children */
//      graphLine(x + node->branchlen*treeScale, yr, x + node->branchlen*treeScale, y);
//      graphColor(leftcolor);
//      graphLine(x + node->branchlen*treeScale, yl, x + node->branchlen*treeScale, y);

      if (node->left->organism != node->right->organism)
        {
          //graphColor(BLACK);
        }
    }
  else 
    {
      /* Sequence name */
      int box ;

      y = bc->tree_y++;

      if (bc->treeColorsOn && node->organism) 
        {
          ALN aln;
          aln.organism = node->organism;

          int ip = 0;
          if (arrayFind(bc->organismArr, &aln, &ip, (void*)organism_order)) 
            {
              //graphColor(arrp(organismArr, ip, ALN)->color);
	    }
	}
      else
        {
          //graphColor(BLACK);
        }

      /* Make clickable box for sequence */
      //      box = graphBoxStart();
//      graphText(node->name, x + node->branchlen*treeScale + 1, y-0.5);
//      graphBoxEnd();
//      graphAssociate(assVoid(100+box), node);

      if (bc->highlightedAln && node->aln == bc->highlightedAln) 
        {
          //graphBoxDraw(box, WHITE, BLACK);
          tree->currentPickedBox = box;
	}
      else if (node->aln) 
        {
          //graphBoxDraw(box, BLACK, node->aln->color);
	}	    

      if (bc->treeShowOrganism && node->organism) 
        {
          //graphText(node->organism, x + node->branchlen*treeScale + 2 + strlen(node->name), y-0.5);
        }
      
      {
        int pos = x + node->branchlen*treeScale + strlen(node->name);
        if (pos > bc->maxTreeWidth) 
          bc->maxTreeWidth = pos;
      }
    }

  /* Horizontal branches */
  //  graphLine(x + node->branchlen*treeScale, y, x, y);

  if (bc->treeShowBranchlen && node->branchlen) 
    {
      char *tmpStr = blxprintf("%.1f", node->branchlen);
      double pos = x + node->branchlen * treeScale * 0.5 - strlen(tmpStr) * 0.5;
      //      graphText(tmpStr, pos, y);
      g_free(tmpStr);
    }

  if (bc->treebootstraps && !node->name && node != bc->treeHead  && !bc->treebootstrapsDisplay) 
    {
      //graphColor(BLUE);
      char *tmpStr = blxprintf("%.0f", node->boot);
      double pos = x + node->branchlen*treeScale - strlen(tmpStr) - 0.5;

      if (pos < 0.0) 
        pos = 0;
      
      printf("%f  %f   \n", node->boot, pos);
      //graphText(tmpStr, pos, y);
      //graphColor(BLACK);
      g_free(tmpStr);
    }

  return y;
}


static void treeDraw(BelvuContext *bc, TreeNode *treeHead)
{
  int i;
  double oldlinew;

  Tree *treeStruct = g_malloc(sizeof(Tree));
  treeStruct->head = treeHead;

  char *title = blxprintf("Belvu %s using %s distances of %s", 
                          bc->treeMethod == NJ ? "Neighbor-joining tree" : "UPGMA tree",
                          bc->treeDistString,
                          bc->Title);

  g_free(title);
  
//  treeGraph = graphCreate (TEXT_FULL_SCROLL, 
//                           title,
//                           0, 0,
//                           (bc->treeMethod == UPGMA ? 130 : 110) / fontwidth * screenWidth, 
//                           (bc->alignArr->len + 7) / fontheight * screenHeight);
//    
//  graphTextFormat(FIXED_WIDTH);
//  graphRegister(PICK, treeboxPick);
//  graphRegister(DESTROY, treeDestroy);

//  graphAssociate(assVoid(1), treestruct);

//  oldlinew = graphLinewidth(treeLinewidth);
//  graphClear();
//  graphMenu(treeMenu);

  bc->maxTreeWidth = 0;
  bc->tree_y = 1;
  treeDrawNodeBox(bc, treeStruct, treeStruct->head, 1.0);
  bc->tree_y = 1;
  treeDrawNode(bc, treeStruct, treeStruct->head, 1.0);
  //graphColor(BLACK);
  //graphTextBounds(bc->maxTreeWidth+2 + (treeShowOrganism ? 25 : 0), nseq+6);

  /* Draw scale */
  if (bc->treeMethod == UPGMA) 
    {
      bc->tree_y += 2;
      //graphLine(1*treeScale, tree_y, 101*treeScale, tree_y);
      for (i=1; i<=101; i += 10) 
        {
          //graphLine(i*treeScale, tree_y, i*treeScale, tree_y+0.5);
          
          if (i < 101) 
            {
              //graphLine((i+5)*treeScale, tree_y, (i+5)*treeScale, tree_y+0.5);
            }
          
          //graphText(messprintf("%d", i-1), (i-0.5)*treeScale, tree_y+1);
        }
    }
  else
    {
      bc->tree_y += 3;
      //graphText("0.1", 5*treeScale, tree_y-1);
      //graphLine(1*treeScale, tree_y, 11*treeScale, tree_y);
      //graphLine(1*treeScale, tree_y-0.5, 1*treeScale, tree_y+0.5);
      //graphLine(11*treeScale, tree_y-0.5, 11*treeScale, tree_y+0.5);
    }
  //graphLinewidth(oldlinew);

  if (bc->treeMethod == NJ) 
    {	
      double lweight, rweight;
      TreeNode *tree = treeStruct->head;

      lweight = treeSize3way(tree->left, tree);
      rweight = treeSize3way(tree->right, tree);

//      graphText((debug ? messprintf("Tree balance = %.1f (%.1f-%.1f)", 
//                                    fabsf(lweight - rweight), lweight, rweight) :
//                 messprintf("Tree balance = %.1f", fabsf(lweight - rweight))),
//                14, tree_y-0.5);
    }

    //    graphRedraw();
}


/* Print tree in New Hampshire format */
void treePrintNH(Tree *tree, TreeNode *node, FILE *file)
 {
   if (!node) 
     return;
    
   if (node->left && node->right) 
     {
       fprintf(file, "(\n");
       treePrintNH(tree, node->left, file);
       fprintf(file, ",\n");
       treePrintNH(tree, node->right, file);
       fprintf(file, ")\n");
       
       if (node != tree->head)	/* Not exactly sure why this is necessary, but njplot crashes otherwise */
         fprintf(file, "%.0f", node->boot+0.5);
     }
   else
     {
       fprintf(file, "%s", node->name);
     }
    
   if (node != tree->head)	/* Not exactly sure why this is necessary, but njplot crashes otherwise */
     fprintf(file, ":%.3f", node->branchlen/100.0);
}


void treeBootstrap(BelvuContext *bc)
{
    Tree *treeStruct = g_malloc(sizeof(Tree));

    separateMarkupLines(bc);
    GArray *alignArrTmp = copyAlignArray(bc->alignArr);

#if defined(__CYGWIN__) || defined(DARWIN)
    srand(time(0));
#else
    srand48(time(0));
#endif

    int iter = 0;
    for (iter = 0; iter < bc->treebootstraps; ++iter) 
      {
	/* Make bootstrapped Alignment */
        int i = 0;
        for (i = 0; i < bc->maxLen; ++i) 
          {
            int src = 0;
            
#if defined(__CYGWIN__) || defined(DARWIN)
	    src = (int)((((double)rand())/RAND_MAX) * bc->maxLen);
#else
	    src = (int)(drand48() * bc->maxLen);
#endif

	    columnCopy(bc->alignArr, i, alignArrTmp, src);
          }


	treeStruct->head = treeMake(bc, 0);

	if (bc->outputBootstrapTrees) 
          {
	    if (bc->treebootstrapsDisplay) 
              {
                treeDraw(bc, treeStruct->head);
              }
	    else
              {
                printf("\n");
                treePrintNH(treeStruct, treeStruct->head, stdout);
                printf(";\n");
              }
          }
	else
          {
	    /* Collect bootstrap statistics */
	    fillBootstrapGroups(bc, treeStruct->head, treeStruct->head, 0);
          }
      }
 
    /* Restore alignment */
    int i = 0;
    for (i = 0; i < bc->alignArr->len; ++i)
      {
	strcpy(g_array_index(bc->alignArr, ALN, i).seq, g_array_index(alignArrTmp, ALN, i).seq);
      }
}



TreeNode *treeMake(BelvuContext *bc, const gboolean doBootstrap)
{
  TreeNode *newnode = NULL ;
  int maxi, maxj;
  double maxid = 0.0, pmaxid, **pairmtx, **Dmtx, **mtx, *src, *trg, 
    *avgdist,		/* vector r in Durbin et al */
    llen, rlen;
  TreeNode **node ;					    /* Array of (primary) nodes.  Value=0 => stale column */
  static BlxHandle handle = 0;
  BlxHandle treeHandle = NULL ;

  if (doBootstrap)
      treeHandle = handle;
    
  /* Allocate memory */
  if (treeHandle)
    handleDestroy(&treeHandle);

  treeHandle = handleCreate();
  pairmtx = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double *));
  Dmtx = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double *));
  node = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(TreeNode *));
  avgdist = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double));

  int i = 0;
  for (i = 0; i < bc->alignArr->len; ++i)
    {
      pairmtx[i] = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double));
      Dmtx[i] = handleAlloc(&treeHandle, bc->alignArr->len*sizeof(double));
      node[i] = handleAlloc(&treeHandle, sizeof(TreeNode));
      node[i]->name =  handleAlloc(&treeHandle, strlen(g_array_index(bc->alignArr, ALN, i).name)+50);

      ALN *aln_i = &g_array_index(bc->alignArr,ALN, i);

      if (!bc->treeCoordsOn) 
	{
	  sprintf(node[i]->name, "%s", aln_i->name);
	}
      else
	{
	  sprintf(node[i]->name, "%s/%d-%d", 
		  aln_i->name,
		  aln_i->start,
		  aln_i->end);
	}
      node[i]->aln = aln_i;
      node[i]->organism = aln_i->organism;
    }
    
	
  if (bc->treeReadDistancesOn)
    {
      treeReadDistances(bc, pairmtx);
    }
  else
    {
      /* Calculate pairwise distance matrix */
      arrayOrder(bc->alignArr);

      for (i = 0; i < bc->alignArr->len - 1; ++i)
        {
          ALN *aln_i = &g_array_index(bc->alignArr, ALN, i);
          
          int j = i+1;
          for (j = i+1; j < bc->alignArr->len; ++j)
            {
              ALN *aln_j = &g_array_index(bc->alignArr, ALN, j);
              pairmtx[i][j] = 100.0 - identity(aln_i->seq, aln_j->seq, bc->penalize_gaps);
              
              if (bc->treeDistCorr == KIMURA) 
                pairmtx[i][j] = treeKimura(pairmtx[i][j]);
              else if (bc->treeDistCorr == JUKESCANTOR) 
                pairmtx[i][j] = treeJUKESCANTOR(pairmtx[i][j]);
              else if (bc->treeDistCorr == STORMSONN) 
                pairmtx[i][j] = treeSTORMSONN(pairmtx[i][j]);
              else if (bc->treeDistCorr == SCOREDIST) 
                pairmtx[i][j] = treeSCOREDIST(aln_i->seq, aln_j->seq, bc);
            }
        }
    }

  if (bc->treePrintDistances) 
    {
      double dist;

      for (i = 0; i < bc->alignArr->len; i++) 
        {
          if (!bc->treeCoordsOn) 
            {
              printf("%s\t", g_array_index(bc->alignArr, ALN, i).name);
            }
          else
            {
              printf("%s/%d-%d\t",
                     g_array_index(bc->alignArr, ALN, i).name,
                     g_array_index(bc->alignArr, ALN, i).start,
                     g_array_index(bc->alignArr, ALN, i).end);
            }
        }
      printf ("\n");
      
      for (i = 0; i < bc->alignArr->len; ++i) 
        {
          int j = 0;
          for (j = 0; j < bc->alignArr->len; ++j) 
            {
              if (i == j) 
                dist = 0.0;
              else if (i < j)
                dist = pairmtx[i][j];
              else
                dist = pairmtx[j][i];
              
              printf("%7.3f\t", dist);
            }
          printf ("\n");
        }
      exit(0);
    }
  

  /* Construct tree */
  int iter = 0;
  for (iter = 0; iter < bc->alignArr->len - 1; ++iter)
    {
      if (bc->treeMethod == NJ)
        {	
          int count = 0;
          
          /* Calculate vector r (avgdist) */
          for (i = 0; i < bc->alignArr->len; ++i)
            {
              if (!node[i]) 
                continue;
              
              avgdist[i] = 0.0;
              
              int j = 0;
              for (count=0, j = 0; j < bc->alignArr->len; j++)
                {
                  if (!node[j]) 
                    continue;
                  
                  avgdist[i] += (i < j ? pairmtx[i][j] : pairmtx[j][i]);
                  count++;
                }
              
              if (count == 2)	/* Hack, to accommodate last node */
                avgdist[i] = 1;
              else
                avgdist[i] /= 1.0*(count - 2);
            }

          /* Calculate corrected matrix Dij */
          if (1 /* gjm */)
            {
              for (i = 0; i < bc->alignArr->len - 1; ++i)
                {
                  if (!node[i])
                    continue;
                  
                  int j = i+1;
                  for (j = i+1; j < bc->alignArr->len; ++j)
                    {
                      if (!node[j]) 
                        continue;
                      
                      Dmtx[i][j] = pairmtx[i][j] - (avgdist[i] + avgdist[j]);
                    }
                }
            }
          else
            {		/* Felsenstein */
              double Q = 0.0;
              
              for (i = 0; i < bc->alignArr->len - 1; ++i)
                {
                  if (!node[i])
                    continue;
                  
                  int j = i+1;
                  for (j = i+1; j < bc->alignArr->len; j++)
                    {
                      if (!node[j]) 
                        continue;
                      
                      Q += pairmtx[i][j];
                    }
                }
              
              for (i = 0; i < bc->alignArr->len - 1; ++i)
                {
                  if (!node[i])
                    continue;
                  
                  int j = i+1;
                  for (j = i+1; j < bc->alignArr->len; j++)
                    {
                      if (!node[j])
                        continue;
                      
                      Dmtx[i][j] = (pairmtx[i][j] +
                                    (2.0*Q)/(count-(count == 2 ? 1 : 2)) - 
                                    avgdist[i] - avgdist[j]) / 2.0;
                    }
                }
            }
          mtx = Dmtx;
        }
      else 
        {
          mtx = pairmtx;
        }
      
#ifdef DEBUG
        {
          printf("Node status, Avg distance:\n");
          for (i = 0; i < bc->alignArr->len; i++)
            printf("%6d ", (node[i] ? 1 : 0)); 
          printf("\n");
          for (i = 0; i < bc->alignArr->len; i++)
            printf("%6.2f ", avgdist[i]); 
          printf("\n\nPairdistances, corrected pairdistances:");
          printMtx(pairmtx);
          printMtx(Dmtx);
          printf("\n");
        }
#endif

      /* Find smallest distance pair in pairmtx */
      maxi = -1;
      maxj = -1;
      maxid = pmaxid = 1000000;
      
      for (i = 0; i < bc->alignArr->len - 1; ++i) 
        {
          if (!node[i]) 
            continue;
          
          int j = i+1;
          for (j = i+1; j < bc->alignArr->len; ++j) 
            {
              if (!node[j]) 
                continue;
              
              /* printf("iter %d, i=%d. j=%d, dist= %f\n", iter, i, j, mtx[i][j]);*/
              
              if (mtx[i][j] < maxid) 
                {
                  maxid = mtx[i][j];
                  pmaxid = pairmtx[i][j];
                  maxi = i;
                  maxj = j;
                }
              else if (bc->treeMethod == NJ && Dmtx[i][j] == maxid && pairmtx[i][j] < pmaxid) 
                {
                  /* To resolve ties - important for tree look! */
                  maxi = i;
                  maxj = j;
                }
            }
        }

      maxid = pairmtx[maxi][maxj]; /* Don't want to point to Dmtx in NJ */

      /* Merge rows & columns of maxi and maxj into maxi
         Recalculate distances to other nodes */
      for (i = 0; i < bc->alignArr->len; ++i)
        {
        if (!node[i]) 
          continue;

        if (i < maxi) 
          trg = &pairmtx[i][maxi];
        else if (i > maxi) 
          trg = &pairmtx[maxi][i];
        else continue;

        if (i < maxj) 
          src = &pairmtx[i][maxj];
        else if (i > maxj) 
          src = &pairmtx[maxj][i];
        else continue;
          
        if (bc->treeMethod == UPGMA) 
          *trg = (*trg + *src) / 2.0;
        else
          *trg = (*trg + *src - maxid) / 2.0;
        }

      /* Create node for maxi and maxj */
      newnode = g_malloc(sizeof(TreeNode));

      if (bc->treeMethod == UPGMA)
        {
          /* subtract lower branch lengths from absolute distance
             Horribly ugly, only to be able to share code UPGMA and NJ */
          TreeNode *tmpnode = node[maxi]->left;

          for (llen = maxid; tmpnode; tmpnode = tmpnode->left)
            llen -= tmpnode->branchlen;

          tmpnode = node[maxj]->right;
          
          for (rlen = maxid; tmpnode; tmpnode = tmpnode->right)
            rlen -= tmpnode->branchlen;
        }
      else
        {
          llen = (maxid + avgdist[maxi] - avgdist[maxj]) / 2.0;
          rlen = maxid - llen;
          
          if (iter == bc->alignArr->len - 2)
            { /* Root node */
              
              /* Not necessary anymore, the tree is re-balanced at the end which calls this too
                 treeBalanceByWeight(node[maxi], node[maxj], &llen, &rlen);*/
      	
              /* Put entire length of root branch in one leg so the rebalancing 
                 will work properly (otherwise it is hard to take this branch into account */
              rlen += llen;
              llen = 0;
            }
        }
      
#ifdef DEBUG
        {
          printf("Iter %d: Merging %d and %d, dist= %f\n", 
                 iter, maxi+1, maxj+1, mtx[maxi][maxj]);
          printf("maxid= %f  llen= %f  rlen= %f\n", maxid, llen, rlen);
          printf("avgdist[left]= %f  avgdist[right]= %f\n\n", 
                 avgdist[maxi], avgdist[maxj]);
        }
#endif

      newnode->left = node[maxi];	
      newnode->left->branchlen = llen;

      newnode->right = node[maxj];
      newnode->right->branchlen = rlen;

      newnode->organism = (node[maxi]->organism == node[maxj]->organism ?
                           node[maxi]->organism : 0);
      
      node[maxi] = newnode;
      node[maxj] = 0;
    }

  fillParents(newnode, newnode->left);
  fillParents(newnode, newnode->right);

#ifdef DEBUG
    treeDraw(newnode) ;
#endif

  if (bc->treeMethod == UPGMA)
    newnode->branchlen = 100 - maxid ;

  if (bc->treeMethod == NJ)
    newnode = treeFindBalance(bc, newnode) ;
  
  fillOrganism(newnode);

  if (doBootstrap && bc->treebootstraps) 
    treeBootstrapStats(bc, newnode);
  
  return newnode ;
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

  bc->treeHead = treeMake(bc, FALSE);

  treeOrder(bc->treeHead, 1); /* Set nr field according to tree order */

  g_array_sort(bc->alignArr, nrorder);
  
  reInsertMarkupLines(bc);
}


void highlightScoreSort(char mode, BelvuContext *bc)
{
  if (!bc->highlightedAln) 
    {
      g_critical("Please highlight a sequence first");
      return;
    }
    
  if (bc->highlightedAln->markup) 
    {
      g_critical("Please do not highlight a markup line");
      return;
    }

  ALN aln;
  alncpy(&aln, bc->highlightedAln);
    
  separateMarkupLines(bc);

  /* if (displayScores) {
     if (!(graphQuery("This will erase the current scores, do you want to continue?")))
	    return;
            }*/
    
  bc->displayScores = 1;

  /* Calculate score relative to highlighted sequence */
  int i = 0;
  
  for (i = 0; i < bc->alignArr->len; ++i)
    {
      if (mode == 'P')
        {
          ALN *curAln = &g_array_index(bc->alignArr, ALN, i);
          curAln->score = score(aln.seq, curAln->seq, bc->penalize_gaps);
        }
      else if (mode == 'I')
        {
          ALN *curAln = &g_array_index(bc->alignArr, ALN, i);
          curAln->score = identity(aln.seq, curAln->seq, bc->penalize_gaps);
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
      g_critical("Cannot find back highlighted seq after sort - probably a bug");
      bc->highlightedAln = 0;
    }
  else
    {
      bc->highlightedAln = &g_array_index(bc->alignArr, ALN, i);
    }
  
  arrayOrder(bc->alignArr);
  bc->alignYStart = 0;
  //belvuRedraw();
}


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
      src = alni->seq;

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

      g_free(alni->seq);
      alni->seq = seq;
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
  SEG	*SegList, *seg;
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
	  resetALN(&aln);

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

	  aln.len = strlen(seq);
	  aln.seq = seq;
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


void checkAlignment(BelvuContext *bc)
{
  int i, g, cres, nres, tmp;

  bc->maxNameLen = bc->maxStartLen = bc->maxEndLen = 0;

  for (i = 0; i < bc->alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);

      if (!alnp->markup) 
        {
          /* Count residues */
          for (cres = g = 0; g < bc->maxLen; g++) 
            {
              if (isAlign(alnp->seq[g]) && !isGap(alnp->seq[g])) 
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
          bc->conservCount[a2b[(unsigned char)(alnp->seq[i])]][i]++;

          if (isalpha(alnp->seq[i]) || alnp->seq[i] == '*')
            bc->conservResidues[i]++;
        }
    }

  return nseqeff;
}


/* Return 1 if c1 has priority over c2, 0 otherwise */
static int colorPriority(BelvuContext *bc, int c1, int c2) 
{
  if (c2 == WHITE) return 1;
  if (c2 == bc->maxbgColor) return 0;
  if (c2 == bc->lowbgColor) {
    if (c1 == bc->lowbgColor) return 0;
    else return 1;
  }
  if (c2 == bc->midbgColor) {
    if (c1 == bc->maxbgColor) return 1;
    else return 0;
  }
  
  g_error("Unknown colour %s", colorNames[c2]);		    /* exits program. */

  return 0 ;
}


void setConservColors(BelvuContext *bc)
{
  int i, j, k, l, colornr, simCount, n, nseqeff;
  double id, maxid;

  if (!bc->conservCount) 
    initConservMtx(bc);

  nseqeff = countResidueFreqs(bc);

  for (i = 0; i < bc->maxLen; ++i)
    for (k = 1; k < 21; ++k)
      bc->colorMap[k][i] = WHITE;

  for (i = 0; i < bc->maxLen; ++i) 
    {
      maxid = -100.0;
      for (k = 1; k < 21; k++) 
        {
          if (bc->color_by_similarity) 
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
                id = 0.0;
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
                    colornr = bc->maxbgColor;
                  else if (id > bc->midSimCutoff) colornr = bc->midbgColor;
                  else colornr = bc->lowbgColor;
		    
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
		
              if (bc->colorByResIdOn) 
                {
                  if (id*100.0 > bc->colorByResIdCutoff)
                    bc->colorMap[k][i] = color[(unsigned char)(b2a[k])];
                  else
                    bc->colorMap[k][i] = WHITE;
		}
              else if (id > bc->lowIdCutoff) 
                {
                  if (id > bc->maxIdCutoff) 
                    colornr = bc->maxbgColor;
                  else if (id > bc->midIdCutoff) 
                    colornr = bc->midbgColor;
                  else
                    colornr = bc->lowbgColor;
		    
                  bc->colorMap[k][i] = colornr;
		    
                  if (bc->id_blosum) 
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
            maxid = id;
	}
      
      bc->conservation[i] = maxid;
    }
}
	    

void initResidueColors(BelvuContext *bc)
{
  bc->colorScheme = COLORBYRESIDUE;

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


void readColorCodes(BelvuContext *bc, FILE *fil, int *colorarr)
{
  char *cp, line[MAXLINE+1], setColor[MAXLINE+1];
  unsigned char c ;
  int i, colornr;

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

  bc->color_by_conserv = 0;
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


/* Removes all columns between from and to, inclusively.
 * Note that from and to are sequence coordinates, 1...maxLen !!!!!  */
static void rmColumn(BelvuContext *bc, int from, int to)
{
  int
    i, j, len = to - from + 1;
  ALN 
    *alni;
  
  fprintf(stderr, "Removing Columns %d-%d.\n", from, to);

  for (i = 0; i < bc->alignArr->len; i++) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);

      /* If N or C terminal trim, change the coordinates */
      if (from == 1)
        for (j = from; j <= to;  j++) 
          {
            /* Only count real residues */
            if (!isGap(alni->seq[j-1]))
              (alni->start < alni->end ? alni->start++ : alni->start--); 
          }

      if (to == bc->maxLen)
        for (j = from; j <= to;  j++) 
          {
            /* Only count real residues */
            if (!isGap(alni->seq[j-1]))
              (alni->start < alni->end ? alni->end-- : alni->end++); 
          }
	
      /* Remove the columns */
      for (j = 0; j < bc->maxLen-to+1 /* Including terminator 0 */;  j++) 
        {
          alni->seq[from+j-1] = alni->seq[to+j];
	}
      j--;
      
      if (alni->seq[from+j-1] || alni->seq[to+j])
        printf("Still a bug in rmColumn(): End=%c, Oldend=%c\n", 
               alni->seq[from+j-1],
               alni->seq[to+j]);
    }

  bc->maxLen -= len;

  bc->saved = 0;
}


void rmEmptyColumns(BelvuContext *bc, double cutoff)
{
  int
    i, j, gaps, totseq, removed=0, oldmaxLen=bc->maxLen;
  ALN 
    *alni;
  char
    c;

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
          
          if (!alni->markup) 
            {
              c = alni->seq[j];
              if (isGap(c) || c == ' ') 
                gaps++;
	    }
	}
      
      if ((double)gaps/totseq >= cutoff) 
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


static void rmFinalise(BelvuContext *bc) 
{
  /*    ruler[maxLen] = 0;*/
  checkAlignment(bc);
  setConservColors(bc);
  //  belvuRedraw();
}


/* Get rid of seqs that start or end with a gap.
 */
void rmPartial(BelvuContext *bc)
{
  int i, n=0;
  ALN *alni;

  for (i = 0; i < bc->alignArr->len; ) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);

      if (isGap(alni->seq[0]) || isGap(alni->seq[bc->maxLen-1])) 
        {
          /* Remove entry */
          n++;

          if (bc->highlightedAln == alni) 
            bc->highlightedAln = 0;
          
          g_array_remove_index(bc->alignArr, alni->nr);
          g_array_sort(bc->alignArr, nrorder);

          bc->saved = 0;
	}
      else i++;
    }
  
  fprintf(stderr, "%d partial sequences removed.  %d seqs left\n\n", n, bc->alignArr->len);
  
  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


/* Get rid of seqs that are too gappy.
 */
void rmGappySeqs(BelvuContext *bc, double cutoff)
{
  int i, j, n=0, gaps;
  ALN *alni;

  for (i = 0; i < bc->alignArr->len; ) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);

      for (gaps = j = 0; j < bc->maxLen; j++)
        if (isGap(alni->seq[j])) 
          gaps++;

      if ((double)gaps/bc->maxLen >= cutoff/100.0) 
        {
          /* Remove entry */
          n++;
          
          if (bc->highlightedAln == alni) 
            bc->highlightedAln = 0;
          
          g_array_remove_index(bc->alignArr, i);
          g_array_sort(bc->alignArr, nrorder);
          
          bc->saved = 0;
	}
      else 
	{
	  i++;
	}
    }

    g_message("%d gappy sequences removed.  %d seqs left\n\n", n, bc->alignArr->len);

    arrayOrder(bc->alignArr);
}


void rmFinaliseGapRemoval(BelvuContext *bc)
{
  if (bc->rmEmptyColumnsOn) 
    rmEmptyColumns(bc, 1.0);
  
  rmFinalise(bc);
}


/* Get rid of seqs that are more than x% identical with others. 
 * Keep the  first one.
 */
void mkNonRedundant(BelvuContext *bc, double cutoff)
{
  int i,j, n=0, overhang;
  ALN *alni, *alnj;
  double id;

  for (i = 0; i < bc->alignArr->len; i++) 
    {
      alni = &g_array_index(bc->alignArr, ALN, i);

      for (j = 0; j < bc->alignArr->len; j++) 
        {
          if (i == j) 
            continue;

          alnj = &g_array_index(bc->alignArr, ALN, j);
	    
          overhang = alnOverhang(alnj->seq, alni->seq);
          id = identity(alni->seq, alnj->seq, bc->penalize_gaps);

          if (!overhang && id > cutoff)
	    {
              fprintf(stderr, "%s/%d-%d and %s/%d-%d are %.1f%% identical. "
                      "The first includes the latter which was removed.\n",
                      alni->name, alni->start, alni->end,
                      alnj->name, alnj->start, alnj->end,
                      id);
		
              /* Remove entry j */
              n++;
		
              if (bc->highlightedAln == alnj) 
                bc->highlightedAln = NULL;
              
              g_array_remove_index(bc->alignArr, alnj->nr);
              g_array_sort(bc->alignArr, nrorder);

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

  fprintf(stderr, "%d sequences removed at the %.0f%% level.  %d seqs left\n\n", n, cutoff, bc->alignArr->len);
  
  arrayOrder(bc->alignArr);
  rmFinaliseGapRemoval(bc);
}


int GCGchecksum(BelvuContext *bc, char *seq)
{
  int  
    check = 0, 
    count = 0, 
    i;
  
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
    i,
    grand_checksum=0;

  for(i=0; i < bc->alignArr->len; ++i) 
    grand_checksum += GCGchecksum(bc, g_array_index(bc->alignArr, ALN, i).seq);

  return (grand_checksum % 10000);
}


void writeMSF(BelvuContext *bc, FILE *pipe) /* c = separator between name and coordinates */
{
  int 
    i, j,
    maxfullnamelen = bc->maxNameLen + 1 + bc->maxStartLen + 1 + bc->maxEndLen,
    paragraph=0, alnstart, alnend, alnlen, linelen=50, blocklen=10;
  
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
              GCGchecksum(bc, alnp->seq),
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
              fprintf(pipe, "%c", alnp->seq[j]);
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
}


/* Utility to ask the user for a file to save to. Returns the file name */
static const char* getSaveFileName(GtkWidget *widget, 
                                   const char *currentName, 
                                   const char *defaultPath,
                                   const char *title)
{
  const char *filename = NULL;
  
  GtkWindow *window = GTK_WINDOW(gtk_widget_get_toplevel(widget));

  GtkWidget *dialog = gtk_file_chooser_dialog_new (title,
                                                   window,
                                                   GTK_FILE_CHOOSER_ACTION_SAVE,
                                                   GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                                   GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
                                                   NULL);
  
  gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (dialog), TRUE);
  
  if (defaultPath)
    {
      gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (dialog), defaultPath);
    }
  
  if (currentName)
    {
      gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER (dialog), currentName);
    }
  
  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
    {
      filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    }
  
  gtk_widget_destroy (dialog);
  return filename;
}


static void saveMSF(BelvuContext *bc)
{
  /* to do: pass parent widget */
  const char *filename = getSaveFileName(NULL, bc->fileName, bc->dirName, "Save as MSF (/) file:");
  
  FILE *fil = fopen(filename, "w");
  
  if (fil)
    {
      writeMSF(bc, fil);
      bc->saved = TRUE;
    }
}


static void readMSF(BelvuContext *bc, FILE *pipe)
{
  char seq[1001], *cp, *cq;
  int i;
  ALN aln;
  char line[MAXLENGTH + 1];

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

          aln.len = len;
          aln.seq = g_malloc(aln.len+1);
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
              ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
              strcat(alnp->seq, seq);
	    }
          else
            {
              g_critical("Cannot find back %s %d %d seq=%s.\n", 
                         aln.name, aln.start, aln.end, seq);
	    }
	}
    }
  
  strcpy(bc->saveFormat, MSFStr);
}


/* Called from readMul.  Does the work to parse a selex file. */
static void readSelex(BelvuContext *bc, FILE *pipe)
{
  char   *cp, *cq, ch = '\0';
  int    len, i, alnstart;
  
  char line[MAXLENGTH+1];

  /* Read raw alignment into stack
   *******************************/
  
//  alnstart = MAXLENGTH;
//  while (!feof (pipe))
//    { 
//      if (!fgets (line, MAXLENGTH, pipe) || (unsigned char)*line == (unsigned char)EOF ) break;
//      /* EOF checking to make acedb calling work */
//      
//      if (!strncmp(line, "PileUp", 6)) {
//        readMSF(bc, pipe);
//        return;
//      }
//      
//      if ((cp = strchr(line, '\n'))) *cp = 0;
//      len = strlen (line);
//      
//      /* Sequences */
//      if (len && *line != '#' && strcmp(line, "//"))
//	{
//          if (!(cq = strchr(line, ' '))) g_error("No spacer between name and sequence in %s", line);
//          
//          /* Find leftmost start of alignment */
//          for (i = 0; line[i] && line[i] != ' '; i++);
//          for (; line[i] && !isAlign(line[i]); i++);
//          if (i < alnstart) alnstart = i;
//          
//          /* Remove optional accession number at end of alignment */
//          /* FOR PRODOM STYLE ALIGNMENTS W/ ACCESSION TO THE RIGHT - MAYBE MAKE OPTIONAL
//           This way it's incompatible with alignments with ' ' gapcharacters */
//          for (; line[i] && isAlign(line[i]); i++);
//          line[i] = 0;
//          
//          pushText(stack, line);
//	}
//      /* Markup line */
//      else if (!strncmp(line, "#=GF ", 5) || 
//               !strncmp(line, "#=GS ", 5)) {
//        pushText(AnnStack, line);
//      }
//      else if (!strncmp(line, "#=GC ", 5) || 
//               !strncmp(line, "#=GR ", 5) || 
//               !strncmp(line, "#=RF ", 5)) {
//        pushText(stack, line);
//      }
//      /* Match Footer  */
//      else if (!strncmp(line, "# matchFooter", 13)) {
//        matchFooter = 1;
//        break;
//      }
//    }
//  
//  /* Read alignment from Stack into Align array
//   ***************************************/
//  
//  /* 
//   * First pass - Find number of unique sequence names and the total sequence lengths 
//   */
//  stackCursor(stack, 0);
//  while ((cp = stackNextText(stack)))
//    {
//      parseMulLine(bc, cp, &aln);
//      
//      /* Sequence length */
//      len = strlen(cp+alnstart);
//      
//      if (arrayFind(Align, &aln, &ip, (void*)alphaorder))
//        arrp(Align, ip, ALN)->len += len;
//      else
//        {
//          aln.len = len;
//          aln.nr = ++nseq;
//          arrayInsert(Align, &aln, (void*)alphaorder);
//          
//          /*printf("\nInserted %s %6d %6d %s %d\n", aln.name, aln.start, aln.end, cp+alnstart, len);
//           fflush(stdout);*/
//	}
//    }
//  
//  /* Find maximum length of alignment; allocate it */
//  for (i = 0; i < nseq; i++)
//    {
//      alnp = arrp(Align, i, ALN);
//      if (alnp->len > maxLen) maxLen = alnp->len;
//    }
//  
//  for (i = 0; i < nseq; i++)
//    {
//      alnp = arrp(Align, i, ALN);
//      alnp->seq = messalloc(maxLen+1);
//    }
//  
//  /* 
//   * Second pass - allocate and read in sequence data 
//   */
//  stackCursor(stack, 0);
//  while ((cp = stackNextText(stack)))
//    {
//      parseMulLine(bc, cp, &aln);
//      
//      if (arrayFind(Align, &aln, &i, (void*)alphaorder))
//        {
//          alnp = arrp(Align, i, ALN);
//          strcat(alnp->seq, cp+alnstart);
//        }
//      else
//        {
//          messout("Cannot find back %s %d %d seq=%s.\n", 
//                  aln.name, aln.start, aln.end, aln.seq);
//        }
//    }
//  
//  /* Parse sequence attributes from #=GS lines */
//  stackCursor(AnnStack, 0);
//  while ((cp = stackNextText(AnnStack)))
//    {
//      char
//      *namep,		/* Seqname */
//      name[MAXNAMESIZE+1],/* Seqname */
//      *labelp,		/* Label (OS) */
//      *valuep;		/* Organism (c. elegans) */
//      
//      if (strncmp(cp, "#=GS", 4))
//        continue;
//      
//      strcpy(name, cp+4);
//      namep = name;
//      while (*namep == ' ') namep++;
//      
//      labelp = namep;
//      while (*labelp != ' ') labelp++;
//      *labelp = 0;		/* Terminate namep */
//      labelp++;		/* Walk across terminator */
//      while (*labelp == ' ') labelp++;
//      
//      valuep = labelp;
//      while (*valuep != ' ') valuep++;
//      while (*valuep == ' ') valuep++;
//      
//      if (!strncasecmp(labelp, "LO", 2))
//        {
//          int i, colornr;
//          
//          for (colornr = -1, i = 0; i < NUM_TRUECOLORS; i++)
//            if (!strcasecmp(colorNames[i], valuep)) colornr = i;
//          
//          if (colornr == -1)
//            {
//              printf("Unrecognized color: %s\n", valuep);
//              colornr = 0;
//            }
//          
//          str2aln(namep, &aln);
//          /* Find the corresponding sequence */
//          if (!arrayFind(Align, &aln, &i, (void*)alphaorder))
//            {
//              messout("Cannot find %s %d %d in alignment", aln.name, aln.start, aln.end);
//              continue;
//            }
//          alnp = arrp(Align, i, ALN);
//          alnp->color = colornr;
//        }
//      
//      else if (!strncasecmp(labelp, OrganismLabel, 2))
//        {
//          
//          /* Add organism to table of organisms.  This is necessary to make all 
//           sequences of the same organism point to the same place and to make a 
//           non-redundant list of present organisms */
//          
//          /* Find string in permanent stack */
//          if (!(valuep = strstr(cp, valuep)))
//            {
//              messout("Failed to parse organism properly");
//              continue;
//            }
//          aln.organism = valuep; /* Find string in permanent stack */
//          arrayInsert(organismArr, &aln, (void*)organism_order);
//          
//          if (strchr(cp, '/') && strchr(cp, '-'))
//            {
//              
//              str2aln(namep, &aln);
//              
//              /* Find the corresponding sequence */
//              if (!arrayFind(Align, &aln, &i, (void*)alphaorder)) {
//                messout("Cannot find %s %d %d in alignment", aln.name, aln.start, aln.end);
//                continue;
//              }
//              alnp = arrp(Align, i, ALN);
//              
//              /* Store pointer to unique organism in ALN struct */
//              aln.organism = valuep;
//              arrayFind(organismArr, &aln, &ip, (void*)organism_order);
//              alnp->organism = arrp(organismArr, ip, ALN)->organism;
//            }
//          else {
//            /* Organism specified for entire sequences.
//             Find all ALN instances of this sequences.
//             */
//            for (i = 0; i < nseq; i++) {
//              alnp = arrp(Align, i, ALN);
//              if (!strcmp(alnp->name, namep)) {
//                
//                /* Store pointer to unique organism in ALN struct */
//                aln.organism = valuep;
//                arrayFind(organismArr, &aln, &ip, (void*)organism_order);
//                alnp->organism = arrp(organismArr, ip, ALN)->organism;
//              }
//            }
//          }
//	}
//    }
//  
//  
//  /* For debugging * /
//   for (i = 0; i < nseq; i++) {
//   alnp = arrp(Align, i, ALN);
//   printf("\n%-10s %4d %4d %s %d %s\n", 
//   alnp->name, alnp->start, alnp->end, alnp->seq, 
//   alnp->len, alnp->organism);
//   }
//   for (i=0; i < arrayMax(organismArr); i++) 
//   printf("%s\n", arrp(organismArr, i, ALN)->organism);
//   */
//  
//  
//  if (!nseq || !maxLen) g_error("Unable to read sequence data");
//  
//  strcpy(saveFormat, MulStr);
}


void writeMul(BelvuContext *bc, FILE *fil)
{
//  static char *cp;
//  int i, w, W;
//
//  /* Write Annotation */
//  if (!stackEmpty(AnnStack))
//    {
//      stackCursor(AnnStack, 0);
//      while ((cp = stackNextText(AnnStack)))
//        {
//          fprintf(fil, "%s\n", cp);
//        }
//    }
//
//  /* Find max width of name column */
//  for (i = w = 0; i < bc->alignArr->len; i++)
//    if ( (W = strlen(writeMulName(g_array_index(bc->alignArr, ALN, i)))) > w) 
//      w = W;
//    
//  /* Write alignment */
//  for (i = 0; i < bc->alignArr->len; i++)
//    fprintf(fil, "%-*s %s\n", w, writeMulName(g_array_index(bc->alignArr, ALN, i)), g_array_index(bc->alignArr, ALN, i)->seq);
//  
//  fprintf(fil, "//\n");
//  
//  fclose(fil);
//  fflush(fil);
//
//  bc->saved = 1;
}


static void saveMul(BelvuContext *bc)
{
  /* to do: pass parent widget */
  const char *filename = getSaveFileName(NULL, bc->fileName, bc->dirName, "Save as Stockholm file:");
  FILE *fil = fopen(filename, "w");

  if (fil)
    {
      writeMul(bc, fil);
    }
}


/* ReadMul 
 * parses an alignment file in mul or selex format and puts it in the Align array
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
void readMul(BelvuContext *bc, FILE *pipe)
{
  char   ch = '\0';
  char line[MAXLENGTH+1];
  
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
      else if (ch != '\n')
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

  return readSelex(bc, pipe);
}


void writeFasta(BelvuContext *bc, FILE *pipe)
{
  int i, n;
  char *cp;

  for (i = 0; i < bc->alignArr->len; ++i) 
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (alnp->markup) 
        continue;
	
      if (bc->saveCoordsOn)
        fprintf(pipe, ">%s%c%d-%d\n", alnp->name, bc->saveSeparator, alnp->start, alnp->end);
      else
        fprintf(pipe, ">%s\n", alnp->name);
	
      for (n=0, cp = alnp->seq; *cp; cp++)
        {
          if (!strcmp(bc->saveFormat, FastaAlnStr)) 
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
}


static void saveFasta(BelvuContext *bc)
{
  /* to do: pass parent widget */
  const char *filename = getSaveFileName(NULL, bc->fileName, bc->dirName, "Save as unaligned Fasta file:");
  
  FILE *fil = fopen(filename, "w");
  
  if (fil)
    {
      writeFasta(bc, fil);
      bc->saved = TRUE;
    }
}


static int getMatchStates(BelvuContext *bc)
{
  int i, j, n, retval=0;

  for (j = 0; j < bc->maxLen; ++j) 
    {
      n = 0;
      for (i = 0; i < bc->alignArr->len; ++i) 
        {
          ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
          
          if (isalpha(alnp->seq[j])) 
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

  double p;
  int i, c[21], n=0, nmat;
  char *cp;

  for (i = 1; i <= 20; i++) 
    c[i] = 0;

  for (i = 0; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      cp = alnp->seq;

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
  static double dist;
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
  
  treeStruct->head = treeMake(bc, 0);
  
  treeTraverseLRfirst(bc, treeStruct->head, subfamilyTrav);
}


void showAnnotation(BelvuContext *bc)
{
//  int maxwidth=0,
//    nlines=0,
//    i=0;
//  char *cp;
//
//  if (stackEmpty(AnnStack)) 
//    return;
//
//  stackCursor(AnnStack, 0);
//  while ((cp = stackNextText(AnnStack))) 
//    {
//      nlines++;
//      if (strlen(cp) > maxwidth) 
//        maxwidth = strlen(cp);
//    }
//
//  if (!graphActivate(annGraph)) 
//    {
//      annGraph = graphCreate (TEXT_FULL_SCROLL, "Annotations", 0, 0,
//                              (maxwidth+1)/fontwidth*screenWidth, 
//                              (nlines+1)/fontheight*screenHeight);
//      graphTextFormat(FIXED_WIDTH);
//    }
//  
//  graphClear();
//  graphTextBounds(maxwidth+1, nlines+1);
//
//  stackCursor(AnnStack, 0);
//
//  while ((cp = stackNextText(AnnStack))) 
//    {
//      graphText(cp, 0, i++);
//    }
//  
//  graphRedraw();
}


void treeDisplay(BelvuContext *bc)
{
  separateMarkupLines(bc);
  bc->treeHead = treeMake(bc, 1);
  treeDraw(bc, bc->treeHead);
  reInsertMarkupLines(bc);
}


static void colorCons(BelvuContext *bc)
{
  setConservColors(bc);

//  menuSetFlags(menuItem(colorMenu, thresholdStr), MENUFLAG_DISABLED);
//  menuUnsetFlags(menuItem(colorMenu, printColorsStr), MENUFLAG_DISABLED);
//  menuUnsetFlags(menuItem(colorMenu, ignoreGapsStr), MENUFLAG_DISABLED);

  bc->colorByResIdOn = FALSE;
  //  belvuRedraw();
}


void colorSim(BelvuContext *bc)
{
  bc->colorScheme = COLORSIM;
  bc->color_by_conserv = 1;
  bc->color_by_similarity = 1;
  colorCons(bc);
}


/* Create the context, which contains all program-wide variables */
BelvuContext* createBelvuContext()
{
  BelvuContext *bc = g_malloc(sizeof *bc);
  
  bc->defaultColors = NULL;
  
  bc->alignArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN), 100);  /* was called 'Align' */
  bc->organismArr = g_array_sized_new(FALSE, FALSE, sizeof(ALN), 100);
  bc->markupAlignArr = NULL;
  bc->bootstrapGroups = NULL;

  bc->highlightedAln = NULL;

  bc->treeHead = NULL;
  bc->treeBestBalancedNode = NULL;

  bc->treeReadDistancesPipe = NULL;

  bc->treeMethod = NJ;
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
  bc->colorScheme = COLORSIM;
  
  bc->maxfgColor = BLACK;
  bc->midfgColor = BLACK,
  bc->lowfgColor = BLACK;
  bc->maxbgColor = CYAN;
  bc->midbgColor = MIDBLUE;
  bc->lowbgColor = LIGHTGRAY;

  bc->treeDistCorr = SCOREDIST;
  bc->treeBestBalance = 0.0;
  bc->treeBestBalance_subtrees = 0.0;
  bc->tree_y = 0.3;
  bc->lowIdCutoff = 0.4;
  bc->midIdCutoff = 0.6;
  bc->maxIdCutoff = 0.8;
  bc->lowSimCutoff = 0.5;
  bc->midSimCutoff = 1.5;
  bc->maxSimCutoff = 3.0;
  bc->colorByResIdCutoff = 20.0;
  bc->mksubfamilies_cutoff = 0.0;

  strcpy(bc->treeDistString, SCOREDISTstr);
  strcpy(bc->treeMethodString, NJstr);
  
  bc->saveSeparator = '/';
  bc->Title[0] = '\0';
  bc->saveFormat[0] = '\0';
  
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
  bc->displayScores = TRUE;
  bc->outputBootstrapTrees = FALSE;
  bc->treeColorsOn = TRUE;
  bc->treeShowOrganism = TRUE;
  bc->treeShowBranchlen = FALSE;
  bc->matchFooter = FALSE;
  bc->saved = TRUE;
  bc->color_by_similarity = TRUE;
  bc->color_by_conserv = TRUE;
  bc->ignoreGapsOn = FALSE;
  bc->colorByResIdOn = FALSE;
  bc->id_blosum = TRUE;
  bc->rmEmptyColumnsOn = TRUE;
  bc->lowercaseOn = FALSE;
  
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


/* This function just returns the value in the b2a array at the given index */
char b2aIndex(const int idx)
{
  return b2a[idx];
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


