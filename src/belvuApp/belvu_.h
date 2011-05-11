/*  File: belvu_.h
 *  Author: Gemma Barson, 2011-03-06
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
 * Description: Internal header for Belvu package-wide code
 *----------------------------------------------------------------------------
 */

#ifndef DEF_BELVU_P_H
#define DEF_BELVU_P_H

#include <gtk/gtk.h>
#include <config.h>
#include <seqtoolsUtils/utilities.h>
#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/version.h>
#include <sys/param.h>


/*            blixem program version and information.                        */
#define BELVU_TITLE   "Belvu"
#define BELVU_DESC    "Multiple alignment visualisation tool."

/* The Seqtools package version should be specified in src/version.m4. autoconf will then set PACKAGE_VERSION in config.h */
#define BELVU_VERSION_STRING	   PACKAGE_VERSION
#define BELVU_PACKAGE_VERSION	   UT_MAKE_VERSION_INFO_STRING(PACKAGE_NAME, PACKAGE_VERSION)
#define BELVU_TITLE_STRING	   UT_MAKE_TITLE_STRING(BELVU_TITLE, BELVU_VERSION_STRING)
#define BELVU_VERSION_COMPILE	   BELVU_VERSION_STRING "  " UT_MAKE_COMPILE_DATE()

#define BELVU_COPYRIGHT_STRING	   UT_MAKE_COPYRIGHT_STRING("2009-2010")
#define BELVU_WEBSITE_STRING	   "http://www.sanger.ac.uk/resources/software/seqtools/"
#define BELVU_LICENSE_STRING	   UT_MAKE_LICENCE_STRING(BELVU_TITLE)

#define BELVU_COMMENTS_STRING()                                \
"("BELVU_TITLE_STRING", "					\
UT_COMPILE_PHRASE " " UT_MAKE_COMPILE_DATE() ")\n\n"            \
AUTHOR_TEXT "\n"


/* _MAX_PATH is 260 in WIN32 but each path component can be max. 256 in size */
#define FIL_BUFFER_SIZE 256
#define DIR_BUFFER_SIZE MAXPATHLEN

#define MAXNAMESIZE  256
#define MAXLENGTH 100000


#define UPGMAstr "UPGMA"
#define NJstr "NJ"

#define MSFStr "MSF"
#define MulStr "Stockholm (Pfam/HMMER)"
#define FastaAlnStr "Aligned Fasta"
#define FastaStr "Unaligned Fasta"

#define SWAPstr "Swap"
#define ROOTstr "Reroot"

#define UNCORRstr "Uncorrected differences"
#define KIMURAstr "Kimura correction"
#define JUKESCANTORstr "Jukes-Cantor correction"
#define STORMSONNstr "Storm & Sonnhammer correction"
#define SCOREDISTstr "Scoredist correction"

#define NOCOLOR -1
#define NN NOCOLOR
#define BG WHITE
#define BOXCOLOR  WHITE /* LIGHTGRAY */



enum Colour    
{
  WHITE, BLACK, LIGHTGRAY, DARKGRAY,
  RED, GREEN, BLUE,
  YELLOW, CYAN, MAGENTA,
  LIGHTRED, LIGHTGREEN, LIGHTBLUE,
  DARKRED, DARKGREEN, DARKBLUE,
  PALERED, PALEGREEN, PALEBLUE,
  PALEYELLOW, PALECYAN, PALEMAGENTA,
  BROWN, ORANGE, PALEORANGE,
  PURPLE, VIOLET, PALEVIOLET,
  GRAY, PALEGRAY,
  CERISE, MIDBLUE,
  NUM_TRUECOLORS,
  TRANSPARENT,	/* pseudocolour only for background */
  FORECOLOR,	/* pseudocolor to force box->fcol after graphColor() */
  BACKCOLOR	/* pseudocolor to force box->bcol after graphColor() */
};


enum { MUL, RAW };		/* Input formats for IN_FORMAT */
enum { dummy, GF, GC, GS, GR };	/* Markup types in Stockholm format */

enum { COLORBYRESIDUE, COLORSCHEMESTANDARD, COLORSCHEMEGIBSON, COLORSCHEMECYS, 
       COLORSCHEMEEMPTY, COLORSIM, COLORID, COLORIDSIM }; /* Modes - used for checkmarks */

enum { NOLL, SORT_ALPHA, SORT_ORGANISM, SORT_TREE, SORT_SCORE, SORT_SIM, SORT_ID }; /* Initial sorting modes */

/* Tree picking methods */
typedef enum _BelvuPickMode
{
  NODESWAP,
  NODEROOT
} BelvuPickMode;


/* Tree building methods */
typedef enum _BelvuBuildMethod
{
  UPGMA,
  NJ
} BelvuBuildMethod;


/* Distance correction methods for tree building */
typedef enum _BelvuDistCorr
{
  UNCORR,
  KIMURA,
  STORMSONN,
  SCOREDIST,
  JUKESCANTOR
} BelvuDistCorr; 


/* The following are used to define default colors for certain types of features in Belvu.
 * One of several different actual colors from the BlxColor struct may be used depending 
 * on state, e.g. we use a different color if "print colors" (i.e. black and 
 * white mode) is on. */
typedef enum _BelvuColorId
  {
    BELCOLOR_MIN,                             /* dummy value so that we don't get a zero ID */
    
    BELCOLOR_BACKGROUND,                      /* background color for widgets */
    BELCOLOR_ALIGN_TEXT,                      /* text color for alignments */
    BELCOLOR_TREE_DEFAULT,                    /* default line color for the tree */
    BELCOLOR_TREE_BOOTSTRAP,
    
    BELCOLOR_NUM_COLORS
  } BelvuColorId;


/* This defines the possible sort orders */
typedef enum _BelvuSortType
  {
    BELVU_SORT_SCORE,          /* Sort by score */
    BELVU_SORT_ALPHA,          /* Sort alphabetically */
    BELVU_SORT_ORGANISM,       /* Sort by organism */
    BELVU_SORT_TREE,           /* Sort by tree order */
    BELVU_SORT_SIM,            /* Sort by similarity to selected sequence */
    BELVU_SORT_ID,             /* Sort by identity to selected sequence */
  } BelvuSortType;


/* This defines the types of color scheme that we deal with */
typedef enum _BelvuSchemeType
  {
    BELVU_SCHEME_TYPE_RESIDUE,          /* Color by residue */
    BELVU_SCHEME_TYPE_CONS              /* Color by conservation */
  } BelvuSchemeType;


/* This defines the color schemes when coloring by residue */
typedef enum _ResidueColorSchemes
  {
    BELVU_SCHEME_NONE,              /* Clean slate (no coloring) */
    BELVU_SCHEME_ERIK,              /* Erik's original scheme */
    BELVU_SCHEME_GIBSON,            /* Toby's */
    BELVU_SCHEME_CGP                /* Cyc/Gly/Pro */
  } ResidueColorSchemes;


/* This defines the color schemes when coloring by conservation */
typedef enum _ConsColorSchemes
  {
    BELVU_SCHEME_BLOSUM,            /* Average similarity by BLOSUM62 */
    BELVU_SCHEME_ID,                /* Percent identity */
    BELVU_SCHEME_ID_BLOSUM          /* Percent ID and BLOSUM62 */
  } ConsColorSchemes;



/* Structs */
typedef struct alnStruct {
  char name[MAXNAMESIZE + 1];
  int  start;
  int  end;
  int  len;			/* Length of this sequence - Do not use!  Use maxLen instead ! */
  char *seq;
  int  nr;			/* Ordinal number in array - must be sorted on this */
  char fetch[MAXNAMESIZE + 11];
  float score;
  int  color;			/* Background color of name */
  int  markup;		/* Markup line */
  gboolean  hide;			/* Hide this line */
  gboolean nocolor;		/* Exclude from coloring */
  char *organism;
} ALN;


typedef struct TreeNodeStruct {

    /* KEEP IN SYNC WITH TREECPY !!!! */

    double dist;			/* Absolute distance position */
    double branchlen;		/* Length of branch to higher node */
    double boot;	/* Bootstrap value */
    struct TreeNodeStruct *left;
    struct TreeNodeStruct *right;
    struct TreeNodeStruct *parent;
    char *name;
    char *organism;
    ALN *aln;
    int box;
    int color;

    /* KEEP IN SYNC WITH TREECPY !!!! */

} TreeNode;


typedef struct TreeStruct
{
   TreeNode *head;
   int lastNodeBox;		/* Last box used by a node - name boxes come next */
   int currentPickedBox;
} Tree;


/* Struct to store bootstrap group */
typedef struct BootstrapGroupStruct
{   
  TreeNode *node;    /* Points to node in original tree (for incrementing) */
  char *s;           /* Sorted concatenation of all sequences in node, to be inserted in list */
} BootstrapGroup;


typedef struct SegStruct
{
  int  start;
  int  end;
  int  qstart;
  int  qend;
  struct SegStruct *next;
} SEG;


/* This enum contains IDs for all the persistent dialogs in the application, and should be used
 * to access a stored dialog in the dialogList array in the BelvuContext. Note that the dialogList
 * array will contain null entries until the dialogs are created for the first time */
typedef enum
  {
    BELDIALOG_NOT_PERSISTENT = 0,   /* Reserved for dialogs that do not have an entry in the array */
    
    BELDIALOG_MAKE_TREE,            /* The make-tree dialog */
    
    BELDIALOG_NUM_DIALOGS           /* The number of dialogs. Must always be the last entry in this enum */
  } BelvuDialogId;


typedef struct BelvuContextStruct
{
  GArray *defaultColors;            /* Default colors used by Belvu */
  
  GArray *alignArr;
  GArray *organismArr;
  GArray *markupAlignArr;
  GArray *bootstrapGroups;

  ALN *highlightedAln;

  TreeNode *treeHead;              /* global current tree head */
  TreeNode *treeBestBalancedNode;

  FILE *treeReadDistancesPipe;

  int IN_FORMAT;
  int maxScoreLen;
  int alignYStart;
  int treebootstraps;              /* Number of bootstrap trees to be made */
  int maxLen;                      /* number of columns in alignment */
  int maxTreeWidth;
  int maxNameLen;                  /* Max string length of any sequence name */
  int maxStartLen;                 /* Max string length of any sequence start */
  int maxEndLen;                   /* Max string length of any sequence end */
  int maxScoreen;                  /* Max string length of any sequence score */
  int colorScheme;                 /* Current colour scheme mode (used for checkmarks) */

  int maxfgColor;
  int midfgColor;
  int lowfgColor;
  int maxbgColor;
  int midbgColor;
  int lowbgColor;
                
  BelvuBuildMethod treeMethod;  /* Default building method for trees */
  BelvuDistCorr treeDistCorr;   /* Default distance correction method for trees */
  BelvuPickMode treePickMode;   /* Default action when picking a node in a tree */
  
  double treeBestBalance;
  double treeBestBalance_subtrees;
  double tree_y;
  double lowIdCutoff;	           /* %id cutoff for lowest colour */
  double midIdCutoff;	           /* %id cutoff for medium colour */
  double maxIdCutoff;	           /* %id cutoff for maximum colour */
  double lowSimCutoff;	           /* %id cutoff for lowest colour */
  double midSimCutoff;	           /* %id cutoff for medium colour */
  double maxSimCutoff;	           /* %id cutoff for maximum colour */
  double colorByResIdCutoff;       /* Colour by residue + id cutoff */
  double mksubfamilies_cutoff; 
  double treeScale;                /* Default scale to use for drawing the tree */
  double treeLineWidth;            /* Default line width of the branch lines in trees */
  
  char saveSeparator;
  char treeDistString[50];
  char treeMethodString[50];
  char Title[256];
  char saveFormat[50];
  char fileName[FIL_BUFFER_SIZE + 1];
  char dirName[DIR_BUFFER_SIZE + 1]; /* Default directory for file browser */
  
  
  int **conservCount;              /* Matrix of conservation values  - 21 x maxLen */
  int **colorMap;                  /* Matrix of conservation colours - 21 x maxLen */
  int *conservResidues;            /* Array of number of residues present in each column */
  double *conservation;            /* The max conservation in each column [0..maxLen] */

  gboolean treeCoordsOn;
  gboolean treeReadDistancesOn;
  gboolean treePrintDistances;
  gboolean penalize_gaps;
  gboolean stripCoordTokensOn;
  gboolean saveCoordsOn;
  gboolean displayScores;
  gboolean outputBootstrapTrees;   /* Output the individual bootstrap trees */
  gboolean treebootstrapsDisplay;  /* Display bootstrap trees on screen */
  gboolean treeColorsOn;           
  gboolean treeShowOrganism;       /* whether to display the organism name in the tree */
  gboolean treeShowBranchlen;      /* whether to display the branch length in the tree */
  gboolean matchFooter;
  gboolean saved;
  gboolean color_by_similarity;    /* FALSE => by id */
  gboolean color_by_conserv;       /* TRUE => by id or similarity; FALSE => by residue  */
  gboolean ignoreGapsOn;
  gboolean colorByResIdOn;         /* colour by residue type above identity cutoff */
  gboolean id_blosum;              /* Paint similar residues too */
  gboolean rmEmptyColumnsOn;
  gboolean lowercaseOn;
  gboolean removingSeqs;	   /* Set to true if in the 'removing sequences' mode */
  
  GtkWidget *dialogList[BELDIALOG_NUM_DIALOGS];   /* Array of all the persistent dialogs in the application */
  
} BelvuContext;


/* Functions */
BelvuContext*                             createBelvuContext();
void                                      destroyBelvuContext(BelvuContext **bc);

void					  setTreeScaleCorr(BelvuContext *bc, const int treeMethod);
void					  setTreeScale(BelvuContext *bc, const double newScale) ;

gint                                      alphaorder(gconstpointer xIn, gconstpointer yIn);
gint                                      organism_order(gconstpointer xIn, gconstpointer yIn);
gint                                      organismorder(gconstpointer xIn, gconstpointer yIn);
gint                                      scoreorder(gconstpointer xIn, gconstpointer yIn);
gint                                      nrorder(gconstpointer xIn, gconstpointer yIn);

void                                      highlightScoreSort(char mode, BelvuContext *bc);
void                                      treeSortBatch(BelvuContext *bc);


void                                      arrayOrder(GArray *alignArr);
gboolean                                  alignFind(GArray *alignArr, ALN *obj, int *idx);
void                                      initAln(ALN *alnp);

void                                      setOrganismColors(GArray *organismArr);

void                                      parseMulLine(BelvuContext *bc, char *line, ALN *aln);

void                                      readMatch(BelvuContext *bc, FILE *fil);                
void                                      checkAlignment(BelvuContext *bc);
void                                      setConservColors(BelvuContext *bc);
void                                      initResidueColors(BelvuContext *bc);
void                                      initMarkupColors(void);              
void                                      readColorCodes(BelvuContext *bc, FILE *fil, int *colorarr);
void                                      mkNonRedundant(BelvuContext *bc, double cutoff);
void                                      rmPartialSeqs(BelvuContext *bc);         
void                                      rmEmptyColumns(BelvuContext *bc, double cutoff);
void                                      rmGappySeqs(BelvuContext *bc, double cutoff);
void                                      rmFinaliseGapRemoval(BelvuContext *bc);
void					  rmOutliers(BelvuContext *bc, const double cutoff);
void					  rmScore(BelvuContext *bc, const double cutoff);

void                                      readMul(BelvuContext *bc, FILE *pipe);
void                                      writeMul(BelvuContext *bc, FILE *fil);
void                                      writeFasta(BelvuContext *bc, FILE *pipe);
void                                      writeMSF(BelvuContext *bc, FILE *pipe);

void                                      separateMarkupLines(BelvuContext *bc);
void                                      reInsertMarkupLines(BelvuContext *bc);
TreeNode*                                 treeMake(BelvuContext *bc, const gboolean doBootstrap);
void                                      treePrintNH(Tree *tree, TreeNode *node, FILE *file);

void                                      outputProbs(BelvuContext *bc, FILE *fil);
void                                      mksubfamilies(BelvuContext *bc, double cutoff);        

void                                      treeDisplay(BelvuContext *bc);       

void                                      colorSim(BelvuContext *bc);          
void                                      showAnnotation(BelvuContext *bc);    

char                                      b2aIndex(const int idx);
int                                       getMarkupColor(const char inputChar);
int                                       getConservColor(BelvuContext *bc, const char inputChar, const int idx);
int                                       getColor(const char inputChar);

int*                                      getColorArray();
int*                                      getMarkupColorArray();

gboolean                                  isGap(char c);
int                                       strcmp_(gconstpointer xIn, gconstpointer yIn);
gboolean                                  arrayFind(GArray *a, void *s, int *ip, int (* orderFunc)(gconstpointer, gconstpointer));
GArray*                                   copyAlignArray(GArray *inputArr);
void                                      columnCopy(GArray *alignArrDest, int destIdx, GArray *alignArrSrc, int srcIdx);
double                                    identity(char *s1, char *s2, const gboolean penalize_gaps);

void                                      convertColorNumToGdkColor(const int colorNum, GdkColor *result);
void                                      drawText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const char *text);
void                                      drawIntAsText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const int value);


#endif /* DEF_BELVU_P_H */
