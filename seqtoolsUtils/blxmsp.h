/*
 *  blxmsp.h
 *  blixem
 *
 *  Created by Gemma Barson on 02/09/2010.
 *  Copyright 2010 Sanger Institute. All rights reserved.
 *
 */

#ifndef _blxmsp_included_
#define _blxmsp_included_

#include <SeqTools/utilities.h>
#include <stdlib.h>


//extern GArray    *fsArr;		   /* in dotter.c */

/* These are used by both blixem and dotter */
#define selectFeaturesStr          "Feature series selection tool"
#define XY_NOT_FILLED -1000        /* Magic value meaning "value not provided" */


/* Main Blixem error domain */
#define BLX_ERROR g_quark_from_string("Blixem")

/* Error codes */
typedef enum
{
  BLX_ERROR_SEQ_SEGMENT,	      /* error finding sequence segment */
  BLX_ERROR_EMPTY_STRING,           /* error code for when user entered a zero-length string */
  BLX_ERROR_STRING_NOT_FOUND,       /* error code for when a search string is not found */
  BLX_ERROR_SEQ_NAME_NOT_FOUND,     /* the sequence name(s) being searched for were not found */
  BLX_ERROR_SEQ_DATA_MISMATCH       /* same sequence was parsed more than once and data does not match */
} BlxError;



/* Supported types of MSP */
typedef enum 
{
  BLXMSP_INVALID,                /* No valid type was set */
  
  BLXMSP_MATCH,                  /* A match (i.e. alignment) */
  BLXMSP_MATCH_SET,              /* The parent of a set of matches. Can be used to specify generic properties such as color. */
  BLXMSP_CDS,                    /* CDS (coding) region of an exon */
  BLXMSP_UTR,                    /* UTR (untranslated) region of an exon */
  BLXMSP_INTRON,                 /* Intron */
  BLXMSP_EXON,			 /* Exon (should appear AFTER CDS and UTR for sorting, as required by constructTranscriptData) */
  BLXMSP_TRANSCRIPT,		 /* Transcript */
  BLXMSP_POLYA_SITE,		 /* polyA tail site */
  BLXMSP_POLYA_SIGNAL,		 /* polyA signal */
  
  BLXMSP_VARIATION,              /* SNP, substitution, deletion, insertion */
  
  BLXMSP_HSP,                    /*  */
  BLXMSP_GSP,                    /*  */
  
  BLXMSP_FS_SEG,                 /* Feature Series Segment */
  BLXMSP_XY_PLOT,                /* x/y coordinates - for plotting feature-series curves */
  
  BLXMSP_NUM_TYPES               /* the number of MSP types. MUST BE LAST IN LIST */
} BlxMspType;


/* Type definition for BlxSequences */
typedef enum
  {
    BLXSEQUENCE_UNSET,
    BLXSEQUENCE_TRANSCRIPT,         /* transcript (i.e. collection of exons and introns) */
    BLXSEQUENCE_MATCH,              /* match sequence (i.e. collection of matches) */
    BLXSEQUENCE_VARIATION           /* variation (i.e. insertion, deletion or substitution) */
  } BlxSequenceType;


/* Structure that contains information about a sequence */
typedef struct _BlxSequence
{
  BlxSequenceType type;            /* What type of collection of MSPs this is */

  char *idTag;			   /* Unique identifier e.g. from ID tag in GFF files */
  
  char *fullName;                  /* full name of the sequence and variant, including prefix characters, e.g. EM:AV274505.2 */
  char *shortName;                 /* short name of the sequence, excluding prefix and variant, e.g. AV274505 */
  char *variantName;               /* short name of the variant, excluding prefix but including variant number, e.g. AV274505.2 */
  
  GString *organism;               /* organism from the EMBL data OS line */
  GString *geneName;               /* gene name from the EMBL data GN line */
  GString *tissueType;             /* tissue type from the /tissue_type attribute from the FT lines in the EMBL file */
  GString *strain;                 /* strain from the /strain attribute from the FT lines in the EMBL file */
  
  BlxStrand strand;                /* which strand of the sequence this is */
  GString *sequence;               /* the actual sequence data */
  gboolean sequenceReqd;           /* whether the sequence data is required (e.g. it is not needed for exons/introns etc.) */
  IntRange qRange;		   /* the extent of the sequence on the ref sequence */ 
  
  GList *mspList;                  /* list of MSPs from this sequence */
} BlxSequence;


typedef struct _FeatureSeries
{
  char        *name;             /* Name of the Feature Series */
  int         on;                /* Flag to show/hide the Feature Series */
  int         order;             /* Order number used for sorting */
  float       x;	           /* Series offset on x axis, to bump series on the screen */
  float       y;	           /* Series offset on y axis */
  int         xy;	           /* Flag for XY plot series */
} FeatureSeries;


/* Shapes of XY curves */
typedef enum
{ 
  BLXCURVE_PARTIAL, 
  BLXCURVE_INTERPOLATE, 
  BLXCURVE_BADSHAPE
} BlxCurveShape;



/* Structure holding information about an alignment */
typedef struct _MSP
{
  struct _MSP       *next;
  GList             *childMsps;    /* Child MSPs of this MSP if it has them, e.g. an exon has CDS and UTR children (part_of relationship). */
  BlxMspType        type;          /* Whether this is a match, exon, SNP etc. */
  gdouble           score;         /* Score as a percentage. Technically this should be a weighted score taking into account gaps, length of the match etc., but for unknown reasons the ID has always been passed instead of score and the ID gets stored in here */
  gdouble           id;            /* Identity as a percentage. A simple comparison of bases within the match, ignoring gaps etc. Currently this is calculated internally by blixem. */
  int               phase;         /* phase: q start coord is offset by this amount to give the first base in the first complete codon (only relevant to CDSs) */
  char		    *url;
  
  char              *qname;        /* For Dotter, the MSP can belong to either sequence */
  IntRange	    qRange;	   /* the range of coords on the ref sequence where the alignment lies */
  BlxStrand         qStrand;       /* which strand on the reference sequence the match is on */
  int               qFrame;        /* which frame on the reference sequence the match is on */
  
  BlxSequence       *sSequence;    /* pointer to a struct holding info about the sequence/strand this match is from */
  char              *sname;        /* sequence name (could be different to the sequence name in the blxSequence e.g. exons have a postfixed 'x') */
  IntRange	    sRange;	   /* the range of coords on the match sequence where the alignment lies */
  
  char              *desc;         /* Optional description text for the MSP */
  char              *source;       /* Optional source text for the MSP */
  GSList            *gaps;         /* Array of "gaps" in this homolgy (this is a bit of a misnomer because the array
                                    * gives the ranges of the bits that align rather than the ranges of the gaps in between */
  
  BlxStyle          *style;        /* Specifies drawing style for this MSP, e.g. fill color and line color */
  
  FeatureSeries     *fs;           /* Feature series that this MSP belongs to */
  int               fsColor;       /* Color to draw this MSP in the feature series */
  BlxCurveShape     fsShape;       /* Shape data for drawing feature series curves, i.e. XY type PARTIAL or INTERPOLATE shapes */
  GArray            *xy;            /* For XY plot feature series */
} MSP ;


/* MSP functions */
gboolean              typeIsExon(const BlxMspType mspType);
gboolean              typeIsIntron(const BlxMspType mspType);
gboolean              typeIsMatch(const BlxMspType mspType);
gboolean              typeIsVariation(const BlxMspType mspType);

int		      mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType);
BlxStrand	      mspGetRefStrand(const MSP const *msp);
BlxStrand	      mspGetMatchStrand(const MSP const *msp);
const char*           mspGetMatchSeq(const MSP const *msp);
const char*	      mspGetSName(const MSP *msp);
char*                 mspGetExonTranscriptName(const MSP *msp);
const IntRange const* mspGetRefCoords(const MSP const *msp);
const IntRange const* mspGetMatchCoords(const MSP const *msp);
int		      mspGetQStart(const MSP const *msp);
int		      mspGetQEnd(const MSP const *msp);
int		      mspGetSStart(const MSP const *msp);
int		      mspGetSEnd(const MSP const *msp);
char*                 mspGetSSeq(const MSP const *msp);
int		      mspGetQRangeLen(const MSP const *msp);
int		      mspGetSRangeLen(const MSP const *msp);
int		      mspGetMatchSeqLen(const MSP const *msp);

const GdkColor*       mspGetColor(const MSP const *msp, 
				  GArray *defaultColors, const 
				  BlxSequence *blxSeq, 
				  const gboolean selected, 
				  const gboolean usePrintColors, 
				  const gboolean fill,
				  const int exonFillColorId,
				  const int exonLineColorId,
				  const int cdsFillColorId,
				  const int cdsLineColorId,
				  const int utrFillColorId,
				  const int utrLineColorId);

char*                 mspGetOrganism(const MSP const *msp);
char*                 mspGetGeneName(const MSP const *msp);
char*                 mspGetTissueType(const MSP const *msp);
char*                 mspGetStrain(const MSP const *msp);
char*                 mspGetCoordsAsString(const MSP const *msp);

gboolean              mspLayerIsVisible(const MSP const *msp);
gboolean	      mspIsExon(const MSP const *msp);
gboolean	      mspIsIntron(const MSP const *msp);
gboolean	      mspIsSnp(const MSP const *msp);
gboolean	      mspIsBlastMatch(const MSP const *msp);
gboolean	      mspIsPolyASite(const MSP const *msp);
gboolean	      mspIsVariation(const MSP const *msp);
gboolean	      mspIsZeroLenVariation(const MSP const *msp);

gboolean              mspHasSName(const MSP const *msp);
gboolean              mspHasSSeq(const MSP  const *msp);
gboolean              mspHasSCoords(const MSP const *msp);
gboolean              mspHasSStrand(const MSP const *msp);
gboolean              mspHasPolyATail(const MSP const *msp, const GList const *polyASiteList);
gboolean              mspCoordInPolyATail(const int coord, const MSP const *msp, const GList const *polyASiteList);

void                  writeBlxSequenceToOutput(FILE *pipe, const BlxSequence *blxSeq, IntRange *reqdRange);
BlxSequence*          readBlxSequenceFromText(char *text, int *numMsps);
void                  writeMspToOutput(FILE *pipe, const MSP const *msp);
void                  readMspFromText(MSP *msp, char *text);

MSP*                  createEmptyMsp(MSP **lastMsp, MSP **mspList);
MSP*                  createNewMsp(GList* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, const BlxMspType mspType, char *source, const gdouble score, const gdouble percentId, const int phase,
				   char *url, char *idTag, char *qName, const int qStart, const int qEnd, const BlxStrand qStrand, const int qFrame, 
				   char *sName, const int sStart, const int sEnd, const BlxStrand sStrand, char *sequence, 
				   GError **error);  

//void                  insertFS(MSP *msp, char *series);

/* Feature series */
gint		      fsSortByNameCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in);
gint		      fsSortByOrderCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in);

/* BlxSequence */
char*		      blxSequenceGetSummaryInfo(const BlxSequence const *blxSeq);
BlxSequence*          createEmptyBlxSequence(const char *fullName, const char *idTag, GError **error);
void                  addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error);
BlxSequence*          addBlxSequence(const char *name, const char *idTag, BlxStrand strand, GList **seqList, char *sequence, MSP *msp, GError **error);
BlxSequence*          findBlxSequence(GList *seqList, const char *reqdName, const char *reqdIdTag, const BlxStrand reqdStrand);
void		      blxSequenceSetName(BlxSequence *seq, const char *fullName);
const char*	      blxSequenceGetFullName(const BlxSequence *seq);
const char*	      blxSequenceGetVariantName(const BlxSequence *seq);
const char*	      blxSequenceGetDisplayName(const BlxSequence *seq);
const char*	      blxSequenceGetShortName(const BlxSequence *seq);
int		      blxSequenceGetLength(const BlxSequence *seq);
char*                 blxSequenceGetSeq(const BlxSequence *seq);
gboolean	      blxSequenceRequiresSeqData(const BlxSequence *seq);
gboolean	      blxSequenceRequiresOptionalData(const BlxSequence *seq);
BlxSequence*          blxSequenceGetVariantParent(const BlxSequence *variant, GList *allSeqs);
char*                 blxSequenceGetInfo(BlxSequence *blxSeq, const gboolean allowNewlines, const gboolean dataLoaded);
int		      blxSequenceGetStart(const BlxSequence *seq);
int		      blxSequenceGetEnd(const BlxSequence *seq);

void		      destroyBlxSequence(BlxSequence *seq);

#endif /* _blxmsp_included_ */
