/*  File: blxmsp.h
 *  Author: Gemma Barson, 2010-09-02
 *  Copyright (c) 2010 - 2012 Genome Research Ltd
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
 *      Ed  fGriffiths      (Sanger Institute, UK)  <edgrif@sanger.ac.uk>
 *      Roy Storey        (Sanger Institute, UK)  <rds@sanger.ac.uk>
 *      Malcolm Hinsley   (Sanger Institute, UK)  <mh17@sanger.ac.uk>
 *
 * Description: Defines the MSP data struct and related functions. I think 
 *              MSP stood for Matching Segment Pair and originally represented
 *              an alignment. However, it is now used for other feature
 *              types too, so should be renamed to "BlxFeature" or something.
 *              Ideally also it would be separated out into a base feature and
 *              derived features for the different features types.
 *----------------------------------------------------------------------------
 */

#ifndef _blxmsp_included_
#define _blxmsp_included_

#include <seqtoolsUtils/utilities.h>
#include <stdlib.h>
#include <gtk/gtk.h>


//extern GArray    *fsArr;                 /* in dotter.c */

/* These are used by both blixem and dotter */
#define selectFeaturesStr          "Feature series selection tool"
#define XY_NOT_FILLED -1000        /* Magic value meaning "value not provided" */


/* Key names for values in data-type stanzas in the config file */
#define SEQTOOLS_BULK_FETCH          "bulk-fetch"
#define SEQTOOLS_USER_FETCH          "user-fetch"


/* Main Blixem error domain */
#define BLX_ERROR g_quark_from_string("Blixem")

/* Error codes */
typedef enum
{
  BLX_ERROR_SEQ_SEGMENT,              /* error finding sequence segment */
  BLX_ERROR_EMPTY_STRING,           /* error code for when user entered a zero-length string */
  BLX_ERROR_STRING_NOT_FOUND,       /* error code for when a search string is not found */
  BLX_ERROR_SEQ_NAME_NOT_FOUND,     /* the sequence name(s) being searched for were not found */
  BLX_ERROR_SEQ_DATA_MISMATCH,      /* same sequence was parsed more than once and data does not match */
  BLX_ERROR_INVALID_COLUMN          /* error when an invalid column is requested */
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
  BLXMSP_EXON,                   /* Exon (should appear AFTER CDS and UTR for sorting, as required by constructTranscriptData) */
  BLXMSP_POLYA_SITE,             /* polyA tail site */
  BLXMSP_POLYA_SIGNAL,           /* polyA signal */
  
  BLXMSP_VARIATION,              /* SNP, substitution, deletion, insertion */
  BLXMSP_SHORT_READ,             /* one fragment of a read-pair */
  
  BLXMSP_HSP,                    /* obsolete? */
  BLXMSP_GSP,                    /* obsolete? */
  
  BLXMSP_FS_SEG,                 /* Feature Series Segment - obsolete? */
  BLXMSP_XY_PLOT,                /* x/y coordinates - for plotting feature-series curves - obsolete? */

  BLXMSP_REGION,                 /* Region */
  BLXMSP_GAP,                    /* Gap, e.g. assembly gap */

  
  
  BLXMSP_NUM_TYPES,               /* the number of valid MSP types - any types following this may be used
                                   * e.g. for parsing, but no real MSP will be created from them */
  
  BLXMSP_TRANSCRIPT              /* Transcript */
} BlxMspType;


/* Type definition for BlxSequences */
typedef enum
  {
    BLXSEQUENCE_UNSET,
    BLXSEQUENCE_TRANSCRIPT,         /* transcript (i.e. collection of exons and introns) */
    BLXSEQUENCE_MATCH,              /* match sequence (i.e. collection of matches) */
    BLXSEQUENCE_VARIATION,          /* variation (i.e. insertion, deletion or substitution) */
    BLXSEQUENCE_READ_PAIR,          /* read pair */
    BLXSEQUENCE_REGION              /* region */
  } BlxSequenceType;


/* This enum provides an ID for each type of data model that we use for the
 * detail-view trees. */
typedef enum
  {
    BLXMODEL_NORMAL,                /* the normal model, where each row contains one feature */
    BLXMODEL_SQUASHED,              /* the "squashed" model, where all MSPs from the same sequence appear on the same row */
    
    BLXMODEL_NUM_MODELS             /* the number of model IDs. MUST BE LAST IN LIST */
  } BlxModelId;



/* Defines a data type for sequences. The data type contains properties applicable
 * to multiple sequences, e.g. which fetch method to use. */
typedef struct _BlxDataType
  {
    GQuark name;           /* the name of the data-type */
    GArray *bulkFetch;     /* list of fetch methods (by name as a GQuark) to use when bulk fetching sequences, in order of priority */
    GArray *userFetch;     /* list of fetch methods (by name as a GQuark) to use when user fetches a sequence, in order of priority */
  } BlxDataType;


/* Structure that contains information about a sequence */
typedef struct _BlxSequence
{
  BlxSequenceType type;            /* What type of collection of MSPs this is */
  BlxDataType *dataType;           /* Optional data type that specifies additional properties for this type of sequence data */

  char *idTag;                     /* Unique identifier e.g. from ID tag in GFF files */
  char *source;                    /* Optional source text for the sequence */

  GQuark fullName;                 /* full name of the sequence, including variant postfix, e.g. AV274505.2 */
  char *shortName;                 /* short name of the sequence, excluding variant, e.g. AV274505 */
  
  GString *organism;               /* organism from the EMBL data OS line */
  GString *geneName;               /* gene name from the EMBL data GN line */
  GString *tissueType;             /* tissue type from the /tissue_type attribute from the FT lines in the EMBL file */
  GString *strain;                 /* strain from the /strain attribute from the FT lines in the EMBL file */
  
  BlxStrand strand;                /* which strand of the sequence this is */
  GString *sequence;               /* the actual sequence data */
  gboolean sequenceReqd;           /* whether the sequence data is required (e.g. it is not needed for exons/introns etc.) */
  IntRange qRangeFwd;              /* the extent of alignments from this sequence on the ref sequence forward strand */ 
  IntRange qRangeRev;              /* the extent of alignments from this sequence on the ref sequence reverse strand */ 
  
  GList *mspList;                  /* list of MSPs from this sequence */
} BlxSequence;


typedef struct _FeatureSeries
{
  char        *name;             /* Name of the Feature Series */
  int         on;                /* Flag to show/hide the Feature Series */
  int         order;             /* Order number used for sorting */
  float       x;                   /* Series offset on x axis, to bump series on the screen */
  float       y;                   /* Series offset on y axis */
  int         xy;                  /* Flag for XY plot series */
} FeatureSeries;


/* Shapes of XY curves */
typedef enum
{ 
  BLXCURVE_PARTIAL, 
  BLXCURVE_INTERPOLATE, 
  BLXCURVE_BADSHAPE
} BlxCurveShape;


/* Structure holding information about a feature (see note at the top of this
 * file about the naming of this struct). */
typedef struct _MSP
{
  gchar* treePaths[BLXMODEL_NUM_MODELS]; /* identifies the row in the tree data model this msp is in (for the model given by the modelId) */

  struct _MSP       *next;
  GList             *childMsps;    /* Child MSPs of this MSP if it has them, e.g. an exon has CDS and UTR children (part_of relationship). */
  
  BlxMspType        type;          /* The type of the MSP, e.g. match, exon, SNP etc. */
  gdouble           score;         /* Score as a percentage. Technically this should be a weighted score taking into account gaps, length of the match etc., but for unknown reasons the ID has always been passed instead of score and the ID gets stored in here */
  gdouble           id;            /* Identity as a percentage. A simple comparison of bases within the match, ignoring gaps etc. Currently this is calculated internally by blixem. */
  int               phase;         /* phase: q start coord is offset by this amount to give the first base in the first complete codon (only relevant to CDSs) */
  char              *url;          /* URL to info about the MSP (e.g. variations have a URL which can be opened by double-clicking the variation) */
  
  char              *qname;        /* For Dotter, the MSP can belong to either sequence */
  IntRange          qRange;        /* the range of coords on the ref sequence where the alignment lies */
  BlxStrand         qStrand;       /* which strand on the reference sequence the match is on */
  int               qFrame;        /* which frame on the reference sequence the match is on */
  
  BlxSequence       *sSequence;    /* pointer to a struct holding info about the sequence/strand this match is from */
  char              *sname;        /* sequence name (could be different to the sequence name in the blxSequence e.g. exons have a postfixed 'x') */
  IntRange          sRange;        /* the range of coords on the match sequence where the alignment lies */
  
  /* The following ranges are all calculated from the above but are
   * cached in the MSP because they are used a lot. Note that these
   * these are not current used by dotter so dotter does not bother to
   * initialise them. */
  IntRange          displayRange;  /* the same range as qRange but in display coords */
  IntRange          fullRange;     /* the full range of display coords to show this match against (includes any unaligned portions of sequence that we're showing) */
  IntRange          fullSRange;    /* the full range of coords on the match sequence that we're showing (including any unaligned portions of sequence) */
  
  char              *desc;         /* Optional description text for the MSP */
  GSList            *gaps;         /* Array of "gaps" in this homolgy (this is a bit of a misnomer because the array
                                    * gives the ranges of the bits that align rather than the ranges of the gaps in between */
  
  BlxStyle          *style;        /* Specifies drawing style for this MSP, e.g. fill color and line color */
  
  /* obsolete? */
  FeatureSeries     *fs;           /* Feature series that this MSP belongs to */
  int               fsColor;       /* Color to draw this MSP in the feature series */
  BlxCurveShape     fsShape;       /* Shape data for drawing feature series curves, i.e. XY type PARTIAL or INTERPOLATE shapes */
  GArray            *xy;            /* For XY plot feature series */
} MSP ;


/* MSP functions */
gboolean              typeIsExon(const BlxMspType mspType);
gboolean              typeIsIntron(const BlxMspType mspType);
gboolean              typeIsTranscript(const BlxMspType mspType);
gboolean              typeIsMatch(const BlxMspType mspType);
gboolean              typeIsVariation(const BlxMspType mspType);
gboolean              typeIsShortRead(const BlxMspType mspType);
gboolean              typeIsRegion(const BlxMspType mspType);
gboolean              typeShownInDetailView(const BlxMspType mspType);
gboolean              blxSequenceShownInDetailView(const BlxSequence *blxSeq);
gboolean              blxSequenceShownInGrid(const BlxSequence *blxSeq);

const char*           mspGetRefName(const MSP const *msp);
int                   mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType);
BlxStrand             mspGetRefStrand(const MSP const *msp);
BlxStrand             mspGetMatchStrand(const MSP const *msp);
const char*           mspGetMatchSeq(const MSP const *msp);
const char*           mspGetSource(const MSP const *msp);
const char*           mspGetSName(const MSP *msp);
const IntRange const* mspGetRefCoords(const MSP const *msp);
const IntRange const* mspGetMatchCoords(const MSP const *msp);
int                   mspGetQStart(const MSP const *msp);
int                   mspGetQEnd(const MSP const *msp);
int                   mspGetSStart(const MSP const *msp);
int                   mspGetSEnd(const MSP const *msp);
char*                 mspGetSSeq(const MSP const *msp);
int                   mspGetQRangeLen(const MSP const *msp);
int                   mspGetSRangeLen(const MSP const *msp);
int                   mspGetMatchSeqLen(const MSP const *msp);

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
gchar*                mspGetTreePath(const MSP const *msp, BlxModelId modelId);

MSP*                  mspArrayIdx(const GArray const *array, const int idx);
gint                  compareFuncMspPos(gconstpointer a, gconstpointer b);
gint                  compareFuncMspArray(gconstpointer a, gconstpointer b);

gboolean              mspLayerIsVisible(const MSP const *msp);
gboolean              mspIsExon(const MSP const *msp);
gboolean              mspIsIntron(const MSP const *msp);
gboolean              mspIsSnp(const MSP const *msp);
gboolean              mspIsBlastMatch(const MSP const *msp);
gboolean              mspIsPolyASite(const MSP const *msp);
gboolean              mspIsVariation(const MSP const *msp);
gboolean              mspIsShortRead(const MSP const *msp);
gboolean              mspIsZeroLenVariation(const MSP const *msp);

gboolean              mspHasSName(const MSP const *msp);
gboolean              mspHasSSeq(const MSP  const *msp);
gboolean              mspHasSCoords(const MSP const *msp);
gboolean              mspHasSStrand(const MSP const *msp);
gboolean              mspHasPolyATail(const MSP const *msp, const GArray const *polyASiteList);
gboolean              mspCoordInPolyATail(const int coord, const MSP const *msp, const GArray const *polyASiteList);

int                   getMaxMspLen();
void                  setMaxMspLen(const int len);

void                  writeBlxSequenceToOutput(FILE *pipe, const BlxSequence *blxSeq, IntRange *range1, IntRange *range2);
BlxSequence*          readBlxSequenceFromText(char *text, int *numMsps);
void                  writeMspToOutput(FILE *pipe, const MSP const *msp);
void                  readMspFromText(MSP *msp, char *text);

void                  destroyMspList(MSP **mspList);
void                  destroyBlxSequenceList(GList **seqList);
void                  destroyMspData(MSP *msp);
MSP*                  createEmptyMsp(MSP **lastMsp, MSP **mspList);
MSP*                  createNewMsp(GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, const BlxMspType mspType, 
                                   BlxDataType *dataType, const char *source, const gdouble score, const gdouble percentId, const int phase,
                                   const char *url, const char *idTag, const char *qName, const int qStart, const int qEnd, 
                                   const BlxStrand qStrand, const int qFrame, const char *sName, const int sStart, const int sEnd, 
                                   const BlxStrand sStrand, char *sequence, GError **error);  
MSP*                  copyMsp(const MSP const *src, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, GError **error);

//void                  insertFS(MSP *msp, char *series);

void                  finaliseBlxSequences(GArray* featureLists[], MSP **mspList, GList **seqList, const int offset, const BlxSeqType seqType, const int numFrames, const IntRange const *refSeqRange, const gboolean calcFrame);
int                   findMspListSExtent(GList *mspList, const gboolean findMin);
int                   findMspListQExtent(GList *mspList, const gboolean findMin, const BlxStrand strand);

/* Feature series */
gint                  fsSortByNameCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in);
gint                  fsSortByOrderCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in);

/* BlxSequence */
char*                 blxSequenceGetSummaryInfo(const BlxSequence const *blxSeq);
BlxSequence*          createEmptyBlxSequence(const char *fullName, const char *idTag, GError **error);
BlxDataType*          createBlxDataType();
void                  destroyBlxDataType(BlxDataType **blxDataType);
const char*           getDataTypeName(BlxDataType *blxDataType);
void                  addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error);
BlxSequence*          addBlxSequence(const char *name, const char *idTag, BlxStrand strand, BlxDataType *dataType, const char *source, GList **seqList, char *sequence, MSP *msp, GError **error);
void                  blxSequenceSetName(BlxSequence *seq, const char *fullName);
const char*           blxSequenceGetFullName(const BlxSequence *seq);
const char*           blxSequenceGetDisplayName(const BlxSequence *seq);
GQuark                blxSequenceGetFetchMethod(const BlxSequence *seq, const gboolean bulk, const int index, const GArray *defaultMethods);
const char*           blxSequenceGetShortName(const BlxSequence *seq);
int                   blxSequenceGetLength(const BlxSequence *seq);
char*                 blxSequenceGetSeq(const BlxSequence *seq);
gboolean              blxSequenceRequiresSeqData(const BlxSequence *seq);
gboolean              blxSequenceRequiresOptionalData(const BlxSequence *seq);
BlxSequence*          blxSequenceGetVariantParent(const BlxSequence *variant, GList *allSeqs);
char*                 blxSequenceGetInfo(BlxSequence *blxSeq, const gboolean allowNewlines, const gboolean dataLoaded);
int                   blxSequenceGetStart(const BlxSequence *seq, const BlxStrand strand);
int                   blxSequenceGetEnd(const BlxSequence *seq, const BlxStrand strand);
const char*           blxSequenceGetSource(const BlxSequence *seq);
char*                 blxSequenceGetOrganism(const BlxSequence *seq);
char*                 blxSequenceGetGeneName(const BlxSequence *seq);
char*                 blxSequenceGetTissueType(const BlxSequence *seq);
char*                 blxSequenceGetStrain(const BlxSequence *seq);
char*                 blxSequenceGetFasta(const BlxSequence *seq);

void                  destroyBlxSequence(BlxSequence *seq);

#endif /* _blxmsp_included_ */
