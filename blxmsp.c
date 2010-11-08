/*
 *  blxmsp.c
 *  blixem
 *
 *  Created by Gemma Barson on 02/09/2010.
 *  Copyright 2010 Sanger Institute. All rights reserved.
 *
 */

#include "SeqTools/blxmsp.h"
#include "SeqTools/utilities.h"
#include <string.h>

static char*                    blxSequenceGetOrganism(const BlxSequence *seq);
static char*                    blxSequenceGetGeneName(const BlxSequence *seq);
static char*                    blxSequenceGetTissueType(const BlxSequence *seq);
static char*                    blxSequenceGetStrain(const BlxSequence *seq);


gboolean typeIsExon(const BlxMspType mspType)
{
  return (mspType == BLXMSP_CDS || mspType == BLXMSP_UTR || mspType == BLXMSP_EXON);
}

gboolean typeIsIntron(const BlxMspType mspType)
{
  return (mspType == BLXMSP_INTRON);
}

gboolean typeIsMatch(const BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH);
}

gboolean typeIsVariation(const BlxMspType mspType)
{
  return (mspType == BLXMSP_VARIATION);
}


gboolean mspIsExon(const MSP const *msp)
{
  return (msp && typeIsExon(msp->type));
}


/* Determine whether the given msp is in a visible layer */
gboolean mspLayerIsVisible(const MSP const *msp)
{
  gboolean result = TRUE;
  
  /* Currently only applicable to exons. Show plain exons OR their CDS/UTR sections,
   * but not both. The plan is to add some options to toggle layers on and off, but for 
   * now just hard code this. */
  if (msp->type == BLXMSP_EXON /* msp->type == BLXMSP_CDS || msp->type == BLXMSP_UTR */)
    {
    result = FALSE;
    }
  
  return result;
}


/* Determine whether the given MSP is in a coding region or untranslated region. For 
 * exons, this is determined by the exon type. For introns, we have to look at the
 * adjacent exons to determine whether to show them as CDS or UTR - only show it as
 * CDS if there is a CDS exon on both sides of the intron. */
static const GdkColor* mspGetIntronColor(const MSP const *msp, 
					 GArray *defaultColors,
					 const BlxSequence *blxSeq,
					 const gboolean selected,
					 const gboolean usePrintColors,
					 const gboolean fill,
					 const int exonFillColorId,
					 const int exonLineColorId,
					 const int cdsFillColorId,
					 const int cdsLineColorId,
					 const int utrFillColorId,
					 const int utrLineColorId)
{
  const GdkColor *result = NULL;
  
  /* Find the nearest exons before and after this MSP */
  const MSP *prevExon = NULL;
  const MSP *nextExon = NULL;
  
  GList *mspItem = blxSeq->mspList;
  for ( ; mspItem; mspItem = mspItem->next)
    {
    const MSP *curMsp = (const MSP *)(mspItem->data);
    
    if (mspIsExon(curMsp) && mspLayerIsVisible(curMsp))
      {
      const int curOffset = mspGetQStart(curMsp) - mspGetQStart(msp);
      
      if (curOffset < 0 && (!prevExon || curOffset > mspGetQStart(prevExon) - mspGetQStart(msp)))
	{
	/* Current MSP is before our MSP and is the smallest offset so far */
	prevExon = curMsp;
	}
      else if (curOffset > 0 && (!nextExon || curOffset < mspGetQStart(nextExon) - mspGetQStart(msp)))
	{
	/* Current MSP is after our MSP and is the smallest offset so far */
	nextExon = curMsp;
	}
      }
    }
  
  gboolean prevIsUtr = prevExon && prevExon->type == BLXMSP_UTR;
  gboolean nextIsUtr = nextExon && nextExon->type == BLXMSP_UTR;
  
  /* if either exon is UTR, the intron is UTR */
  if (prevIsUtr)
    {
    result = mspGetColor(prevExon, defaultColors, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (nextIsUtr)
    {
    result = mspGetColor(nextExon, defaultColors, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (prevExon)
    {
    /* Both exons (or the sole exon, if only one exists) are CDS, so the intron is CDS */
    result = mspGetColor(prevExon, defaultColors, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (nextExon)
    {
    /* This is the only exon and it is CDS, so make the intron CDS */
    result = mspGetColor(nextExon, defaultColors, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else
    {
    /* No exon exists adjacent to this intron: default to generic exon color for want of anything better to do. */
    if (fill)
      result = getGdkColor(exonFillColorId, defaultColors, selected, usePrintColors);
    else
      result = getGdkColor(exonLineColorId, defaultColors, selected, usePrintColors);
    }
  
  return result;
}

gboolean mspIsIntron(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_INTRON);
}

gboolean mspIsBlastMatch(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_MATCH);
}

gboolean mspIsPolyASite(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_POLYA_SITE);
}

gboolean mspIsVariation(const MSP const *msp)
{
  return (msp && typeIsVariation(msp->type));
}

gboolean mspIsZeroLenVariation(const MSP const *msp)
{
  return (mspIsVariation(msp) && msp->sSequence && msp->sSequence->sequence && msp->sSequence->sequence->str && msp->sSequence->sequence->str[0] == '-');
}

/* Whether the msp has a Target name (i.e. Subject sequence name) */
gboolean mspHasSName(const MSP const *msp)
{
  return TRUE;
}

/* Whether the MSP requires subject sequence coords to be set. Only matches 
 * and exons have coords on the subject sequence. (to do: is this optional for exons?) */
gboolean mspHasSCoords(const MSP const *msp)
{
  return mspIsExon(msp) || mspIsBlastMatch(msp);
}

/* Whether the MSP requires subject sequence strand to be set. Only matches 
 * require strand on the subject sequence, although exons may have them set. */
gboolean mspHasSStrand(const MSP const *msp)
{
  return mspIsBlastMatch(msp);
}

/* Whether the MSP requires the actual sequence data for the subject sequence. Only
 * matches require the sequence data. */
gboolean mspHasSSeq(const MSP  const *msp)
{
  return mspIsBlastMatch(msp);
}


/***********************************************************
 *		MSP data access functions		   * 
 ***********************************************************/

/* Get the range of coords of the alignment on the reference sequence */
const IntRange const *mspGetRefCoords(const MSP const *msp)
{
  return &msp->qRange;
}

/* Get the range of coords of the alignment on the match sequence */
const IntRange const *mspGetMatchCoords(const MSP const *msp)
{
  return &msp->sRange;
}

/* Return the length of the range of alignment coords on the ref seq */
int mspGetQRangeLen(const MSP const *msp)
{
  return msp->qRange.max - msp->qRange.min + 1;
}

/* Return the length of the range of alignment coords on the match seq */
int mspGetSRangeLen(const MSP const *msp)
{
  return msp->sRange.max - msp->sRange.min + 1;
}

/* Get the start (5 prime) coord of the alignment on the reference sequence. This is
 * the lowest value coord if the strand is forwards or the highest if it is reverse. */
int mspGetQStart(const MSP const *msp)
{
  return (mspGetRefStrand(msp) == BLXSTRAND_REVERSE ? msp->qRange.max : msp->qRange.min);
}

/* Get the end (3 prime) coord of the alignment on the reference sequence. This is
 * the highest value coord if the strand is forwards or the lowest if it is reverse. */
int mspGetQEnd(const MSP const *msp)
{
  return (mspGetRefStrand(msp) == BLXSTRAND_REVERSE ? msp->qRange.min : msp->qRange.max);
}

/* Get the start coord of the alignment on the match sequence. This is
 * the lowest value coord if the match strand is the same direction as the ref seq strand,
 * or the highest value coord otherwise. */
int mspGetSStart(const MSP const *msp)
{
  return (mspGetMatchStrand(msp) == mspGetRefStrand(msp) ? msp->sRange.min : msp->sRange.max);
}

/* Get the end coord of the alignment on the match sequence. This is
 * the highest value coord if the match strand is in the same direction as the ref seq strand, 
 * or the lowest value coord otherwise. */
int mspGetSEnd(const MSP const *msp)
{
  return (mspGetMatchStrand(msp) == mspGetRefStrand(msp) ? msp->sRange.max : msp->sRange.min);
}

/* Get the match sequence data for the given msp, or NULL if none exists */
char* mspGetSSeq(const MSP const *msp)
{
  return (msp && msp->sSequence ? blxSequenceGetSeq(msp->sSequence) : NULL);
}

/* Return the match sequence name. (Gets it from the BlxSequence if the MSP itself doesn't have
 * a name) */
const char *mspGetSName(const MSP const *msp)
{
  const char *result = NULL;
  
  if (msp)
    {
      if (msp->sname)
        {
          result = msp->sname;
        }
      else if (msp->sSequence)
        {
          result = blxSequenceGetFullName(msp->sSequence);
        }
    }
  
  return result;
}


/* Get the transcript name for an old-style exon/intron. These were postfixed with 'x' or 'i' 
 * to indicate exon and intron; we need to remove this postfix in order to find the real transcript
 * name so that we can group exons and introns from the same transcript together in the same BlxSequence.
 * If not an exon or intron, just returns a copy of the msp name. The result should be free'd with g_free. */
char* mspGetExonTranscriptName(const MSP *msp)
{
  char *name = g_strdup(msp->sname);
  
  if (name)
    {
    int i = strlen(name) - 1;
    
    if (mspIsExon(msp) && (name[i] == 'x' || name[i] == 'X'))
      {
      name[i] = '\0';
      }
    else if (mspIsIntron(msp) && (name[i] == 'i' || name[i] == 'I'))
      {
      name[i] = '\0';
      }
    }
  
  return name;
}

/* Return the length of the match sequence that the given MSP lies on */
int mspGetMatchSeqLen(const MSP const *msp)
{
  return blxSequenceGetLength(msp->sSequence);
}

/* Return the reading frame of the ref sequence that the given MSP is a match against */
int mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType)
{
  int result = UNSET_INT;
  
  if (seqType == BLXSEQ_DNA)
    {
    /* Ignore the frame in  the msp. For DNA matches we only have one frame on each strand. */
    result = 1;
    }
  else
    {
    result = msp->qFrame;
    }
  
  return result;
}

/* Return the strand of the ref sequence that the given MSP is a match against */
BlxStrand mspGetRefStrand(const MSP const *msp)
{
  return msp->qStrand;
}

/* Return the strand of the match sequence that the given MSP is a match on */
BlxStrand mspGetMatchStrand(const MSP const *msp)
{
  return msp->sSequence ? msp->sSequence->strand : BLXSTRAND_NONE;
}

/* Get the match sequence for the given MSP */
const char* mspGetMatchSeq(const MSP const *msp)
{
  return blxSequenceGetSeq(msp->sSequence);
}


/* Return the color matching the given properties from the given style */
static const GdkColor *styleGetColor(const BlxStyle *style, const gboolean selected, const gboolean usePrintColors, const gboolean fill, const gboolean utr)
{
  const GdkColor *result = NULL;
  
  if (fill)
    {
      if (utr)
        result = blxColorGetColor(&style->fillColorUtr, selected, usePrintColors);
      else
        result = blxColorGetColor(&style->fillColor, selected, usePrintColors);
    }
  else
    {
      if (utr)
        result = blxColorGetColor(&style->lineColorUtr, selected, usePrintColors);
      else
        result = blxColorGetColor(&style->lineColor, selected, usePrintColors);
    }
  return result;
}


/* Get the color for drawing the given MSP (If 'selected' is true, returns
 * the color when the MSP is selected.). Returns the fill color if 'fill' is 
 * true, otherwise the line color */
const GdkColor* mspGetColor(const MSP const *msp, 
			    GArray *defaultColors, 
			    const BlxSequence *blxSeq,
			    const gboolean selected, 
			    const gboolean usePrintColors, 
			    const gboolean fill,
			    const int exonFillColorId,
			    const int exonLineColorId,
			    const int cdsFillColorId,
			    const int cdsLineColorId,
			    const int utrFillColorId,
			    const int utrLineColorId)
{
  const GdkColor *result = NULL;
  
  if (msp->style)
    {
      result = styleGetColor(msp->style, selected, usePrintColors, fill, msp->type == BLXMSP_UTR);
    }
  
  if (!result)
    {
    /* Use the default color for this MSP's type */
    switch (msp->type)
      {
	case BLXMSP_EXON:
	result = getGdkColor(fill ? exonFillColorId : exonLineColorId, defaultColors, selected, usePrintColors);
	break;
	
	case BLXMSP_CDS:
	result = getGdkColor(fill ? cdsFillColorId : cdsLineColorId, defaultColors, selected, usePrintColors);
	break;
	
	case BLXMSP_UTR:
	result = getGdkColor(fill ? utrFillColorId : utrLineColorId, defaultColors, selected, usePrintColors);
	break;
	
	case BLXMSP_INTRON:
	result = mspGetIntronColor(msp, defaultColors, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
	break;
	
	default:
	break;
      };
    }
  
  return result;
}


/* Get functions from msps for various sequence properties. Result is owned by the sequence
 * and should not be free'd. Returns NULL if the property is not set. */
char *mspGetOrganism(const MSP const *msp)
{
  return (msp && msp->sSequence && msp->sSequence->organism ? msp->sSequence->organism->str : NULL);
}

char *mspGetGeneName(const MSP const *msp)
{
  return (msp && msp->sSequence && msp->sSequence->geneName ? msp->sSequence->geneName->str : NULL);
}

char *mspGetTissueType(const MSP const *msp)
{
  return (msp && msp->sSequence && msp->sSequence->tissueType ? msp->sSequence->tissueType->str : NULL);
}

char *mspGetStrain(const MSP const *msp)
{
  return (msp && msp->sSequence && msp->sSequence->strain ? msp->sSequence->strain->str : NULL);
}


/* Return the coords of an MSP as a string. The result should be free'd with g_free */
char *mspGetCoordsAsString(const MSP const *msp)
{
  char *result = NULL;
  
  if (msp)
    {
    GString *resultStr = g_string_new("");
    
    g_string_append_printf(resultStr, "%d - %d [%d - %d]", msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
    
    result = resultStr->str;
    g_string_free(resultStr, FALSE);
    }
  
  return result;
}


/* Returns true if a feature-series by the given name exists in the feature-series array and
 * and, if so, sets index_out with its index. */
static gboolean fsArrayFindByName(GArray *fsArray, FeatureSeries *fs, int *index_out)
{
  gboolean result = FALSE;
  
  if (fsArray)
    {
    int i = 0;
    for ( ; i < fsArray->len; ++i)
      {
      FeatureSeries *compareFs = &g_array_index(fsArray, FeatureSeries, i);
      
      if (!fsSortByNameCompareFunc(fs, compareFs))
	{
	result = TRUE;
	*index_out = i;
	break;
	}
      }
    }
  
  return result;
}


/* Comparison function to sort two Feature Series by the order number stored in the FeatureSeries
 * struct. Returns -1 if the first item is before the second, 1 if the second is first, or 0 if 
 * they are equal.  */
gint fsSortByOrderCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in)
{
  int result = 0;
  
  FeatureSeries *fs1 = (FeatureSeries *)fs1_in;
  FeatureSeries *fs2 = (FeatureSeries *)fs2_in;
  
  if (fs1->order < fs2->order)
    result = -1;
  else if (fs1->order > fs2->order)
    result = 1;
  
  return result;
}


/* Comparison function to sort two Features Series by name. */
gint fsSortByNameCompareFunc(gconstpointer fs1_in, gconstpointer fs2_in)
{
  FeatureSeries *fs1 = (FeatureSeries *)fs1_in;
  FeatureSeries *fs2 = (FeatureSeries *)fs2_in;
  
  /*printf("%s - %s : %d\n", fs1->name, fs2->name,  strcmp(fs1->name, fs2->name));*/
  
  return strcmp(fs1->name, fs2->name);
}


/* Insert the given MSP into the Feature Series of the given name. If the Feature Series
 * does not exist yet, create it. Also, if the Feature Series array does not exist yet, create
 * that too. */
void insertFS(MSP *msp, char *series)
{
  if (!fsArr) 
    {
    fsArr = g_array_sized_new(TRUE, FALSE, sizeof(FeatureSeries), 50);
    }
  
  static int orderNum = 0; /* will increment this each time we add a feature series to the array */
  
  FeatureSeries *fs = g_malloc(sizeof(FeatureSeries));
  fs->on = 1;
  fs->y = 0.0;
  fs->xy = (msp->type == BLXMSP_XY_PLOT ? 1 : 0);
  
  fs->name = g_strdup(series);
  
  int i;
  if (fsArrayFindByName(fsArr, fs, &i))
    {
    msp->fs = &g_array_index(fsArr, FeatureSeries, i);
    g_free(fs->name);
    g_free(fs);
    }
  else
    {
    /* Remember the order we added them so we can sort by it again later. */
    orderNum++;
    fs->order = orderNum;
    
    g_array_append_val(fsArr, *fs);
    //      g_array_sort(fsArr, fsSortByNameCompareFunc);
    
    msp->fs = fs;
    }
}



/* Returns true if there is a polyA site at the 3' end of this MSP's alignment range. The input
 * list should be a list containing all polya sites (and only polya sites) */
gboolean mspHasPolyATail(const MSP const *msp, const GList const *polyASiteList)
{
  gboolean found = FALSE;
  
  /* Only matches have polyA tails. */
  if (mspIsBlastMatch(msp))
    {
      /* For now, loop through all poly A sites and see if the site coord matches the 3' end coord of
       * the alignment. If speed proves to be an issue we could do some pre-processing to link MSPs 
       * to relevant polyA signals/sites so that we don't have to loop each time we want to check. */
      const GList *item = polyASiteList;
      
      for ( ; !found && item; item = item->next)
          {
            const MSP const *curPolyASite = (const MSP const*)(item->data);
            const int qEnd = mspGetQEnd(msp);
            
            if (mspGetRefStrand(msp) == BLXSTRAND_FORWARD)
              {
                found = (qEnd == curPolyASite->qRange.min);
              }
            else
              {
                found = (qEnd == curPolyASite->qRange.min + 1);
              }
        }
    }
  
  return found;
}


/* Returns true if the given MSP coord (in ref seq nucleotide coords) is inside a polyA tail, if
 * this MSP has one. */
gboolean mspCoordInPolyATail(const int coord, const MSP const *msp, const GList *polyASiteList)
{
  gboolean result = mspHasPolyATail(msp, polyASiteList);
  
  /* See if the coord is outside the 3' end of the alignment range (i.e. is greater than the
   * max coord if we're on the forward strand or less than the min coord if on the reverse). */
  result &= ((mspGetRefStrand(msp) == BLXSTRAND_FORWARD && coord > msp->qRange.max) ||
             (mspGetRefStrand(msp) == BLXSTRAND_REVERSE && coord < msp->qRange.min));
  
  return result;
}


/***********************************************************
 *		      BlxSequence			   * 
 ***********************************************************/

/* Append the contents of the second GString to the first GString, if the second is non-null,
 * including the given separator before the second GString is appended. Pass the separator
 * as an empty string if no separation is required. */
static void appendTextIfNonNull(GString *gstr, const char *separator, const GString *text)
{
  if (text && text->str)
    {
    g_assert(separator);
    g_string_append_printf(gstr, "%s%s", separator, text->str);
    }
}

/* Return summary info about a given BlxSequence (e.g. for displaying in the status bar). The
 * result should be free'd with g_free. */
char* blxSequenceGetSummaryInfo(const BlxSequence const *blxSeq)
{
  char *result = NULL;
  
  if (blxSeq)
    {
    GString *resultStr = g_string_new("");
    const char *separator = ";  ";
    
    g_string_append_printf(resultStr, "%s", blxSequenceGetFullName(blxSeq));
    appendTextIfNonNull(resultStr, separator, blxSeq->organism);
    appendTextIfNonNull(resultStr, separator, blxSeq->geneName);
    appendTextIfNonNull(resultStr, separator, blxSeq->tissueType);
    appendTextIfNonNull(resultStr, separator, blxSeq->strain);
    
    result = resultStr->str;
    g_string_free(resultStr, FALSE);
    }
  
  return result;
}


/* Returns the pointer to the BlxSequence that matches the given strand and sequence name/tag.
 * Returns NULL if no match was found. */
BlxSequence *findBlxSequence(GList *seqList, const char *reqdName, const char *reqdIdTag, const BlxStrand reqdStrand)
{
  BlxSequence *result = NULL;
  
  /* Loop through all sequences in the list */
  GList *listItem = seqList;
  
  for ( ; listItem; listItem = listItem->next)
    {
    BlxSequence *currentSeq = (BlxSequence*)(listItem->data);
    
    if (currentSeq->strand == reqdStrand &&
	( (reqdName && currentSeq->fullName && !strcmp(currentSeq->fullName, reqdName)) ||
	 (reqdIdTag && currentSeq->idTag && !strcmp(currentSeq->idTag, reqdIdTag)) ))
      {
      result = currentSeq;
      break;
      }
    }
  
  return result;
}



/* Return the full name of a BlxSequence (including prefix and variant) */
const char *blxSequenceGetFullName(const BlxSequence *seq)
{
  const char *result = NULL;
  
  if (seq)
    {
      if (seq->fullName)
        result = seq->fullName;
      else if (seq->idTag)
        result = seq->idTag;
      else
        g_error("Sequence does not have a name specified.\n");
    }
  
  return result;
}

/* Return the variant name of a BlxSequence (excludes prefix but includes variant) */
const char *blxSequenceGetVariantName(const BlxSequence *seq)
{
  return seq->variantName;
}

/* Return the display name of a BlxSequence (same as full name for now) */
const char *blxSequenceGetDisplayName(const BlxSequence *seq)
{
  return seq->fullName;
}

/* Return the short name of a BlxSequence (excludes prefix and variant number) */
const char *blxSequenceGetShortName(const BlxSequence *seq)
{
  return seq->shortName;
}

/* Return the length of the given blxsequence's sequence data */
int blxSequenceGetLength(const BlxSequence *seq)
{
  return (seq && seq->sequence ? seq->sequence->len : 0);
}

/* Get the start extent of the sequence on the ref sequence */
int blxSequenceGetStart(const BlxSequence *seq)
{
  return seq->qRange.min;
}

/* Get the end extend of the sequence on the ref sequence */
int blxSequenceGetEnd(const BlxSequence *seq)
{
  return seq->qRange.max;
}

/* Get the sequence data for the given blxsequence */
char *blxSequenceGetSeq(const BlxSequence *seq)
{
  return (seq && seq->sequence ? seq->sequence->str : NULL);
}

gboolean blxSequenceRequiresSeqData(const BlxSequence *seq)
{
  return (seq && (seq->type == BLXSEQUENCE_MATCH || seq->type == BLXSEQUENCE_VARIATION));
}

gboolean blxSequenceRequiresOptionalData(const BlxSequence *seq)
{
  return (seq && seq->type == BLXSEQUENCE_MATCH);
}

static char *blxSequenceGetOrganism(const BlxSequence *seq)
{
  return (seq && seq->organism ? seq->organism->str : "");
}

static char *blxSequenceGetGeneName(const BlxSequence *seq)
{
  return (seq && seq->geneName ? seq->geneName->str : "");
}

static char *blxSequenceGetTissueType(const BlxSequence *seq)
{
  return (seq && seq->tissueType ? seq->tissueType->str : "");
}

static char *blxSequenceGetStrain(const BlxSequence *seq)
{
  return (seq && seq->strain ? seq->strain->str : "");
}


/* Return all the stored info about a blx sequenece (description, organism, tissue type etc.) 
 * in a single string. The result should be free'd by the caller using g_free. If 'allowNewlines'
 * is true the data is separated with newline characters, otherwise with tabs (i.e. pass as false
 * if returned string is to be shown as a single line). The dataLoaded flag should be passed as
 * false if the optional data has not been loaded yet. */
char *blxSequenceGetInfo(BlxSequence *blxSeq, const gboolean allowNewlines, const gboolean dataLoaded)
{
  char *result = NULL;
  
  if (blxSeq)
    {
    GString *resultStr = g_string_new("");
    char separator = allowNewlines ? '\n' : '\t';
    char strand = blxSeq->strand == BLXSTRAND_REVERSE ? '-' : '+';
    char unloadedStr[] = "(optional data not loaded)";
    
    g_string_append_printf(resultStr, "SEQUENCE NAME:\t%s%c%c", blxSeq->fullName, strand, separator);
    g_string_append_printf(resultStr, "ORGANISM:\t\t\t%s%c", !dataLoaded ? unloadedStr : blxSequenceGetOrganism(blxSeq), separator);
    g_string_append_printf(resultStr, "GENE NAME:\t\t\t%s%c", !dataLoaded ? unloadedStr : blxSequenceGetGeneName(blxSeq), separator);
    g_string_append_printf(resultStr, "TISSUE TYPE:\t\t%s%c", !dataLoaded ? unloadedStr : blxSequenceGetTissueType(blxSeq), separator);
    g_string_append_printf(resultStr, "STRAIN:\t\t\t\t%s%c", !dataLoaded ? unloadedStr : blxSequenceGetStrain(blxSeq), separator);
    
    /* Loop through the MSPs and show their info */
    g_string_append(resultStr, "ALIGNMENTS:\t\t");
    GList *mspItem = blxSeq->mspList;
    
    for ( ; mspItem; mspItem = mspItem->next)
      {
      const MSP const *msp = (const MSP const *)(mspItem->data);
      g_string_append_printf(resultStr, "%s\t", mspGetCoordsAsString(msp));
      }
    
    /* Add a final separator after the alignments line*/
    g_string_append_printf(resultStr, "%c", separator);
    
    result = resultStr->str;
    g_string_free(resultStr, FALSE);
    }
  
  return result;
}

/* Find a BlxSequence by its full name (must be an exact match but is case insensitive) */
static BlxSequence* blxSequenceFindByName(const char *name, GList *allSeqs)
{
  BlxSequence *result = NULL;
  GList *listItem = allSeqs;
  
  for ( ; listItem; listItem = listItem->next)
    {
    BlxSequence *curSeq = (BlxSequence*)(listItem->data);
    
    if (stringsEqual(curSeq->fullName, name, FALSE))
      {
      result = curSeq;
      break;
      }
    }
  
  return result;
}


/* Get the "parent" sequence of the given protein variant. Assumes the variant
 * contains a dash '-' in the name followed by the variant number as a digit. 
 * e.g. SW:P51531-2.2. The function looks for a sequence with the same name but
 * with this dash and the following digit(s) (up to the end of the name or the '.'
 * if there is one) removed. Returns NULL if no parent was found. */
BlxSequence* blxSequenceGetVariantParent(const BlxSequence *variant, GList *allSeqs)
{
  BlxSequence *result = NULL;
  
  const char *variantName = blxSequenceGetFullName(variant);
  
  if (variantName)
    {
      char *parentName = g_strdup(variantName);
      char *insertPoint = strchr(parentName, '-');
      
      if (insertPoint)
        {
          /* Replace '-' by terminating char, in case there's nothing else to copy in. */
          *insertPoint = '\0';
          
          /* The insert point is where we'll copy into. Create another pointer that we'll increment
           * until we find a '.' and then we'll copy from that point. */
          char *copyPoint = insertPoint;
          ++copyPoint;
          
          gboolean foundRestartPoint = FALSE; /* set to true when we find where to start copying from again */
          
          while (copyPoint && *copyPoint != '\0')
              {
              if (foundRestartPoint)
                {
                  *insertPoint = *copyPoint;
                  ++insertPoint;
                }
              else if (*copyPoint == '.')
                {
                  foundRestartPoint = TRUE;
                  *insertPoint = *copyPoint;
                  ++insertPoint;
                }
          
              ++copyPoint;
            }
        
          *insertPoint = '\0';
        
          result = blxSequenceFindByName(parentName, allSeqs);
          g_free(parentName);
        }
    }
  
  return result;
}


/* Frees all memory used by a BlxSequence */
void destroyBlxSequence(BlxSequence *seq)
{
  if (seq)
    {
    g_free(seq->fullName);
    g_free(seq->variantName);
    g_free(seq->shortName);
    
    if (seq->sequence)      g_string_free(seq->sequence, TRUE);
    if (seq->organism)      g_string_free(seq->organism, TRUE);
    if (seq->geneName)      g_string_free(seq->geneName, TRUE);
    if (seq->tissueType)    g_string_free(seq->tissueType, TRUE);
    if (seq->strain)        g_string_free(seq->strain, TRUE);
    
    g_free(seq);
    }
}


/* Set the name of a BlxSequence (also sets the short name etc. Does nothing if the 
 * name is already set.) */
void blxSequenceSetName(BlxSequence *seq, const char *fullName)
{  
  if (fullName && !seq->fullName)
    {
      seq->fullName = fullName ? g_strdup(fullName) : NULL;
      
      /* The variant name: just cut off the prefix chars. We can use a pointer into
       * the original string. */
      seq->variantName = g_strdup(getSeqVariantName(fullName));
      
      /* The short name: cut off the prefix chars (before the ':') and the variant
       * number (after the '.'). Need to duplicate the string to change the end of it. */
      seq->shortName = g_strdup(seq->variantName);
      char *cutPoint = strchr(seq->shortName, '.');
      
      if (cutPoint)
        {
          *cutPoint = '\0';
        }
    }
}


/* Utility to create a BlxSequence with the given name. */
BlxSequence* createEmptyBlxSequence(const char *fullName, const char *idTag, GError **error)
{
  BlxSequence *seq = g_malloc(sizeof(BlxSequence));
  
  seq->type = BLXSEQUENCE_UNSET;
  seq->idTag = idTag ? g_strdup(idTag) : NULL;
  
  seq->fullName = NULL;
  seq->variantName = NULL;
  seq->shortName = NULL;
  blxSequenceSetName(seq, fullName);
  
  seq->mspList = NULL;
  seq->sequence = NULL;
  seq->sequenceReqd = FALSE;
  
  seq->organism = NULL;
  seq->geneName = NULL;
  seq->tissueType = NULL;
  seq->strain = NULL;
  
  return seq;
}

/* returns true if the given msp should be output */
static gboolean outputMsp(const MSP const *msp, IntRange *validRange)
{
  return ((msp->type == BLXMSP_FS_SEG || mspIsExon(msp) || mspIsIntron(msp) || mspIsBlastMatch(msp))
           && rangesOverlap(&msp->qRange, validRange));
}

static int countMspsToOutput(const BlxSequence const *blxSeq, IntRange *validRange)
{
  int numMsps = 0;
  
  GList *mspItem = blxSeq->mspList;
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *msp = (const MSP const *)(mspItem->data);
      
      if (outputMsp(msp, validRange))
        ++numMsps;
    }
  
  return numMsps;
}


/* write data from the given blxsequence to the given output pipe. if validRange is given, only 
 * blxsequences that overlap that range will be output */
void writeBlxSequenceToOutput(FILE *pipe, const BlxSequence *blxSeq, IntRange *validRange)
{
  gboolean outputSeq = (blxSeq && (blxSeq->type == BLXSEQUENCE_TRANSCRIPT || blxSeq->type == BLXSEQUENCE_MATCH));
  
  int numMsps = countMspsToOutput(blxSeq, validRange);
  outputSeq &= numMsps > 0; /* only output the sequence if it has some valid msps */
  
  if (outputSeq)
    {
      fprintf(pipe, "%d %d %d",
              blxSeq->type,
              blxSeq->strand,
              numMsps); /* output number of msps so we know how many to read in */
      
      stringProtect(pipe, blxSeq->fullName);
      stringProtect(pipe, blxSeq->idTag);
      
      fputc('\n', pipe);
      
      /* now output the msps */
      GList *mspItem = blxSeq->mspList;
      for ( ; mspItem; mspItem = mspItem->next)
        {
          const MSP const *msp = (const MSP const *)(mspItem->data);
          
          if (outputMsp(msp, validRange))
            {
              writeMspToOutput(pipe, msp);
            }
        }
    }
}


/* increment the given char pointer until it is not whitespace */
static void nextChar(char **curChar)
{
  while (**curChar == ' ')
    ++(*curChar);
}


/* read in a blxsequence from the given text. */
BlxSequence* readBlxSequenceFromText(char *text, int *numMsps)
{
  DEBUG_ENTER("readBlxSequenceFromText(text=%s)", text);

  char *curChar = text;
  
  GError *error = NULL;
  BlxSequence *blxSeq = createEmptyBlxSequence(NULL, BLXSEQUENCE_UNSET, &error);
  
  if (error)
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_CRITICAL);
      DEBUG_EXIT("readBlxSequenceFromText returning NULL");
      return NULL;
    }
  
  nextChar(&curChar);
  
  blxSeq->type = strtol(curChar, &curChar, 10);
  nextChar(&curChar);

  blxSeq->strand = strtol(curChar, &curChar, 10);
  nextChar(&curChar);

  *numMsps = strtol(curChar, &curChar, 10);
  nextChar(&curChar);
  
  blxSeq->fullName = stringUnprotect(&curChar, NULL);
  blxSeq->idTag = stringUnprotect(&curChar, NULL);
  
  DEBUG_EXIT("readBlxSequenceFromText returning numMsps=%d", *numMsps);
  return blxSeq;
}


/* write data from the given msp to the given output pipe. */
void writeMspToOutput(FILE *pipe, const MSP const *msp)
{
  fprintf(pipe, "%d %f %d %d %d %d %d %d %d %d", 
          msp->type,
          msp->score, 
          msp->phase,
          msp->fsColor, 
          msp->qRange.min,
          msp->qRange.max,
          msp->sRange.min,
          msp->sRange.max,
          msp->qStrand,
          msp->qFrame);
  stringProtect(pipe, msp->qname);
  stringProtect(pipe, msp->sname);
  stringProtect(pipe, msp->desc);
  fputc('\n', pipe);
}


/* Read in an msp from the given text. the text is in the format as written by writeMspToFile  */
void readMspFromText(MSP *msp, char *text)
{
  DEBUG_ENTER("readMspFromText(text=%s)", text);

  char *curChar = text;

  nextChar(&curChar);
  msp->type = strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->score = g_strtod(curChar, &curChar);

  nextChar(&curChar);
  msp->phase = strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->fsColor = strtol(curChar, &curChar, 10);

  /* ref seq range */
  nextChar(&curChar);
  const int qStart = strtol(curChar, &curChar, 10); 
  nextChar(&curChar);
  const int qEnd = strtol(curChar, &curChar, 10); 
  intrangeSetValues(&msp->qRange, qStart, qEnd);

  /* match seq range */
  nextChar(&curChar);
  const int sStart = strtol(curChar, &curChar, 10); 
  nextChar(&curChar);
  const int sEnd = strtol(curChar, &curChar, 10); 
  intrangeSetValues(&msp->sRange, sStart, sEnd);
  
  nextChar(&curChar);
  msp->qStrand = strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->qFrame = strtol(curChar, &curChar, 10);

  nextChar(&curChar);
  msp->qname = stringUnprotect(&curChar, NULL);
  msp->sname = stringUnprotect(&curChar, NULL);
  msp->desc = stringUnprotect(&curChar, NULL);

  DEBUG_EXIT("readMspFromText returning");
}


/* Allocate memory for an MSP and initialise all its fields to relevant 'empty' values.
 * Add it into the given MSP list and make the end pointer ('lastMsp') point to the new
 * end of the list. NB the msplist is obsolete and is being replaced by featureLists, but we don't
 * add it to a feature list yet because we don't know its type. Returns a pointer to the newly-created MSP */
MSP* createEmptyMsp(MSP **lastMsp, MSP **mspList)
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  msp->next = NULL;
  msp->childMsps = NULL;
  msp->type = BLXMSP_INVALID;
  msp->score = 0.0;
  msp->id = 0.0;
  msp->phase = 0;
  msp->url = NULL;
  
  msp->qname = NULL;
  msp->qFrame = UNSET_INT;
  
  msp->qRange.min = 0;
  msp->qRange.max = 0;
  
  msp->sSequence = NULL;
  msp->sname = NULL;
  
  msp->sRange.min = 0;
  msp->sRange.max = 0;
  
  msp->desc = NULL;
  msp->source = NULL;
  
  msp->style = NULL;
  
  msp->fs = NULL;
  msp->fsColor = 0;
  msp->fsShape = BLXCURVE_BADSHAPE;
  
  msp->xy = NULL;
  msp->gaps = NULL;
  
#ifdef ACEDB
  msp->key = 0;
#endif
  
  /* Add it to the list */
  if (!*mspList) 
    {
      /* Nothing in the list yet: make this the first entry */
      *mspList = msp;
    }
  
  if (*lastMsp)
    {
      /* Tag it on to the end of the list */
      (*lastMsp)->next = msp;
    }
  
  /* Make the 'lastMsp' pointer point to the new end of the list */
  *lastMsp = msp;
  
  return msp;
}
