/*  File: blxmsp.c
 *  Author: Gemma Barson, 2010-09-02
 *  Copyright (c) 2010 Genome Research Ltd
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
 * Description: See blxmsp.h
 *----------------------------------------------------------------------------
 */

#include <seqtoolsUtils/blxmsp.h>
#include <seqtoolsUtils/utilities.h>
#include <string.h>

#define ALL_READ_PAIRS_FWD      FALSE


/* Globals */
static int g_MaxMspLen = 0;     /* max length in display coords of all MSPs in the detail-view */


/* Local function declarations */
static char*                    blxSequenceGetOrganism(const BlxSequence *seq);
static char*                    blxSequenceGetGeneName(const BlxSequence *seq);
static char*                    blxSequenceGetTissueType(const BlxSequence *seq);
static char*                    blxSequenceGetStrain(const BlxSequence *seq);
static const char*              blxSequenceGetSource(const BlxSequence *seq);



/* Get/set the max MSP length */
int getMaxMspLen()
{
  return g_MaxMspLen;
}

void setMaxMspLen(const int len)
{
  g_MaxMspLen = len;
}


/* Type determination methods */
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

gboolean typeIsShortRead(const BlxMspType mspType)
{
  return (mspType == BLXMSP_SHORT_READ);
}

/* This returns true if the given type should be shown in the detail-view */
gboolean typeShownInDetailView(const BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH || mspType == BLXMSP_SHORT_READ || mspType == BLXMSP_CDS || mspType == BLXMSP_UTR);
}

/* This returns true if the given sequence should be shown in the detail-view */
gboolean blxSequenceShownInDetailView(const BlxSequence *blxSeq)
{
  return (blxSeq->type == BLXSEQUENCE_MATCH || blxSeq->type == BLXSEQUENCE_READ_PAIR || blxSeq->type == BLXSEQUENCE_TRANSCRIPT);
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
  return (msp && (msp->type == BLXMSP_MATCH || msp->type == BLXMSP_SHORT_READ));
}

gboolean mspIsPolyASite(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_POLYA_SITE);
}

gboolean mspIsVariation(const MSP const *msp)
{
  return (msp && typeIsVariation(msp->type));
}

gboolean mspIsShortRead(const MSP const *msp)
{
  return (msp && typeIsShortRead(msp->type));
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
  BlxStrand result = msp->qStrand;
  
  if (ALL_READ_PAIRS_FWD && mspIsShortRead(msp))
    result = BLXSTRAND_FORWARD;
  
  return result;
}

/* Return the strand of the match sequence that the given MSP is a match on */
BlxStrand mspGetMatchStrand(const MSP const *msp)
{
  BlxStrand result = (msp->sSequence ? msp->sSequence->strand : BLXSTRAND_NONE);

  /* If we're displaying all read-pairs on the forward strand, return forward */
  if (ALL_READ_PAIRS_FWD && mspIsShortRead(msp))
    result = BLXSTRAND_FORWARD;
  
  return result;
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
	
        /* to do: mspGetIntronColor() is non-trivial, because it has to work out the color
         * from the adjacent exons. Since mspGetColor() is called many times on re-draw, it
         * would be better to work out whether an intron is CDS or UTR during initialisation 
         * and use different types (e.g. BLXMSP_INTRON_CDS) that can be queried here to quickly
         * determine what color to use. */
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


/* Return the path in the given tree model that this MSP lies in */
gchar* mspGetTreePath(const MSP const *msp, BlxModelId modelId)
{
  return msp->treePaths[modelId];
}


///* Returns true if a feature-series by the given name exists in the feature-series array and
// * and, if so, sets index_out with its index. */
//static gboolean fsArrayFindByName(GArray *fsArray, FeatureSeries *fs, int *index_out)
//{
//  gboolean result = FALSE;
//  
//  if (fsArray)
//    {
//    int i = 0;
//    for ( ; i < fsArray->len; ++i)
//      {
//      FeatureSeries *compareFs = &g_array_index(fsArray, FeatureSeries, i);
//      
//      if (!fsSortByNameCompareFunc(fs, compareFs))
//	{
//	result = TRUE;
//	*index_out = i;
//	break;
//	}
//      }
//    }
//  
//  return result;
//}


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
//void insertFS(MSP *msp, char *series)
//{
//  if (!fsArr) 
//    {
//    fsArr = g_array_sized_new(TRUE, FALSE, sizeof(FeatureSeries), 50);
//    }
//  
//  static int orderNum = 0; /* will increment this each time we add a feature series to the array */
//  
//  FeatureSeries *fs = g_malloc(sizeof(FeatureSeries));
//  fs->on = 1;
//  fs->y = 0.0;
//  fs->xy = (msp->type == BLXMSP_XY_PLOT ? 1 : 0);
//  
//  fs->name = g_strdup(series);
//  
//  int i;
//  if (fsArrayFindByName(fsArr, fs, &i))
//    {
//    msp->fs = &g_array_index(fsArr, FeatureSeries, i);
//    g_free(fs->name);
//    g_free(fs);
//    }
//  else
//    {
//    /* Remember the order we added them so we can sort by it again later. */
//    orderNum++;
//    fs->order = orderNum;
//    
//    g_array_append_val(fsArr, *fs);
//    //      g_array_sort(fsArr, fsSortByNameCompareFunc);
//    
//    msp->fs = fs;
//    }
//}
//


/* Get the MSP at the given index in the given array. Returns null if out of bounds */
MSP* mspArrayIdx(const GArray const *array, const int idx)
{
  MSP *msp = NULL;
  
  if (idx >= 0 && idx < array->len)
    msp = g_array_index(array, MSP*, idx);
  
  return msp;
}


/* Returns true if there is a polyA site at the 3' end of this MSP's alignment range. The input
 * list should be a list containing all polya sites (and only polya sites) */
gboolean mspHasPolyATail(const MSP const *msp, const GArray const *polyASiteList)
{
  gboolean found = FALSE;
  
  /* Only matches have polyA tails. */
  if (mspIsBlastMatch(msp))
    {
      /* For now, loop through all poly A sites and see if the site coord matches the 3' end coord of
       * the alignment. If speed proves to be an issue we could do some pre-processing to link MSPs 
       * to relevant polyA signals/sites so that we don't have to loop each time we want to check. */
      int i = 0;
      MSP *curPolyASite = mspArrayIdx(polyASiteList, i);
      
      for ( ; !found && curPolyASite; curPolyASite = mspArrayIdx(polyASiteList, ++i))
        {
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
gboolean mspCoordInPolyATail(const int coord, const MSP const *msp, const GArray *polyASiteList)
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

/* Append the contents of the given text to the GString, if the text is non-null,
 * including the given separator before the text is appended. Pass the separator
 * as an empty string if no separation is required. */
static void appendTextIfNonNull(GString *gstr, const char *separator, const char *text)
{
  if (text)
    {
      g_assert(separator);
      g_string_append_printf(gstr, "%s%s", separator, text);
    }
}

/* Append the contents of the second GString to the first GString, if the second is non-null,
 * including the given separator before the second GString is appended. Pass the separator
 * as an empty string if no separation is required. */
static void appendStringIfNonNull(GString *gstr, const char *separator, const GString *text)
{
  if (text)
    appendTextIfNonNull(gstr, separator, text->str);
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
      appendTextIfNonNull(resultStr, separator, blxSequenceGetSource(blxSeq));
      appendStringIfNonNull(resultStr, separator, blxSeq->organism);
      appendStringIfNonNull(resultStr, separator, blxSeq->geneName);
      appendStringIfNonNull(resultStr, separator, blxSeq->tissueType);
      appendStringIfNonNull(resultStr, separator, blxSeq->strain);
      
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

/* Return the Source text of a BlxSequence, if it has one (note that it gets
 * this from the first MSP and does no checking whether other MSPs have the 
 * same source or not). */
static const char *blxSequenceGetSource(const BlxSequence *seq)
{
  const char *result = NULL;
  
  if (seq && seq->mspList && g_list_length(seq->mspList) > 0 && seq->mspList->data)
    {
      const MSP const *msp = (const MSP const*)(seq->mspList->data);
      result = msp->source;
    }
  
  return result;
}


/* Return the variant name of a BlxSequence (excludes prefix but includes variant) */
const char *blxSequenceGetVariantName(const BlxSequence *seq)
{
  /* Only applicable for matches; for anything else, return the full name */
  return (seq->type == BLXSEQUENCE_MATCH ? seq->variantName : seq->fullName);
}

/* Return the display name of a BlxSequence (same as full name for now) */
const char *blxSequenceGetDisplayName(const BlxSequence *seq)
{
  return seq->fullName;
}

/* Return the short name of a BlxSequence (excludes prefix and variant number) */
const char *blxSequenceGetShortName(const BlxSequence *seq)
{
  /* Only applicable to matches */
  return (seq->type == BLXSEQUENCE_MATCH ? seq->shortName : seq->fullName);
}

/* Return the length of the given blxsequence's sequence data */
int blxSequenceGetLength(const BlxSequence *seq)
{
  return (seq && seq->sequence ? seq->sequence->len : 0);
}

/* Get the start extent of the alignments from this match sequence on the
 * given reference sequence strand */
int blxSequenceGetStart(const BlxSequence *seq, const BlxStrand strand)
{
  g_assert(strand == BLXSTRAND_FORWARD || strand == BLXSTRAND_REVERSE);
  return (strand == BLXSTRAND_FORWARD ? seq->qRangeFwd.min : seq->qRangeRev.min);
}

/* Get the end extend of the sequence on the ref sequence */
int blxSequenceGetEnd(const BlxSequence *seq, const BlxStrand strand)
{
  g_assert(strand == BLXSTRAND_FORWARD || strand == BLXSTRAND_REVERSE);
  return (strand == BLXSTRAND_FORWARD ? seq->qRangeFwd.max : seq->qRangeRev.max);
}

/* Get the sequence data for the given blxsequence */
char *blxSequenceGetSeq(const BlxSequence *seq)
{
  return (seq && seq->sequence ? seq->sequence->str : NULL);
}

gboolean blxSequenceRequiresSeqData(const BlxSequence *seq)
{
  return (seq && (seq->type == BLXSEQUENCE_MATCH || seq->type == BLXSEQUENCE_VARIATION || seq->type == BLXSEQUENCE_READ_PAIR));
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
      
      /* To do: variant name and short name are only applicable to matches so 
       * ideally we wouldn't even attempt to calculate them for other types; 
       * however, when we create an empty blxsequence we don't know the type... */
      
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


/* Compare the start position in the ref seq of two MSPs. Returns a negative value if a < b; zero
 * if a = b; positive value if a > b. Secondarily sorts by type in the order that types appear in 
 * the BlxMspType enum. */
gint compareFuncMspPos(gconstpointer a, gconstpointer b)
{
  gint result = 0;

  const MSP const *msp1 = (const MSP const*)a;
  const MSP const *msp2 = (const MSP const*)b;
  
  if (msp1->qRange.min == msp2->qRange.min)
    {
      /* Sort by type. Lower type numbers should appear first. */
      result = msp2->type - msp1->type;
    }
  else 
    {
      result = msp1->qRange.min -  msp2->qRange.min;
    }
  
  return result;
}


/* Same as compareFuncMspPos but accepts pointers to MSP pointers (which is 
 * what the GArray of MSPs holds). */
gint compareFuncMspArray(gconstpointer a, gconstpointer b)
{
  gint result = 0;
  
  const MSP const *msp1 = *((const MSP const**)a);
  const MSP const *msp2 = *((const MSP const**)b);
  
  if (msp1->qRange.min == msp2->qRange.min)
    {
      /* Sort by type. Lower type numbers should appear first. */
      result = msp2->type - msp1->type;
    }
  else 
    {
      result = msp1->qRange.min -  msp2->qRange.min;
    }
  
  return result;
}


/* Determine the type of BlxSequence that an MSP belongs to */
static BlxSequenceType getBlxSequenceTypeForMsp(const MSP const *msp)
{
  BlxSequenceType result = BLXSEQUENCE_UNSET;
  
  if (msp->type == BLXMSP_MATCH)
    {
      result = BLXSEQUENCE_MATCH;
    }
  else if (mspIsExon(msp) || mspIsIntron(msp))
    {
      result = BLXSEQUENCE_TRANSCRIPT;
    }
  else if (mspIsVariation(msp))
    {
      result = BLXSEQUENCE_VARIATION;
    }
  else if (mspIsShortRead(msp))
    {
      result = BLXSEQUENCE_READ_PAIR;
    }

  return result;
}


/* Add or create a BlxSequence struct, creating the BlxSequence if one does not
 * already exist for the MSP's sequence name. Seperate BlxSequence structs are created
 * for the forward and reverse strands of the same sequence. The passed-in sequence 
 * should always be forwards, and we reverse complement it here if we need the 
 * reverse strand. Returns the new BlxSequence */
BlxSequence* addBlxSequence(const char *name, const char *idTag, BlxStrand strand, GList **seqList, char *sequence, MSP *msp, GError **error)
{
  BlxSequence *blxSeq = NULL;
  
  if (name || idTag)
    {
      /* If this is an exon or intron the match strand is not applicable. The exon should 
       * be in the same direction as the ref seq, so use the ref seq strand. */
      if (msp && (mspIsExon(msp) || mspIsIntron(msp)))
        {
          strand = msp->qStrand;
        }
    
      /* See if this strand for this sequence already exists. Horrible hack for backwards compatibility:
       * if the msp is an exon/intron, cut off the old-style 'x' or 'i' postfix from the name, if it has one. */
      char *seqName = msp ? mspGetExonTranscriptName(msp) : g_strdup(name);
      blxSeq = findBlxSequence(*seqList, seqName, idTag, strand);
      
      if (!blxSeq)
        {
          /* Create a new BlxSequence, and take ownership of the passed in sequence (if any) */
          blxSeq = createEmptyBlxSequence(seqName, idTag, NULL);
          *seqList = g_list_prepend(*seqList, blxSeq);
          blxSeq->strand = strand;
        }
      
      if (seqName && !blxSeq->fullName)
	{
	  /* It's possible that the BlxSequence was created without a name if we found an
	   * unnamed child exon before we found the parent transcript, so set the name if we have it. */
	  blxSequenceSetName(blxSeq, seqName);
	}
      
      if (msp)
        {
          /* Add the MSP to the BlxSequence's list. Keep it sorted by position. */
          blxSeq->mspList = g_list_insert_sorted(blxSeq->mspList, msp, compareFuncMspPos);
          msp->sSequence = blxSeq;
          
          if (blxSeq->type == BLXSEQUENCE_UNSET)
            {
              blxSeq->type = getBlxSequenceTypeForMsp(msp);
            }
          else if (blxSeq->type != getBlxSequenceTypeForMsp(msp))
            {
              g_warning("Adding MSP of type %d to parent of type %d (expected parent type to be %d)\n", msp->type, blxSeq->type, getBlxSequenceTypeForMsp(msp));
            }
        }
      
      /* Add the sequence data */
      addBlxSequenceData(blxSeq, sequence, error);
      
      g_free(seqName);
    }
  else
    {
      g_set_error(error, BLX_ERROR, 1, "Sequence name or parent ID must be set.\n");
    }
  
  return blxSeq;
}



/* Add the given sequence data to a BlxSequence. Validates that the existing sequence data is
 * either null or is the same as the new data; sets the given error if not. We claim ownership
 * of the given sequence data (either the BlxSequence owns it, or we delete it if it is not 
 * required). The given sequence should always be the forward strand; we complement it ourselves
 * here if this BlxSequence requires the reverse strand. */
void addBlxSequenceData(BlxSequence *blxSeq, char *sequence, GError **error)
{
  if (!sequence)
    {
      return;
    }
  
  gboolean sequenceUsed = FALSE;
  
  if (blxSeq && blxSequenceRequiresSeqData(blxSeq))
    {
      if (!blxSeq->sequence)
        {
          /* Sequence does not yet exist, so add it */
          blxSeq->sequence = g_string_new(sequence);
          sequenceUsed = TRUE;
          }
      else if (error && *error)
        {
          /* Sequence already exists. Validate that it's the same as the existing one. */
          if (!stringsEqual(sequence, blxSeq->sequence->str, FALSE))
            {
              g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_DATA_MISMATCH, "Sequence data for '%s' does not match previously-found data.\n", blxSeq->fullName);
            }
        }
    }      
  
  if (!sequenceUsed)
    {
      g_free(sequence);
    }
}


/* returns true if the given msp should be output when piping features to dotter */
static gboolean outputMsp(const MSP const *msp, IntRange *range1, IntRange *range2)
{
  return ((msp->type == BLXMSP_FS_SEG || mspIsExon(msp) || mspIsIntron(msp) || mspIsBlastMatch(msp)) && 
          (rangesOverlap(&msp->qRange, range1) || rangesOverlap(&msp->qRange, range2))
         );
}


/* Counts the number of msps that should be output when piping to dotter */
static int countMspsToOutput(const BlxSequence const *blxSeq, IntRange *range1, IntRange *range2)
{
  int numMsps = 0;
  
  GList *mspItem = blxSeq->mspList;
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *msp = (const MSP const *)(mspItem->data);
      
      if (outputMsp(msp, range1, range2))
        ++numMsps;
    }
  
  return numMsps;
}


/* write data from the given blxsequence to the given output pipe. if the ranges
 * are given, only outputs blxsequences that overlap either range */
void writeBlxSequenceToOutput(FILE *pipe, const BlxSequence *blxSeq, IntRange *range1, IntRange *range2)
{
  gboolean outputSeq = (blxSeq && (blxSeq->type == BLXSEQUENCE_TRANSCRIPT || blxSeq->type == BLXSEQUENCE_MATCH));
  
  int numMsps = countMspsToOutput(blxSeq, range1, range2);
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
          
          if (outputMsp(msp, range1, range2))
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
  fprintf(pipe, "%d %f %f %d %d %d %d %d %d %d", 
          msp->type,
          msp->score, 
          msp->id,
          msp->phase,
//          msp->fsColor, 
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
  msp->id = g_strtod(curChar, &curChar);

  nextChar(&curChar);
  msp->phase = strtol(curChar, &curChar, 10);

//  nextChar(&curChar);
//  msp->fsColor = strtol(curChar, &curChar, 10);

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


/* Insert the given MSP into the given list */
static void insertMsp(MSP *msp, MSP **mspList, MSP **lastMsp)
{
  /* Add it to the list */
  if (!*mspList) 
    {
      /* Nothing in list yet: make this the first entry */
      *mspList = msp;
    }
  
  if (*lastMsp)
    {
      /* Tag it on to the end of the list */
      (*lastMsp)->next = msp;
    }
  
  /* Make the 'lastMsp' pointer point to the new end of the list */
  *lastMsp = msp;
}


/* Allocate memory for an MSP and initialise all its fields to relevant 'empty' values.
 * Add it into the given MSP list and make the end pointer ('lastMsp') point to the new
 * end of the list. NB the msplist is obsolete and is being replaced by featureLists, but we don't
 * add it to a feature list yet because we don't know its type. Returns a pointer to the newly-created MSP */
MSP* createEmptyMsp(MSP **lastMsp, MSP **mspList)
{
  MSP *msp = (MSP *)g_malloc(sizeof(MSP));
  
  int i = 0;
  for ( ; i < BLXMODEL_NUM_MODELS; ++i)
    msp->treePaths[i] = NULL;

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
  msp->displayRange.min = 0;
  msp->displayRange.max = 0;
  msp->fullRange.min = 0;
  msp->fullRange.max = 0;
  
  msp->sSequence = NULL;
  msp->sname = NULL;
  
  msp->sRange.min = 0;
  msp->sRange.max = 0;
  msp->fullSRange.min = 0;
  msp->fullSRange.max = 0;
  
  msp->desc = NULL;
  msp->source = NULL;
  
  msp->style = NULL;
  
  msp->fs = NULL;
  msp->fsColor = 0;
  msp->fsShape = BLXCURVE_BADSHAPE;
  
  msp->xy = NULL;
  msp->gaps = NULL;
  
  insertMsp(msp, mspList, lastMsp);
  
  return msp;
}


static void freeStringPointer(char **ptr)
{
  if (*ptr)
    {
      g_free(*ptr);
      *ptr = NULL;
    }
}


/* Free all of the memory used by an MSP */
void destroyMspData(MSP *msp)
{
  freeStringPointer(&msp->qname);
  freeStringPointer(&msp->sname);
  freeStringPointer(&msp->desc);
  freeStringPointer(&msp->source);
  freeStringPointer(&msp->url);
  
  if (msp->gaps)
    {
      g_slist_free(msp->gaps);
      msp->gaps = NULL;
    }
  
  if (msp->xy)
    {
      g_array_free(msp->xy, TRUE);
      msp->xy = NULL;
    }
}


/* Allocates memory for an MSP and initialise all its fields to the given values.
 * For backwards compatibility, adds it into the given MSP list and makes the end pointer
 * ('lastMsp') point to the new end of the list. We will hopefully get rid of mspList eventually
 * and replace it by featureLists. The new msp is added to the relevant list in the featureLists
 * array according to its type. Returns a pointer to the newly-created MSP. Also creates a BlxSequence
 * for this MSP's sequence name (or adds the MSP to the existing one, if it exists already), 
 * and adds that BlxSequence to the given seqList. Takes ownership of 'sequence'. */
MSP* createNewMsp(GArray* featureLists[],
                  MSP **lastMsp, 
                  MSP **mspList,
                  GList **seqList,
                  const BlxMspType mspType,
		  const char *source,
                  const gdouble score,
                  const gdouble percentId,
                  const int phase,
                  const char *url,
		  const char *idTag,
                  const char *qName,
                  const int qStart,
                  const int qEnd,
                  const BlxStrand qStrand,
                  const int qFrame,
                  const char *sName,
                  const int sStart,
                  const int sEnd,
                  BlxStrand sStrand,
                  char *sequence,
                  GError **error)
{
  MSP *msp = createEmptyMsp(lastMsp, mspList);
  
  msp->type = mspType;
  msp->score = score; 
  msp->id = percentId; 
  msp->phase = phase;
  msp->source = g_strdup(source);
  msp->url = g_strdup(url);
  
  msp->qname = qName ? g_strdup(qName) : NULL;
  
  msp->qFrame = qFrame;
  msp->qStrand = qStrand;
  
  msp->sname = sName ? g_ascii_strup(sName, -1) : NULL;
  
  intrangeSetValues(&msp->qRange, qStart, qEnd);  
  intrangeSetValues(&msp->sRange, sStart, sEnd);

  /* Add it to the relevant feature list. Use prepend because it is quicker */
  featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);
  
  /* For exons and introns, the s strand is not applicable. We always want the exon
   * to be in the same direction as the ref sequence, so set the match seq strand to be 
   * the same as the ref seq strand */
  if (mspIsExon(msp) || mspIsIntron(msp))
    {
      sStrand = qStrand;
    }
  
  /* For matches, exons and introns, add (or add to if already exists) a BlxSequence */
  if (typeIsExon(mspType) || typeIsIntron(mspType) || typeIsMatch(mspType) || typeIsShortRead(mspType) || typeIsVariation(mspType))
    {
      addBlxSequence(msp->sname, idTag, sStrand, seqList, sequence, msp, error);
    }

  if (error && *error)
    {
      prefixError(*error, "Error creating MSP (ref seq='%s' [%d - %d], match seq = '%s' [%d - %d]). ",
                  qName, qStart, qEnd, sName, sStart, sEnd);
    }
  
  return msp;
}


/* Exons and UTRs don't have phase, but we want to display them in the same reading frame
 * as the CDS, if there is one, so if we have a corresponding CDS copy its reading frame data. */
static void copyCdsReadingFrame(MSP *exon, MSP *cds, MSP *utr)
{
  if (cds && exon)
    {
      exon->qFrame = cds->qFrame;
      exon->phase = cds->phase;
    }
  
  if (cds && utr)
    {
      utr->qFrame = cds->qFrame;
      utr->phase = cds->phase;
    }
}


/* Utility used by constructTranscriptData to create a missing exon/cds/utr given 
 * two others out of the three - i.e. if we have an overlapping exon and cds we can
 * construct the corresponding utr. If created, the new msp is added to the given 
 * BlxSequence and the  MSP list. If a CDS is given and no UTR exists, assume the exon
 * spans the entire CDS (and similarly if a UTR is given but no CDS exists) */
static void createMissingExonCdsUtr(MSP **exon, MSP **cds, MSP **utr, 
                                    BlxSequence *blxSeq, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList, 
                                    GError **error)
{
  BlxMspType newType = BLXMSP_INVALID; /* only set this if we need to construct an MSP */
  MSP **ptrToUpdate = NULL;
  int newStart = UNSET_INT;
  int newEnd = UNSET_INT;
  int newFrame = 0;
  BlxStyle *newStyle = NULL;
  char *qname = NULL;
  
  if (!*exon && (*cds || *utr))
    {
      ptrToUpdate = exon;
      newType = BLXMSP_EXON; /* always create something */
      
      /* The exon spans both the cds and utr */
      if (*cds && *utr)
        {
          newStart = min((*cds)->qRange.min, (*utr)->qRange.min);
          newEnd = max((*cds)->qRange.max, (*utr)->qRange.max);
          newFrame = (*cds)->qFrame;
          newStyle = (*cds)->style;
          qname = (*cds)->qname ? (*cds)->qname : (*utr)->qname;
        }
      else if (*cds)
        {
          newStart = (*cds)->qRange.min;
          newEnd = (*cds)->qRange.max;
          newFrame = (*cds)->qFrame;
          newStyle = (*cds)->style;
          qname = (*cds)->qname;
        }
      else
        {
          newStart = (*utr)->qRange.min;
          newEnd = (*utr)->qRange.max;
          newFrame = (*utr)->qFrame;
          newStyle = (*utr)->style;
          qname = (*utr)->qname;
        }
    }
  else if (!*cds && *exon && *utr)
    {
      ptrToUpdate = cds;
      newStyle = (*exon)->style;
      qname = (*exon)->qname ? (*exon)->qname : (*utr)->qname;
      newFrame = (*exon)->qFrame;
      
      /* The cds is the range of the exon that is not in the utr */
      if ((*utr)->qRange.max < (*exon)->qRange.max)
        {
          newType = BLXMSP_CDS;
          newStart = (*utr)->qRange.max + 1;
          newEnd = (*exon)->qRange.max;
        }
      else if ((*exon)->qRange.min < (*utr)->qRange.min)
        {
          newType = BLXMSP_CDS;
          newStart = (*exon)->qRange.min;
          newEnd = (*utr)->qRange.min - 1;
        }
    }
  else if (!*utr && *exon && *cds)
    {
      ptrToUpdate = utr;
      newFrame = (*cds)->qFrame;
      newStyle = (*cds)->style;
      qname = (*cds)->qname ? (*cds)->qname : (*exon)->qname;
      
      /* The utr is the range of the exon that is not in the cds */
      if ((*exon)->qRange.min < (*cds)->qRange.min)
        {
          newType = BLXMSP_UTR;
          newStart = (*exon)->qRange.min;
          newEnd = (*cds)->qRange.min - 1;
        }
      else if ((*cds)->qRange.max < (*exon)->qRange.max)
        {
          newType = BLXMSP_UTR;
          newStart = (*cds)->qRange.max + 1;
          newEnd = (*exon)->qRange.max;
        }
    }
  else if (*exon)
    {
      /* We just have an exon. Assume it is all utr. */
      newType = BLXMSP_UTR;
      newStart = (*exon)->qRange.min;
      newEnd = (*exon)->qRange.max;
      newFrame = (*exon)->qFrame;
      newStyle = (*exon)->style;
      qname = (*exon)->qname;
      ptrToUpdate = utr;
    }
  
  if (newType != BLXMSP_INVALID)
    {
      /* Create the new exon/cds/utr */
      DEBUG_OUT("Creating MSP for transcript '%s' of type %d.\n", blxSeq->fullName, newType);
      
      GError *tmpError = NULL;
      
      *ptrToUpdate = createNewMsp(featureLists, lastMsp, mspList, seqList, newType, NULL, UNSET_INT, UNSET_INT, UNSET_INT, NULL, blxSeq->idTag,
                                  qname, newStart, newEnd, blxSeq->strand, newFrame, blxSeq->fullName,
                                  UNSET_INT, UNSET_INT, blxSeq->strand, NULL, &tmpError);
      
      (*ptrToUpdate)->style = newStyle;
      
      if (tmpError)
        {
          prefixError(tmpError, "Error constructing missing exon/cds/utr [type='%d']", newType);
          g_propagate_error(error, tmpError);
        }
    }
  
  /* We should now have all the bits. Set the relationship data. */
  if (*exon && (*cds || *utr))
    {
      if (*cds)
        (*exon)->childMsps = g_list_append((*exon)->childMsps, *cds);
      
      if (*utr)
        (*exon)->childMsps = g_list_append((*exon)->childMsps, *utr);
    }
}


/* Construct any missing transcript data, i.e.
 *   - if we have a transcript and exons we can construct the introns;
 *   - if we have exons and CDSs we can construct the UTRs */
static void constructTranscriptData(BlxSequence *blxSeq, GArray* featureLists[], MSP **lastMsp, MSP **mspList, GList **seqList)
{
  GError *tmpError = NULL;
  
  const MSP *prevMsp = NULL;
  const MSP *prevExon = NULL;
  
  MSP *curExon = NULL;
  MSP *curCds = NULL;
  MSP *curUtr = NULL;
  
  /* Loop through all MSPs on this sequence (which must be sorted by position on the
   * ref seq - createNewMsp automatically sorts them for us) and create any missing 
   * exon/cds/utr/introns. */
  GList *mspItem = blxSeq->mspList;
  gboolean finished = FALSE;
  
  while (!finished)
    {
      MSP *msp = mspItem ? (MSP*)(mspItem->data) : NULL;
      
      /* Only consider exons and introns */
      if (mspIsExon(msp) || mspIsIntron(msp) || !msp)
        {
          /* See if there was a gap between this exon and the previous one. There's a gap if
           * we have two exons with space between them, or if we're at the first or last exon 
           * and there's a gap to the end of the transcript. */
          gboolean foundGap = FALSE;
          
          if (msp && 
              ((prevMsp && mspIsExon(prevMsp) && !rangesOverlap(&prevMsp->qRange, &msp->qRange)) ||
               (!prevMsp && blxSequenceGetStart(blxSeq, blxSeq->strand) < msp->qRange.min)))
            {
              foundGap = TRUE;
            }
          
          if (foundGap || msp == NULL)
            {
              /* We've found a gap between exons, or reached the end. First, see if the current exon/cds or utr
               * is missing and construct it if possible. Also do this if we're at the last MSP. */
              createMissingExonCdsUtr(&curExon, &curCds, &curUtr, blxSeq, featureLists, lastMsp, mspList, seqList, &tmpError);
              reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
              copyCdsReadingFrame(curExon, curCds, curUtr);
              
              IntRange newRange = {UNSET_INT, UNSET_INT};
              
              if (prevExon && curExon && !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  /* Create an intron to span the gap */
                  newRange.min = prevExon->qRange.max + 1;
                  newRange.max = curExon->qRange.min - 1;
                }
              else if (!prevExon && curExon && blxSequenceGetStart(blxSeq, blxSeq->strand) < curExon->qRange.min && 
		       !mspIsIntron(msp) && !mspIsIntron(prevMsp))
                {
                  /* Create an intron at the start */
                  newRange.min = blxSequenceGetStart(blxSeq, blxSeq->strand);
                  newRange.max = curExon->qRange.min - 1;
                }
              else if (msp == NULL && curExon && blxSequenceGetEnd(blxSeq, blxSeq->strand) > curExon->qRange.max &&
		       !mspIsIntron(prevMsp))
                {
                  /* Create an intron at the end */
                  newRange.min = curExon->qRange.max + 1;
                  newRange.max = blxSequenceGetEnd(blxSeq, blxSeq->strand);
                }
              
              if (curExon && newRange.min != UNSET_INT && newRange.max != UNSET_INT)
                {
                  createNewMsp(featureLists, lastMsp, mspList, seqList, BLXMSP_INTRON, curExon->source, 
                               curExon->score, curExon->id, 0, g_strdup(curExon->url), blxSeq->idTag, 
                               curExon->qname, newRange.min, newRange.max, blxSeq->strand, curExon->qFrame, 
                               blxSeq->fullName, UNSET_INT, UNSET_INT, blxSeq->strand, NULL, &tmpError);
                  
                  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
                }
              
              /* We're moving past the gap, so reset the pointers */
              prevExon = curExon;
              curExon = NULL;
              curCds = NULL;
              curUtr = NULL;
            }
          
          if (msp && msp->type == BLXMSP_EXON)
            curExon = msp;
          else if (msp && msp->type == BLXMSP_CDS)
            curCds = msp;
          else if (msp && msp->type == BLXMSP_UTR)
            curUtr = msp;
          
          /* Remember the last MSP we saw */
          prevMsp = msp;
          
          /* Proceed to the next MSP. We allow an extra loop with a NULL mspItem at the end, and then finish. */
          if (mspItem)
            mspItem = mspItem->next;
          else
            finished = TRUE;
        } 
      else
        {
          /* Something that's not an exon/intron. Skip this BlxSequence. */
          finished = TRUE;
        }
    }
}


/* Adjust the MSP's q coords (as parsed from the input file) by the given offset i.e. to convert
 * them to real coords */
static void adjustMspCoordsByOffset(MSP *msp, const int offset)
{
  if (offset != 0)
    {
      /* Convert the input coords (which are 1-based within the ref sequence section
       * that we're dealing with) to "real" coords (i.e. coords that the user will see). */
      msp->qRange.min += offset;
      msp->qRange.max += offset;
      
      /* Gap coords are also 1-based, so convert those too */
      GSList *rangeItem = msp->gaps;
      
      for ( ; rangeItem; rangeItem = rangeItem->next)
        {
          CoordRange *curRange = (CoordRange*)(rangeItem->data);
          curRange->qStart += offset;
          curRange->qEnd += offset;
        }
    }
}


/* Find the first/last base in an entire sequence, if not already set. For transcripts, the
 * range should already be set from the parent transcript item. For matches, get the start/end
 * of the first/last MSP in the sequence */
static void findSequenceExtents(BlxSequence *blxSeq)
{
  blxSeq->qRangeFwd.min = findMspListQExtent(blxSeq->mspList, TRUE, BLXSTRAND_FORWARD);
  blxSeq->qRangeFwd.max = findMspListQExtent(blxSeq->mspList, FALSE, BLXSTRAND_FORWARD);
  blxSeq->qRangeRev.min = findMspListQExtent(blxSeq->mspList, TRUE, BLXSTRAND_REVERSE);
  blxSeq->qRangeRev.max = findMspListQExtent(blxSeq->mspList, FALSE, BLXSTRAND_REVERSE);
}


/* Should be called after all parsed data has been added to a BlxSequence. Calculates summary
 * data and the introns etc. */
void finaliseBlxSequences(GArray* featureLists[], MSP **mspList, GList **seqList, const int offset)
{
  /* Loop through all MSPs and adjust their coords by the offest. Also find
   * the last MSP in the list. */
  MSP *msp = *mspList;
  MSP *lastMsp = msp;

  while (msp)
    {
      adjustMspCoordsByOffset(msp, offset);
      msp = msp->next;
      if (msp)
        lastMsp = msp;
    }
  
  /* Loop through all BlxSequences */
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);

      /* So far we only have the forward strand version of each sequence. We must complement any 
       * that need the reverse strand */
      if (blxSeq && 
          blxSeq->strand == BLXSTRAND_REVERSE && 
          (blxSeq->type == BLXSEQUENCE_MATCH || (blxSeq->type == BLXSEQUENCE_READ_PAIR && !ALL_READ_PAIRS_FWD)) && 
          blxSeq->sequence && 
          blxSeq->sequence->str)
        {
          blxComplement(blxSeq->sequence->str);
        }
      
      findSequenceExtents(blxSeq);
      constructTranscriptData(blxSeq, featureLists, &lastMsp, mspList, seqList);
    }
}


/* Given a list of msps, find the min or max q coords (according to the given flag)
 * for msps on the given q strand (or on both strands if strand is BLXSTRAND_NONE. */
int findMspListQExtent(GList *mspList, const gboolean findMin, const BlxStrand strand)
{
  int result = UNSET_INT;
  gboolean first = TRUE;
  
  GList *mspItem = mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *msp = (const MSP const*)(mspItem->data);
    
      if (msp->qStrand == strand || strand == BLXSTRAND_NONE)
	{      
	  if (first)
	    {
	      result = findMin ? msp->qRange.min : msp->qRange.max;
	      first = FALSE;
	    }
	  else if (findMin && msp->qRange.min < result)
	    {
	      result = msp->qRange.min;
	    }
	  else if (!findMin && msp->qRange.max > result)
	    {
	      result = msp->qRange.max;
	    }
	}
    }
  
  return result;
}


/* Given a list of msps, find the min or max s coords, according to the given flag. */
int findMspListSExtent(GList *mspList, const gboolean findMin)
{
  int result = UNSET_INT;
  gboolean first = TRUE;
  
  GList *mspItem = mspList;
  
  for ( ; mspItem; mspItem = mspItem->next)
    {
      const MSP const *msp = (const MSP const*)(mspItem->data);
      
      if (first)
	{
	  result = findMin ? msp->sRange.min : msp->sRange.max;
	  first = FALSE;
	}
      else if (findMin && msp->sRange.min < result)
	{
	  result = msp->sRange.min;
	}
      else if (!findMin && msp->sRange.max > result)
	{
	  result = msp->sRange.max;
	}
    }
  
  return result;
}


