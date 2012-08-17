/*  File: blxmsp.c
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


static const char* blxSequenceGetValueAsString(const BlxSequence *seq, const int columnId);




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

gboolean typeIsTranscript(const BlxMspType mspType)
{
  return (mspType == BLXMSP_TRANSCRIPT);
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

gboolean typeIsRegion(const BlxMspType mspType)
{
  return (mspType == BLXMSP_REGION);
}

/* This returns true if the given type should be shown in the detail-view */
gboolean typeShownInDetailView(const BlxMspType mspType)
{
  return (mspType == BLXMSP_MATCH || mspType == BLXMSP_SHORT_READ || mspType == BLXMSP_CDS || mspType == BLXMSP_UTR);
}

/* This returns true if the given sequence should be shown in the detail-view */
gboolean blxSequenceShownInDetailView(const BlxSequence *blxSeq)
{
  return (blxSeq->type == BLXSEQUENCE_MATCH || blxSeq->type == BLXSEQUENCE_TRANSCRIPT);
}

/* This returns true if the given sequence should be shown in the big picture grids */
gboolean blxSequenceShownInGrid(const BlxSequence *blxSeq)
{
  return (blxSeq->type == BLXSEQUENCE_MATCH || blxSeq->type == BLXSEQUENCE_READ_PAIR);
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
                                         const int defaultColorId,
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
      result = mspGetColor(prevExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (nextIsUtr)
    {
      result = mspGetColor(nextExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (prevExon)
    {
      /* Both exons (or the sole exon, if only one exists) are CDS, so the intron is CDS */
      result = mspGetColor(prevExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
    }
  else if (nextExon)
    {
      /* This is the only exon and it is CDS, so make the intron CDS */
      result = mspGetColor(nextExon, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
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
  gboolean result = mspIsVariation(msp);
  
  if (result)
    {
      const char *seq = mspGetMatchSeq(msp);
      result = (seq && seq[0] == '-');
    }

  return result;
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


/* Return the match sequence name. (Gets it from the BlxSequence if the MSP itself doesn't have
 * a name) */
const char *mspGetSName(const MSP const *msp)
{
  const char *result = NULL;
  
  if (msp)
    {
      if (msp->sname && msp->sname[0] != 0)
        {
          result = msp->sname;
        }
      else if (msp->sSequence)
        {
          result = blxSequenceGetName(msp->sSequence);
        }
    }
  
  return result;
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

/* Return the reference sequence name that this msp aligns against */
const char* mspGetRefName(const MSP const *msp)
{
  return msp->qname;
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
  return (msp ? blxSequenceGetSequence(msp->sSequence) : NULL);
}

/* Get the source of the MSP */
const char* mspGetSource(const MSP const *msp)
{
  return blxSequenceGetSource(msp->sSequence);
}


/* Return the color matching the given properties from the given style */
static const GdkColor *styleGetColor(const BlxStyle *style, 
                                     const gboolean selected, 
                                     const gboolean usePrintColors,
                                     const gboolean fill, 
                                     const gboolean utr,
                                     GArray *defaultColors,
                                     const int defaultColorId)
{
  const GdkColor *result = NULL;
  const BlxColor *blxColor = NULL;
  
  if (fill)
    {
      if (utr)
        blxColor = &style->fillColorUtr;
      else
        blxColor = &style->fillColor;
    }
  else
    {
      if (utr)
        blxColor = &style->lineColorUtr;
      else
        blxColor = &style->lineColor;
    }

  /* If it's transparent, use the background color instead, unless 
   * selected is true in which case we need to use the highlight color */
  if (blxColor && blxColor->transparent && !selected)
    result = getGdkColor(defaultColorId, defaultColors, selected, usePrintColors);
  else if (blxColor)
    result = blxColorGetColor(blxColor, selected, usePrintColors);

  return result;
}


/* Get the color for drawing the given MSP (If 'selected' is true, returns
 * the color when the MSP is selected.). Returns the fill color if 'fill' is 
 * true, otherwise the line color */
const GdkColor* mspGetColor(const MSP const *msp, 
			    GArray *defaultColors, 
                            const int defaultColorId,
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
      result = styleGetColor(msp->style, selected, usePrintColors, fill, msp->type == BLXMSP_UTR, defaultColors, defaultColorId);
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
          result = mspGetIntronColor(msp, defaultColors, defaultColorId, blxSeq, selected, usePrintColors, fill, exonFillColorId, exonLineColorId, cdsFillColorId, cdsLineColorId, utrFillColorId, utrLineColorId);
	break;
	
	default:
	break;
      };
    }
  
  return result;
}


/* Get functions from msps for various sequence properties. Result is owned by the sequence
 * and should not be free'd. Returns NULL if the property is not set. */
const char *mspGetOrganism(const MSP const *msp)
{
  return (msp ? blxSequenceGetOrganism(msp->sSequence) : NULL);
}

const char *mspGetGeneName(const MSP const *msp)
{
  return (msp ? blxSequenceGetGeneName(msp->sSequence) : NULL);
}

const char *mspGetTissueType(const MSP const *msp)
{
  return (msp ? blxSequenceGetTissueType(msp->sSequence) : NULL);
}

const char *mspGetStrain(const MSP const *msp)
{
  return (msp ? blxSequenceGetStrain(msp->sSequence) : NULL);
}


/* Return the coords of an MSP as a string. The result should be free'd with g_free */
char *mspGetCoordsAsString(const MSP const *msp)
{
  char *result = NULL;
  
  if (msp)
    {
    GString *resultStr = g_string_new("");
    
    g_string_append_printf(resultStr, "%d - %d [%d - %d]", msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
    
    result = g_string_free(resultStr, FALSE);
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

/* Return summary info about a given BlxSequence (e.g. for displaying in the status bar). The
 * result should be free'd with g_free. */
char* blxSequenceGetSummaryInfo(const BlxSequence const *blxSeq, GList *columnList)
{
  char *result = NULL;
  
  if (blxSeq)
    {
      GString *resultStr = g_string_new("");
      const char *separator = "";

      /* Loop through all the columns, appending any that should be shown
       * to the result string */
      GList *item;

      for ( ; item; item = item->next)
        {
          BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
          
          if (columnInfo->showSummary)
            {
              const char *valueText = blxSequenceGetValueAsString(blxSeq, columnInfo->columnId);

              if (valueText)
                {
                  appendTextIfNonNull(resultStr, separator, valueText);

                  /* The first item has no separator, but subsequence items do,
                   * so set the separator after the first item has been added. */
                  separator = ";  ";
                }
            }
        }
      
      result = g_string_free(resultStr, FALSE);
    }
  
  return result;
}


/* Return the full name of a BlxSequence (including prefix and variant) */
const char *blxSequenceGetName(const BlxSequence *seq)
{
  const char *result = NULL;
  
  if (seq)
    {
      result = blxSequenceGetValueAsString(seq, BLXCOL_SEQNAME);
      
      if (!result && seq->idTag)
        result = seq->idTag;

      if (!result)
        g_error("Sequence does not have a name specified.\n");
    }
  
  return result;
}

/* Return the Source text of a BlxSequence, if it has one (note that it gets
 * this from the first MSP and does no checking whether other MSPs have the 
 * same source or not). */
const char *blxSequenceGetSource(const BlxSequence *seq)
{
  const char *result = NULL;
  
  if (seq)
    result = blxSequenceGetValueAsString(seq, BLXCOL_SOURCE);
  
  return result;
}

gboolean blxSequenceGetLinkFeatures(const BlxSequence *seq, const gboolean defaultLinkFeatures)
{
  gboolean result = defaultLinkFeatures;
  
  if (seq && seq->dataType)
    result = seq->dataType->linkFeaturesByName;
  
  return result;
}

/* Return the fetch method of a BlxSequence. If 'bulk' is true,
 * get the bulk-fetch method, otherwise the user-fetch method.
 * 'index' indicates which method to choose if multiple methods
 * are available; 0 is the first (preferred) method, 1 the second 
 * etc. If no fetch method is set, return the given default method 
 * instead (for index==0 only). Pass 'defaultMethod' as 0 if N/A. */
GQuark blxSequenceGetFetchMethod(const BlxSequence *seq, 
                                 const gboolean bulk,
                                 const gboolean optionalColumns,
                                 const int index,
                                 const GArray *defaultMethods)
{
  GQuark result = 0;

  if (seq && seq->dataType)
    {
      GArray *array = NULL;

      if (bulk && optionalColumns)
        array = seq->dataType->optionalFetch;
      else if (bulk)
        array = seq->dataType->bulkFetch;
      else
        array = seq->dataType->userFetch;
      
      if (array && index >= 0 && index < array->len)
        result = g_array_index(array, GQuark, index);
    }

  if (!result && defaultMethods && index >= 0 && index < defaultMethods->len) 
    {
      result = g_array_index(defaultMethods, GQuark, index);
    }
  
  return result;
}

/* Return the length of the given blxsequence's sequence data */
int blxSequenceGetLength(const BlxSequence *seq)
{
  int result = 0;
  
  if (seq)
    {
      const char *sequence = blxSequenceGetSequence(seq);

      if (sequence)
        result = strlen(sequence);
    }
  
  return result;
}

/* Get the start extent of the alignments from this match sequence on the
 * given reference sequence strand */
int blxSequenceGetStart(const BlxSequence *seq, const BlxStrand strand)
{
  return (strand == BLXSTRAND_REVERSE ? seq->qRangeRev.min : seq->qRangeFwd.min);
}

/* Get the end extend of the sequence on the ref sequence */
int blxSequenceGetEnd(const BlxSequence *seq, const BlxStrand strand)
{
  return (strand == BLXSTRAND_FORWARD ? seq->qRangeRev.max : seq->qRangeFwd.max);
}

/* Get the sequence data for the given blxsequence */
const char *blxSequenceGetSequence(const BlxSequence *seq)
{
  const char *result = blxSequenceGetValueAsString(seq, BLXCOL_SEQUENCE);
  return result;
}

/* This returns true if the given sequence object requires sequence data
 * (i.e. the actual dna or peptide sequence string) */
gboolean blxSequenceRequiresSeqData(const BlxSequence *seq)
{
  return (seq && 
          (seq->type == BLXSEQUENCE_MATCH || 
           seq->type == BLXSEQUENCE_VARIATION || 
           seq->type == BLXSEQUENCE_READ_PAIR || 
           seq->type == BLXSEQUENCE_REGION));
}

/* This returns true if the given sequence object requires optional data such
 * as organism and tissue-type, if this data is available. */
gboolean blxSequenceRequiresOptionalData(const BlxSequence *seq)
{
  return (seq && seq->type == BLXSEQUENCE_MATCH);
}

/* Get the value for the given column */
GValue* blxSequenceGetValue(const BlxSequence *seq, const int columnId)
{
  GValue *result = NULL;
  
  if (seq)
    result = &g_array_index(seq->values, GValue, columnId);

  return result;
}

void blxSequenceSetValue(const BlxSequence *seq, const int columnId, GValue *value_in)
{
  if (seq && seq->values)
    {
      GValue *value = blxSequenceGetValue(seq, columnId);
      g_value_reset(value);
      g_value_copy(value_in, value);
    }
}

/* Set the value for the given column. The given string is converted to
 * the relevant value type */
void blxSequenceSetValueFromString(const BlxSequence *seq, const int columnId, const char *inputStr)
{
  if (inputStr && seq && seq->values)
    {
      GValue *value = blxSequenceGetValue(seq, columnId);
      g_value_reset(value);

      if (G_VALUE_HOLDS_STRING(value))
        {
          g_value_take_string(value, g_strdup(inputStr));
        }
      else if (G_VALUE_HOLDS_INT(value))
        {
          const int tmp = atoi(inputStr);
          g_value_set_int(value, tmp);
        }
      else if (G_VALUE_HOLDS_DOUBLE(value))
        {
          const gdouble tmp = g_strtod(inputStr, NULL);
          g_value_set_double(value, tmp);
        }
      else
        {
          g_warning("Tried to set value of unknown type for column '%d'\n", columnId);
        }
    }

}

/* Get the string value for the given column. Returns null if the
 * value is not a string type. */
static const char* blxSequenceGetValueAsString(const BlxSequence *seq, const int columnId)
{
  const char *result = NULL;
  
  GValue *value = blxSequenceGetValue(seq, columnId);

  if (value)
    {
      if (G_VALUE_HOLDS_STRING(value))
        result = g_value_get_string(value);
    }
  
  /* Return null if it's an empty value (i.e. if it's unset) */
  if (result && *result == 0)
    result = NULL;

  return result;
}

const char *blxSequenceGetOrganism(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_ORGANISM);
}

const char *blxSequenceGetGeneName(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_GENE_NAME);
}

const char *blxSequenceGetTissueType(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_TISSUE_TYPE);
}

const char *blxSequenceGetStrain(const BlxSequence *seq)
{
  return blxSequenceGetValueAsString(seq, BLXCOL_STRAIN);
}

/* Return the sequence data as a string in fasta format.
 * Result should be free'd with g_free, unless the name
 * or sequence couldn't be found, in which case the result
 * is null */
char *blxSequenceGetFasta(const BlxSequence *seq)
{
  char *result = NULL;
  
  if (seq)
    {
      const char *name = blxSequenceGetName(seq);
      const char *sequence = blxSequenceGetSequence(seq);
      
      if (name && sequence)
        {
          result = g_strdup_printf(">%s\n%s", name, sequence);
        }
    }
  
  return result;
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
    
    g_string_append_printf(resultStr, "SEQUENCE NAME:\t%s%c%c", blxSequenceGetName(blxSeq), strand, separator);
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
    
    result = g_string_free(resultStr, FALSE);
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
      const char *curName = blxSequenceGetName(curSeq);
    
      if (name && stringsEqual(curName, name, FALSE))
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
  
  const char *variantName = blxSequenceGetName(variant);
  
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


/* Destroy all of the BlxSequences */
void destroyBlxSequenceList(GList **seqList)
{
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);
      destroyBlxSequence(blxSeq);
    }
  
  g_list_free(*seqList);
  *seqList = NULL;
}


/* Frees all memory used by a BlxSequence */
void destroyBlxSequence(BlxSequence *seq)
{
  if (seq)
    {
      /* Free all column values */
      int i = 0;      
      for ( ; i < seq->values->len; ++i)
        {
          GValue *value = &g_array_index(seq->values, GValue, i);
          g_value_unset(value);
        }
      
      if (seq->idTag)
        g_free(seq->idTag);
      
      g_free(seq);
    }
}


/* Set the value for a particular column (by column name) */
void blxSequenceSetColumn(BlxSequence *seq, const char *colName, const char *value, GList *columnList)
{
  if (!colName || !value)
    return;
  
  gboolean found = FALSE;

  /* Loop through the column list and find the one with this name */
  GList *item = columnList;
  for ( ; item && !found; item = item->next)
    {
      BlxColumnInfo *columnInfo = (BlxColumnInfo*)(item->data);
      
      if (stringsEqual(colName, columnInfo->title, FALSE))
        {
          blxSequenceSetValueFromString(seq, columnInfo->columnId, value);
          found = TRUE;
        }
    }
  
  if (!found)
    g_warning("Unrecognised column '%s'\n", colName);
}


/* Utility to create a BlxSequence with the given name. */
BlxSequence* createEmptyBlxSequence()
{
  BlxSequence *seq = g_malloc(sizeof(BlxSequence));
  
  seq->type = BLXSEQUENCE_UNSET;
  seq->dataType = NULL;
  seq->idTag = NULL;
  seq->mspList = NULL;
  seq->sequenceReqd = FALSE;
  seq->values = NULL;

  return seq;
}



BlxDataType* createBlxDataType()
{
  BlxDataType *result = g_malloc(sizeof *result);
  
  result->name = 0;
  result->bulkFetch = NULL;
  result->userFetch = NULL;
  result->linkFeaturesByName = FALSE;
  
  return result;
}

void destroyBlxDataType(BlxDataType **blxDataType)
{
  if (!blxDataType)
    return;
  
  g_free((*blxDataType));
  *blxDataType = NULL;
}

/* The following functions get the data type values as strings */
const char* getDataTypeName(BlxDataType *blxDataType)
{
  return blxDataType ? g_quark_to_string(blxDataType->name) : NULL;
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
      
      stringProtect(pipe, blxSequenceGetName(blxSeq));
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
  BlxSequence *blxSeq = createEmptyBlxSequence();
  
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
  
  char *fullName = stringUnprotect(&curChar, NULL);
  blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQNAME, fullName);
  g_free(fullName);

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


/* Destroy all of the MSPs in the given list */
void destroyMspList(MSP **mspList)
{
  /* Free the allocated sequences and names */
  MSP *msp = NULL;
  for (msp = *mspList; msp; msp = msp->next)
    {
      destroyMspData(msp);
    }
  
  /* Now free the MSPs themselves. */
  MSP *fmsp = NULL;
  for (msp = *mspList; msp; )
    {
      fmsp = msp;
      msp = msp->next;
      g_free(fmsp);
    }
  
  *mspList = NULL;
  
  return ;
}


/* Free all of the memory used by an MSP */
void destroyMspData(MSP *msp)
{
  freeStringPointer(&msp->qname);
  freeStringPointer(&msp->sname);
  freeStringPointer(&msp->desc);
  
  if (msp->gaps)
    {
      /* free the child msp list */
      if (msp->childMsps)
        {
          g_list_free(msp->childMsps);
          msp->childMsps = NULL;
        }

      /* free memory allocated for the gap ranges */
      GSList *item = msp->gaps;
      for ( ; item; item = item->next)
        g_free(item->data);
      
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
                  GList *columnList,
                  const BlxMspType mspType,
                  BlxDataType *dataType,
		  const char *source,
                  const gdouble score,
                  const gdouble percentId,
                  const int phase,
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
                  const gboolean linkFeaturesByName,
                  const GQuark filename,
                  GError **error)
{
  MSP *msp = createEmptyMsp(lastMsp, mspList);
  
  msp->type = mspType;
  msp->score = score; 
  msp->id = percentId; 
  msp->phase = phase;
  msp->filename = filename;
  
  msp->qname = qName ? g_strdup(qName) : NULL;
  
  msp->qFrame = qFrame;
  msp->qStrand = qStrand;
  
  msp->sname = sName ? g_strdup(sName) : NULL;
  //g_ascii_strup(sName, -1)
  
  intrangeSetValues(&msp->qRange, qStart, qEnd);  
  intrangeSetValues(&msp->sRange, sStart, sEnd);

  /* For exons and introns, the s strand is not applicable. We always want the exon
   * to be in the same direction as the ref sequence, so set the match seq strand to be 
   * the same as the ref seq strand */
  if (mspIsExon(msp) || mspIsIntron(msp))
    {
      sStrand = qStrand;
    }
  
  /* For matches, exons and introns, add a new (or add to an existing) BlxSequence */
  if (typeIsExon(mspType) || typeIsIntron(mspType) || 
      typeIsMatch(mspType) || typeIsShortRead(mspType) || 
      typeIsVariation(mspType) || typeIsRegion(mspType))
    {
      addBlxSequence(msp->sname, idTag, sStrand, dataType, source, seqList, columnList, sequence, msp, linkFeaturesByName, error);
    }

  /* Add it to the relevant feature list. */
  featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);

  if (error && *error)
    {
      prefixError(*error, "Error creating MSP (ref seq='%s' [%d - %d], match seq = '%s' [%d - %d]). ",
                  qName, qStart, qEnd, sName, sStart, sEnd);
    }
  
  return msp;
}


/* Make a copy of an MSP */
MSP* copyMsp(const MSP const *src,
             GArray* featureLists[],             
             MSP **lastMsp, 
             MSP **mspList,
             GList **seqList,
             GError **error)
{
  MSP *msp = createEmptyMsp(lastMsp, mspList);
  
  msp->type = src->type;
  msp->score = UNSET_INT; 
  msp->id = UNSET_INT; 
  msp->phase = src->phase;
  
  msp->qname = src->qname ? g_strdup(src->qname) : NULL;
  
  msp->qFrame = src->qFrame;
  msp->qStrand = src->qStrand;
  
  msp->sname = src->sname ? g_strdup(src->sname) : NULL;
  
  intrangeSetValues(&msp->qRange, src->qRange.min, src->qRange.max);  
  intrangeSetValues(&msp->sRange, src->sRange.min, src->sRange.max);
  
  /* For matches, exons and introns, add (or add to if already exists) a BlxSequence */
  if (src->sSequence)
    {
      src->sSequence->mspList = g_list_insert_sorted(src->sSequence->mspList, msp, compareFuncMspPos);
      msp->sSequence = src->sSequence;
    }

  /* Add it to the relevant feature list. */
  featureLists[msp->type] = g_array_append_val(featureLists[msp->type], msp);

  if (error && *error)
    {
      prefixError(*error, "Error creating MSP (ref seq='%s' [%d - %d], match seq = '%s' [%d - %d]). ",
                  src->qname, src->qRange.min, src->qRange.max, src->sname, src->sRange.min, src->sRange.max);
    }
  
  return msp;
}



/* Set the given child list in the given exon. Takes ownership of the child list
 * (and frees it if exon is null).
 * Exons and UTRs don't have phase, but we want to display them in the same reading frame
 * as the CDS in the same exon, if there is one; this function copies it to its siblings. */
static void setExonChildList(MSP *exon, GList *childList)
{
  if (!exon)
    {
      g_list_free(childList);
      return;
    }

  /* Take ownership of the child list */
  if (exon->childMsps)
    g_list_free(exon->childMsps);
  
  exon->childMsps = childList;
  
  /* Loop through and see if there's a CDS */
  int frame = UNSET_INT;
  int phase = UNSET_INT;
  gboolean found = FALSE;
  
  GList *childItem = exon->childMsps;
  for ( ; childItem; childItem = childItem->next)
    {
      MSP *msp = (MSP*)(childItem->data);
      
      if (msp->type == BLXMSP_CDS)
        {
          frame = msp->qFrame;
          phase = msp->phase;
          found = TRUE;
          break;
        }
    }
  
  if (found)
    {
      /* Update the exon */
      exon->qFrame = frame;
      exon->phase = phase;
      
      /* Loop through and update the other msps */
      for (childItem = exon->childMsps; childItem; childItem = childItem->next)
        {
          MSP *msp = (MSP*)(childItem->data);
          msp->qFrame = frame;
          msp->phase = phase;
        }
    }
}


/* Utility to create the msp for the createMissingExonCdsUtr function. */
static MSP* createMissingMsp(const BlxMspType newType,
                             const int newStart,
                             const int newEnd,
                             const char *qname,
                             const int newFrame,
                             BlxStyle *newStyle,
                             BlxSequence *blxSeq, 
                             GArray* featureLists[], 
                             MSP **lastMsp, 
                             MSP **mspList, 
                             GList **seqList,
                             GList *columnList,
                             const gboolean linkFeatures,
                             GError **error)
{
  MSP *result = NULL;
  
  if (newType != BLXMSP_INVALID)
    {
      /* Create the new exon/cds/utr */
      DEBUG_OUT("Creating MSP for transcript '%s' of type %d.\n", blxSequenceGetName(blxSeq), newType);
      
      GError *tmpError = NULL;
      
      result = createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, newType, NULL, blxSequenceGetSource(blxSeq),
                            UNSET_INT, UNSET_INT, UNSET_INT, blxSeq->idTag,
                            qname, newStart, newEnd, blxSeq->strand, newFrame, blxSequenceGetName(blxSeq),
                            UNSET_INT, UNSET_INT, blxSeq->strand, NULL,
                            blxSequenceGetLinkFeatures(blxSeq, linkFeatures), 0, &tmpError);
      
      result->style = newStyle;
      
      if (tmpError)
        {
          prefixError(tmpError, "Error constructing missing exon/cds/utr [type='%d']", newType);
          g_propagate_error(error, tmpError);
        }
    }
  
  return result;
}


/* We have an exon and one or more child CDS/UTRs. Check whether there is a gap
 * at the startor end of the exon and, if so, construct a CDS/UTR to fill it. */
static void createMissingCdsUtr(MSP *exon,
                                GList **childList,
                                BlxSequence *blxSeq, 
                                GArray* featureLists[], 
                                MSP **lastMsp, 
                                MSP **mspList, 
                                GList **seqList, 
                                GList *columnList,
                                const gboolean linkFeatures,
                                GError **error)
{
  MSP *startMsp = (MSP*)(g_list_first(*childList)->data);
  MSP *endMsp = (MSP*)(g_list_last(*childList)->data);
  
  if (exon->qRange.min < startMsp->qRange.min)
    {
      const BlxMspType type = (startMsp->type == BLXMSP_CDS ? BLXMSP_UTR : BLXMSP_CDS);
      MSP *result = createMissingMsp(type, exon->qRange.min, startMsp->qRange.min - 1, exon->qname, exon->qFrame, exon->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
      *childList = g_list_append(*childList, result);
    }

  if (exon->qRange.max > endMsp->qRange.max)
    {
      const BlxMspType type = (endMsp->type == BLXMSP_CDS ? BLXMSP_UTR : BLXMSP_CDS);
      MSP *result = createMissingMsp(type, endMsp->qRange.max + 1, exon->qRange.max, exon->qname, exon->qFrame, exon->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
      *childList = g_list_append(*childList, result);
    }
}


/* Create a UTR that spans the given exon, and add it to the given childList */
static void createMissingUtr(MSP *exon,
                             GList **childList,
                             BlxSequence *blxSeq, 
                             GArray* featureLists[], 
                             MSP **lastMsp, 
                             MSP **mspList, 
                             GList **seqList, 
                             GList *columnList,
                             const gboolean linkFeatures,
                             GError **error)
{
  MSP *result = createMissingMsp(BLXMSP_UTR, exon->qRange.min, exon->qRange.max, exon->qname, exon->qFrame, exon->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
  *childList = g_list_append(*childList, result);
}


/* Create an exon that spans the given child CDSs/UTRs. Must not be called with
 * an empty childList. */
static MSP* createMissingExon(GList *childList,
                              BlxSequence *blxSeq, 
                              GArray* featureLists[], 
                              MSP **lastMsp, 
                              MSP **mspList, 
                              GList **seqList, 
                              GList *columnList,
                              const gboolean linkFeatures,
                              GError **error)
{
  /* Get the max and min extent of the children. They should be in order
   * of increasing coords and should not overlap etc. */
  MSP *startMsp = (MSP*)(g_list_first(childList)->data);
  MSP *endMsp = (MSP*)(g_list_last(childList)->data);
  
  MSP *result = createMissingMsp(BLXMSP_EXON, startMsp->qRange.min, endMsp->qRange.max, startMsp->qname, startMsp->qFrame, startMsp->style, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
  return result;
}


/* Utility used by constructTranscriptData to create a missing exon/cds/utr given 
 * two others out of the three - i.e. if we have an overlapping exon and cds we can
 * construct the corresponding utr. If created, the new msp is added to the given 
 * BlxSequence and the  MSP list. If a CDS is given and no UTR exists, assume the exon
 * spans the entire CDS (and similarly if a UTR is given but no CDS exists) */
static void createMissingExonCdsUtr(MSP **exon, 
                                    GList **childList,
                                    BlxSequence *blxSeq, 
                                    GArray* featureLists[], 
                                    MSP **lastMsp, 
                                    MSP **mspList, 
                                    GList **seqList, 
                                    GList *columnList,
                                    const gboolean linkFeatures,
                                    GError **error)
{
  if (*exon && g_list_length(*childList) > 0)
    {
      /* We have an exon and one or more CDS/UTRs. Check if there any other
       * child CDS/UTRs that we need to create. */
      createMissingCdsUtr(*exon, childList, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
    }
  else if (*exon)
    {
      /* We have an exon but no children. Assume the exon is a single UTR */
      createMissingUtr(*exon, childList, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
    }
  else if (g_list_length(*childList) > 0)
    {
      /* We have children but no exon; create the exon from the children */
      *exon = createMissingExon(*childList, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, error);
    }
}


/* Construct any missing transcript data, i.e.
 *   - if we have a transcript and exons we can construct the introns;
 *   - if we have exons and CDSs we can construct the UTRs */
static void constructTranscriptData(BlxSequence *blxSeq, 
                                    GArray* featureLists[], 
                                    MSP **lastMsp,
                                    MSP **mspList,
                                    GList **seqList,
                                    GList *columnList,
                                    const gboolean linkFeatures)
{
  GError *tmpError = NULL;
  
  const MSP *prevMsp = NULL;
  const MSP *prevExon = NULL;
  
  MSP *curExon = NULL;          /* the current exon we're looking at */
  GList *curChildMsps = NULL;   /* the child CDS/UTRs of the current exon */
  
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

          if (msp && prevMsp)
            {
              /* We have a previous msp; there's a gap if the previous msp is an exon
               * and if there is a gap between it and the current exon */
              foundGap = mspIsExon(prevMsp) && !rangesOverlap(&prevMsp->qRange, &msp->qRange) && !rangesAdjacent(&prevMsp->qRange, &msp->qRange);
            }
          else if (msp)
            {
              /* No previous msp; check if we're at the first exon and if so whether
               * there's a gap between it and the start of the transcript */
              foundGap = blxSequenceGetStart(blxSeq, blxSeq->strand) < msp->qRange.min;
            }
          
          if (foundGap || msp == NULL)
            {
              /* We've found a gap between exons, or reached the end. First, see if the current exon/cds or utr
               * is missing and construct it if possible. Also do this if we're at the last MSP. */
              createMissingExonCdsUtr(&curExon, &curChildMsps, blxSeq, featureLists, lastMsp, mspList, seqList, columnList, linkFeatures, &tmpError);
              reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
              
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
                  createNewMsp(featureLists, lastMsp, mspList, seqList, columnList, BLXMSP_INTRON, NULL, blxSequenceGetSource(blxSeq), 
                               curExon->score, curExon->id, 0, blxSeq->idTag, 
                               curExon->qname, newRange.min, newRange.max, blxSeq->strand, curExon->qFrame, 
                               blxSequenceGetName(blxSeq), UNSET_INT, UNSET_INT, blxSeq->strand, NULL, 
                               blxSequenceGetLinkFeatures(blxSeq, linkFeatures), 0, &tmpError);
                  
                  reportAndClearIfError(&tmpError, G_LOG_LEVEL_CRITICAL);
                }
              
              /* We're done with this exon, so set the exon's list of child msps
               * and reset the pointers */
              setExonChildList(curExon, curChildMsps);

              prevExon = curExon;
              curExon = NULL;
              curChildMsps = NULL;
            }
          
          if (msp && msp->type == BLXMSP_EXON)
            curExon = msp;
          else if (msp && (msp->type == BLXMSP_CDS || msp->type == BLXMSP_UTR))
            curChildMsps = g_list_append(curChildMsps, msp);
          
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


/* Get the offset required from the given base to give the coord that is the first base in 
 * the first codon of reading frame 1 (or the last codon in reading frame 3 for the reverse strand) */
static int getOffsetToCodonStart(const int coord, const int numFrames, const BlxStrand strand)
{
  int offset = 0;
  
  if (strand == BLXSTRAND_FORWARD)
    {
      /* If the sequence is 
       *     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ...
       * then the calculated base below will be
       *     0, 1, 2, 0, 1, 2, 0, 1, 2, 0,  1, ...
       * so we have to offset the coords by the following values to get to a base1 coord
       *     0, 2, 1, 0, 2, 1, 0, 2, 1, 0,  2, ...
       */
      
      int base = ((coord - 1) % numFrames); /* 0, 1 or 2 */
      offset = numFrames - base;            /* 3, 2 or 1 */
      
      if (offset >= numFrames)
	offset -= numFrames;                /* 0, 2 or 1 */
    }
  else
    {
      /* I'm not sure if there is a convention for where the reading frame starts in the reverse
       * strand, so I've made up my own convention. It essentially means that we want the first 
       * coord in the reversed sequence to be the last base in the last reading frame, i.e. base 3 in frame 3.
       *
       * If the forward strand coords are
       *             1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 
       * then the have the following base numbers for each frame
       * frame 1:    1, 2, 3, 1, 2, 3, 1, 2, 3, 1,  2
       * frame 2:    3, 1, 2, 3, 1, 2, 3, 1, 2, 3,  1
       * frame 3:    2, 3, 1, 2, 3, 1, 2, 3, 1, 2,  3 (I've included up to coord 11 so we end at base 3 in frame 3)
       *
       * Looking at frame 1 in reverse we have
       *             11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
       *              2,  1, 3, 2, 1, 3, 2, 1, 3, 2, 1
       *
       * Base 3 in frame 3 is base 2 in frame 1, so this is where we want to start. The offsets we need
       * to get to a frame1-base2 coord are:
       *              0,  2, 1, 0, 2, 1, 0, 2, 1, 0, 2
       */
      
      int base = (coord + 1) % numFrames;   /* 0, 2 or 1 */
      offset = base;
    }
  
  return offset;
}


/* Calculate the reference sequence reading frame that the given MSP belongs to, if not
 * already set. Requires either the phase to be set, or the frame to already be set; otherwise
 * assumes a phase of 0 and gives a warning. */
static void calcReadingFrame(MSP *msp, const BlxSeqType seqType, const int numFrames, const IntRange const *refSeqRange)
{
  /* For matches and exons, calculate frame if the phase is known, because the old code that 
   * used to pass the reading frame in exblx files seemed to occasionally pass an incorrect reading frame. */
  if (!mspIsIntron(msp))
    {
      /* Get the first coord of the first complete codon. This is the start coord (5' end) of the match
       * plus (or minus) the phase (if non-zero), which initially gets stored in the qFrame field in the MSP... */
      const int direction = (msp->qStrand == BLXSTRAND_FORWARD ? 1 : -1);
      const int coord = mspGetQStart(msp) + (direction * msp->phase);
      
      /* Find the reading frame that this coord belongs in. This is the same as the base number within
       * reading frame 1. */
      int frame = UNSET_INT;
      const gboolean invertCoords = (mspGetRefStrand(msp) == BLXSTRAND_REVERSE);
      
      convertDnaIdxToDisplayIdx(coord, seqType, 1, numFrames, invertCoords, refSeqRange, &frame);
      
      if (frame > 0)
	{
	  if (msp->qFrame > 0 && mspIsExon(msp))
	    {
	      /* Frame is already set. This means we have the old exblx file format, which does not provide phase for exons
	       * so we must use the reading frame that was passed in the file. However, this used to assume that the first base
	       * in the ref seq was base 1 in reading frame 1, which is not always the case (we now calculate the correct reading frame
	       * based on mod3 of the coord). Therefore, we must offset the given frame if the first coord is not base1 in frame1. */
	      const int startCoord = (msp->qStrand == BLXSTRAND_REVERSE ? refSeqRange->max : refSeqRange->min);
	      const int offset = getOffsetToCodonStart(startCoord, numFrames, msp->qStrand);
	      msp->qFrame += offset;
	      
	      if (msp->qFrame > numFrames)
		{
		  msp->qFrame -= numFrames;
		}
	      
	      if (msp->qFrame != frame && seqType == BLXSEQ_PEPTIDE)
		{
		  g_warning("MSP '%s' (q=%d-%d; s=%d-%d) has reading frame '%d' but calculated frame was '%d'\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max, msp->qFrame, frame);
		}
	    }
	  else
	    {
	      /* We either have the new file format which provides phase or it's not an exon so phase 
	       * is not applicable, so we can trust the calculated value. */
	      msp->qFrame = frame;
	    }
	}
      
      if (msp->qFrame == UNSET_INT)
	{
	  g_warning("Reading frame could not be calculated for MSP '%s' (q=%d-%d; s=%d-%d) - setting to 1.\n", mspGetSName(msp), msp->qRange.min, msp->qRange.max, msp->sRange.min, msp->sRange.max);
	  msp->qFrame = 1;
	}
    }
}


/* Should be called after all parsed data has been added to a BlxSequence. Calculates summary
 * data and the introns etc. */
void finaliseBlxSequences(GArray* featureLists[], 
			  MSP **mspList, 
			  GList **seqList, 
                          GList *columnList,
			  const int offset,
			  const BlxSeqType seqType, 
			  const int numFrames,
			  const IntRange const *refSeqRange,
			  const gboolean calcFrame,
                          const gboolean linkFeatures)
{
  /* Loop through all MSPs and adjust their coords by the offest, then calculate their reading
   * frame. Also find the last MSP in the list. */
  MSP *msp = *mspList;
  MSP *lastMsp = msp;

  while (msp)
    {
      adjustMspCoordsByOffset(msp, offset);
    
      if (calcFrame)
	calcReadingFrame(msp, seqType, numFrames, refSeqRange);
    
      msp = msp->next;
      if (msp)
        lastMsp = msp;
    }
  
  /* Loop through all BlxSequences */
  GList *seqItem = *seqList;
  
  for ( ; seqItem; seqItem = seqItem->next)
    {
      /* So far we only have the forward strand version of each sequence. We must complement any 
       * that need the reverse strand */
      BlxSequence *blxSeq = (BlxSequence*)(seqItem->data);

      if (blxSeq && 
          blxSeq->strand == BLXSTRAND_REVERSE && 
          (blxSeq->type == BLXSEQUENCE_MATCH || (blxSeq->type == BLXSEQUENCE_READ_PAIR && !ALL_READ_PAIRS_FWD)))
        {
          char *sequence = g_strdup(blxSequenceGetSequence(blxSeq));
          blxComplement(sequence);
          blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQUENCE, sequence);
        }
      
      findSequenceExtents(blxSeq);
      constructTranscriptData(blxSeq, featureLists, &lastMsp, mspList, seqList, columnList, linkFeatures);
    }

  /* Sort msp arrays by start coord (only applicable to msp types that
   * appear in the detail-view because the order is only applicable when
   * filtering detail-view rows) */
  int typeId = 0;
  for ( ; typeId < BLXMSP_NUM_TYPES; ++typeId)
    {
      if (typeShownInDetailView(typeId))
        g_array_sort(featureLists[typeId], compareFuncMspArray);
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



/* Determine the type of BlxSequence that a particular MSP type belongs to */
static BlxSequenceType getBlxSequenceTypeForMsp(const BlxMspType mspType)
{
  BlxSequenceType result = BLXSEQUENCE_UNSET;
  
  if (mspType == BLXMSP_MATCH)
    {
      result = BLXSEQUENCE_MATCH;
    }
  else if (typeIsExon(mspType) || typeIsIntron(mspType) || mspType == BLXMSP_TRANSCRIPT)
    {
      result = BLXSEQUENCE_TRANSCRIPT;
    }
  else if (typeIsVariation(mspType))
    {
      result = BLXSEQUENCE_VARIATION;
    }
  else if (typeIsShortRead(mspType))
    {
      result = BLXSEQUENCE_READ_PAIR;
    }
  else if (mspType == BLXMSP_REGION)
    {
      result = BLXSEQUENCE_REGION;
    }

  return result;
}


/* Utility to get the unique key for a given text string and strand */
static GQuark getLookupKey(const char *text, const BlxStrand strand)
{
  char *keyStr = g_strdup_printf("%s%c", text, (strand == BLXSTRAND_FORWARD ? '+' : '-'));
  GQuark key = g_quark_from_string(keyStr);
  g_free(keyStr);
  
  return key;
}


/* Utility to find a blxsequence with the given name/id/strand in the 
 * given hash table. Returns null if it is not there. */
static BlxSequence* findBlxSequence(GHashTable *lookupTable,
                                    const char *name, 
                                    const char *idTag,
                                    const BlxStrand strand,
                                    const gboolean linkFeaturesByName)
{
  BlxSequence *result = NULL;

  if (idTag)
    {
      /* We compare on the id tag and strand, so combine these into a single
       * string and convert it to a quark for quicker comparisions.
       * (This is the key for the hash table.) */
      GQuark key = getLookupKey(idTag, strand);
      result = (BlxSequence*)g_hash_table_lookup(lookupTable, GINT_TO_POINTER(key));
    }

  if (!result && name && linkFeaturesByName)
    {
      /* No id tag is given but we are asked to link features with the same
       * name, so do the comparison using name and strand */
      GQuark key = getLookupKey(name, strand);
      result = (BlxSequence*)g_hash_table_lookup(lookupTable, GINT_TO_POINTER(key));
    }

  return result;
}


/* Sort two columns by index */
gint columnIdxCompareFunc(gconstpointer a, gconstpointer b)
{
  BlxColumnInfo *col1 = (BlxColumnInfo*)a;
  BlxColumnInfo *col2 = (BlxColumnInfo*)b;
  return col1->columnIdx - col2->columnIdx;
}

/* Sort two columns by ID */
static gint columnIdCompareFunc(gconstpointer a, gconstpointer b)
{
  BlxColumnInfo *col1 = (BlxColumnInfo*)a;
  BlxColumnInfo *col2 = (BlxColumnInfo*)b;
  return col1->columnId - col2->columnId;
}


/* Create a new sequence with the given name and/or ID tag */
static BlxSequence* createBlxSequence(const char *name,
                                      const char *idTag, 
                                      const BlxStrand strand,
                                      BlxDataType *dataType,
                                      const char *source,
                                      GList *columnList)
{
  BlxSequence *seq = createEmptyBlxSequence();

  /* Set the hard-coded values */
  seq->strand = strand;
  seq->dataType = dataType;

  if (idTag)
    seq->idTag = g_strdup(idTag);

  /* Allocate the list of values and initialise the value types
   * to the relevant type for that column */
  //const int numColumns = g_list_length(columnList);
  seq->values = g_array_new(FALSE, FALSE, sizeof(GValue));

  /* Sort the columns by column ID (NOT column index) so that we
   * can easily index on the ID when looking up values in the array */
  columnList = g_list_sort(columnList, columnIdCompareFunc);
  
  GList *item = columnList;
  for ( ; item; item = item->next)
    {
      /* Note that we index on columnId, not columnIdx. This is so that
       * we can quickly index on the columnId enum. The columnIdx represents
       * the display order, which is not relevant here. */
      BlxColumnInfo * columnInfo = (BlxColumnInfo*)(item->data);
      GType type = columnInfo->type;
      
      /* Bit of a hack; the sequence column type is pointer because we
       * pass a pointer to the msp to the cell renderer, but here we want to
       * store the actual sequence string */
      if (columnInfo->columnId == BLXCOL_SEQUENCE)
        type = G_TYPE_STRING;
      
      GValue value = {0};
      g_value_init(&value, type);

      g_array_append_val(seq->values, value);
    }

  /* Set the given name and source */
  blxSequenceSetValueFromString(seq, BLXCOL_SEQNAME, name);
  blxSequenceSetValueFromString(seq, BLXCOL_SOURCE, source);

  /* Make sure we change back to the original sort order (sorted by index) */
  columnList = g_list_sort(columnList, columnIdxCompareFunc);

  return seq;
}


/* Add or create a BlxSequence struct, creating the BlxSequence if one does not
 * already exist for the MSP's sequence name. Seperate BlxSequence structs are created
 * for the forward and reverse strands of the same sequence. The passed-in sequence 
 * should always be forwards, and we reverse complement it here if we need the 
 * reverse strand. Returns the new BlxSequence */
BlxSequence* addBlxSequence(const char *name, 
			    const char *idTag, 
			    BlxStrand strand,
			    BlxDataType *dataType,
                            const char *source,
			    GList **seqList, 
                            GList *columnList,
			    char *sequence, 
			    MSP *msp, 
                            const gboolean linkFeaturesByName,
			    GError **error)
{
  BlxSequence *blxSeq = NULL;
  
  /* Put all blxseqs in a hash table indexed on a quark of the name
   * so that we can quickly check if the same one already exists */
  static GHashTable *lookupTable =  NULL;

  if (!lookupTable)
    lookupTable = g_hash_table_new(g_direct_hash, g_direct_equal);
    
  if (name || idTag)
    {
      /* If this is an exon or intron the match strand is not applicable. The exon should 
       * be in the same direction as the ref seq, so use the ref seq strand. */
      if (msp && (mspIsExon(msp) || mspIsIntron(msp)))
        {
          strand = msp->qStrand;
        }
    
      /* See if this sequence already exists. This matches on name (if linkFeaturesByName is
       * true) or on tag, and strand. */
      blxSeq = findBlxSequence(lookupTable, name, idTag, strand, linkFeaturesByName);
      
      if (!blxSeq)
        {
          /* Create a new BlxSequence, and take ownership of the passed in sequence (if any) */
          blxSeq = createBlxSequence(name, idTag, strand, dataType, source, columnList);
          
          /* Add it to the return sequence list */
          *seqList = g_list_prepend(*seqList, blxSeq);

          /* Add an entry to the lookup table (add an entry for both id and name, if given,
           * because the next feature may have only one or the other set.) */
          if (idTag)
            g_hash_table_insert(lookupTable, GINT_TO_POINTER(getLookupKey(idTag, strand)), blxSeq);

          if (name && linkFeaturesByName)
            g_hash_table_insert(lookupTable, GINT_TO_POINTER(getLookupKey(name, strand)), blxSeq);
        }
      else
        {
          if (dataType && !blxSeq->dataType)
            blxSeq->dataType = dataType;
          else if (dataType && blxSeq->dataType != dataType)
            g_warning("Duplicate sequences have different data types [name=%s, ID=%s, strand=%d, orig type=%s, new type=%s].\n", name, idTag, strand, g_quark_to_string(blxSeq->dataType->name), g_quark_to_string(dataType->name));
          
          const char *oldSource = blxSequenceGetSource(blxSeq);

          if (source && !oldSource)
            blxSequenceSetValueFromString(blxSeq, BLXCOL_SOURCE, source);
          else if (source && !stringsEqual(oldSource, source, FALSE))
            g_warning("Duplicate sequences have different sources [name=%s, ID=%s, strand=%d, orig source=%s, new source=%s].\n", name, idTag, strand, oldSource, source);
        }
      
      if (name && !blxSequenceGetName(blxSeq))
	{
	  /* It's possible that the BlxSequence was created without a name if we found an
	   * unnamed child exon before we found the parent transcript, so set the name if we have it. */
          blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQNAME, name);

          /* Add an entry to the lookup table keyed on name, now that we know it */
          if (linkFeaturesByName)
            g_hash_table_insert(lookupTable, GINT_TO_POINTER(getLookupKey(name, strand)), blxSeq);
	}
      
      if (msp)
        {
          /* Add the MSP to the BlxSequence's list. Keep it sorted by position. */
          blxSeq->mspList = g_list_insert_sorted(blxSeq->mspList, msp, compareFuncMspPos);
          msp->sSequence = blxSeq;
          
          if (blxSeq->type == BLXSEQUENCE_UNSET)
            {
              blxSeq->type = getBlxSequenceTypeForMsp(msp->type);
            }
          else if (blxSeq->type != getBlxSequenceTypeForMsp(msp->type))
            {
              g_warning("Adding MSP of type %d to parent of type %d (expected parent type to be %d)\n", 
			msp->type, blxSeq->type, getBlxSequenceTypeForMsp(msp->type));
            }
        }
      
      /* Add the sequence data */
      addBlxSequenceData(blxSeq, sequence, error);
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

  const char *oldSequence = blxSequenceGetSequence(blxSeq);
  
  if (blxSeq && blxSequenceRequiresSeqData(blxSeq))
    {
      if (!oldSequence)
        {
          /* Sequence does not yet exist, so add it */
          blxSequenceSetValueFromString(blxSeq, BLXCOL_SEQUENCE, sequence);
          sequenceUsed = TRUE;
        }
      else if (error && *error)
        {
          /* Sequence already exists. Validate that it's the same as the existing one. */
          if (!stringsEqual(sequence, oldSequence, FALSE))
            {
              g_set_error(error, BLX_ERROR, BLX_ERROR_SEQ_DATA_MISMATCH, "Sequence data for '%s' does not match previously-found data.\n", blxSequenceGetName(blxSeq));
            }
        }
    }      
  
  if (!sequenceUsed)
    {
      g_free(sequence);
    }
}

