/*
 *  utilities.c
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */

#include "SeqTools/utilities.h"

static CallbackData*	  widgetGetCallbackData(GtkWidget *widget);

/* Functions to get/set/destroy a GdkDrawable in a widget property, if a different drawable
 * than the window is required to be drawn to. The main purpose of this is for printing
 * (and only widgets that have this drawable set are printed) */
static void onDestroyCustomWidget(GtkWidget *widget)
{
  GdkDrawable *drawable = widgetGetDrawable(widget);
  if (drawable)
    {
      g_object_unref(drawable);
      drawable = NULL;
      g_object_set_data(G_OBJECT(widget), "drawable", NULL);
    }
  
  CallbackData *callbackData = widgetGetCallbackData(widget);
  if (callbackData)
    {
      g_free(callbackData);
      g_object_set_data(G_OBJECT(widget), "callbackData", NULL);
    }
}

GdkDrawable* widgetGetDrawable(GtkWidget *widget)
{
  return widget ? (GdkDrawable*)(g_object_get_data(G_OBJECT(widget), "drawable")) : NULL;
}

void widgetSetDrawable(GtkWidget *widget, GdkDrawable *drawable)
{
  if (widget)
    { 
      /* Delete the old one first, if there is one */
      GdkDrawable *oldDrawable = widgetGetDrawable(widget);
      if (oldDrawable)
	{
	  g_object_unref(oldDrawable);
	  oldDrawable = NULL;
	}

      g_object_set_data(G_OBJECT(widget), "drawable", drawable);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyCustomWidget), NULL);
    }
}


/* Call this function to clear the cached drawable for the widget. This means that the
 * next time that expose is called, the bitmap will be redrawn from scratch. */
void widgetClearCachedDrawable(GtkWidget *widget)
{  
  widgetSetDrawable(widget, NULL);
}


/* Functions to get/set a boolean property that is set to true if this widget has
 * been hidden by the user. (We can't just use gtk_widget_hide because we do a 
 * show_all when we remove and re-add widgets to the parent when we toggle strands. 
 * Note that if the property is not set the result is FALSE, so by default widgets
 * are not hidden. */
gboolean widgetGetHidden(GtkWidget *widget)
{
  gboolean result = FALSE;
  
  if (widget)
    {	
      result = (gboolean)GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "hidden"));
    }
  
  return result;
}

void widgetSetHidden(GtkWidget *widget, const gboolean hidden)
{
  if (widget)
    {
      /* Set the property */
      g_object_set_data(G_OBJECT(widget), "hidden", GINT_TO_POINTER((gint)hidden));
      
      /* Show/hide the widget */
      if (hidden)
	{
	  gtk_widget_hide_all(widget);
	}
      else
	{
	  gtk_widget_show_all(widget);
	}	
    }
}


/* Hides the given widget if it has been set as hidden by the user. Recurses
 * over all children. */
void hideUserHiddenWidget(GtkWidget *widget, gpointer callbackData)
{
  if (widgetGetHidden(widget))
    {
      gtk_widget_hide_all(widget);
    }
  
  if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), hideUserHiddenWidget, NULL);
    }
}


/* Create a label with the given properties. If 'showWhenPrinting' is false 
 * this label will not appear when printing */
GtkWidget* createLabel(const char *text, 
		       const gdouble xalign,
		       const gdouble yalign,
		       const gboolean enableCopyPaste,
		       const gboolean showWhenPrinting)
{
  GtkWidget *label = NULL;
  
  if (text)
    {
      label = gtk_widget_new(GTK_TYPE_LABEL, 
                             "label", text, 
                             "xalign", xalign,
                             "yalign", yalign,
                             "selectable", enableCopyPaste,
                             NULL);
    }
  else
    {
      label = gtk_widget_new(GTK_TYPE_LABEL, 
                             "xalign", xalign, 
                             "yalign", yalign, 
                             "selectable", enableCopyPaste,
                             NULL);
    }
  
  if (label && showWhenPrinting)
    {
      /* Connect to the expose event handler that will create the drawable object required for printing */
      g_signal_connect(G_OBJECT(label), "expose-event", G_CALLBACK(onExposePrintableLabel), NULL);
    }

  return label;
}



/* Expose-event handler for labels that are required to be shown during printing. */
gboolean onExposePrintableLabel(GtkWidget *label, GdkEventExpose *event, gpointer callbackData)
{
  if (!label || !GTK_IS_LABEL(label))
    {
      g_critical("Invalid widget type when printing label [%p]\n", label);
    }
  
  /* Only widgets that have a pixmap set will be shown in print output */
  GtkWidget *parent = gtk_widget_get_parent(label);
  GdkDrawable *drawable =  gdk_pixmap_new(parent->window, label->allocation.width, label->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(label, drawable);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GtkStyle *style = gtk_widget_get_style(label);
  GdkColor *bgColor = &style->bg[GTK_STATE_NORMAL];
  gdk_gc_set_foreground(gc, bgColor);
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, label->allocation.width, label->allocation.height);

  PangoLayout *layout = gtk_label_get_layout(GTK_LABEL(label));
  if (layout)
    {
      gtk_paint_layout(label->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, label, NULL, 0, 0, layout);
    }
  
  return FALSE; /* let the default handler continue */
}


/* Utility to return the length of the given range */
int getRangeLength(const IntRange const *range)
{
  return (range->max - range->min + 1);
}


/* Utility to return the centre value of the given range (rounded down if an odd number) */
int getRangeCentre(const IntRange const *range)
{
  return range->min + ((range->max - range->min) / 2);
}


/* Simple utility to determine whether the given value is within the given range.
 * Returns false if the given value is an unset int or the given range is null */
gboolean valueWithinRange(const int value, const IntRange const *range)
{
  return (value != UNSET_INT && range != NULL && value >= range->min && value <= range->max);
}


/* Utility to bounds-limit the given value to within the given range */
void boundsLimitValue(int *value, const IntRange const *range)
{
  if (*value < range->min)
    {
      *value = range->min;
    }
  
  if (*value > range->max)
    {
      *value = range->max;
    }
}


/* Utility to calculate how many digits are in an integer (including the '-'
 * sign if the integer is negative) */
int numDigitsInInt(int val)
{
  int count = 0;

  if (val == 0)
    {
      count = 1;
    }
  else
    {
      if (val < 0)
	{
	  /* Add one for the '-' sign, then treat it as a positive number */
	  ++count;
	  val *= -1;
	}

      while (val > 0)
	{
	  ++count;
	  val /= 10;
	}
    }
      
  return count;
}


/* Create a GdkColor from one of our internally-defined colors, e.g. BLX_YELLOW.
 * At the moment this is just a wrapper around gdk_color_parse. Returns true if ok, 
 * false if there was an error. */
gboolean parseBlxColor(const char *color, GdkColor *result)
{
  gboolean ok = gdk_color_parse(color, result);

  if (ok)
    {
      gboolean failures[1];
      gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), result, 1, TRUE, TRUE, failures);
      
      if (numFailures > 0)
	ok = FALSE;
    }

  if (!ok)
    {
      g_warning("Error parsing color string '%s'\n", color);
    }
  
  return ok;
}


/* Increase the brightness of a color by the given factor */
static void increaseColorBrightness(GdkColor *origColor, const double factor, GdkColor *result)
{
  result->pixel = 0;
  const int maxRgb = 65535;
  
//  int maxVal = origColor->red > origColor->green ? origColor->red : origColor->green;
//  
//  if (origColor->blue > maxVal)
//    maxVal = origColor->blue;
//  
//  int offset = maxVal * (factor - 1);
//
//  result->red	= ((int)origColor->red + offset > maxRgb)   ? maxRgb : (origColor->red + offset);
//  result->green = ((int)origColor->green + offset > maxRgb) ? maxRgb : (origColor->green + offset);
//  result->blue	= ((int)origColor->blue + offset > maxRgb)  ? maxRgb : (origColor->blue + offset);

  int newRed = origColor->red * factor;
  int newGreen = origColor->green * factor;
  int newBlue = origColor->blue * factor;
  
  gboolean bustMax = FALSE;

  if (newRed > maxRgb)
    {
      bustMax = TRUE;
      newRed = maxRgb;
    }

  if (newGreen > maxRgb)
    {
      bustMax = TRUE;
      newGreen = maxRgb;
    }

  if (newBlue > maxRgb)
    {
      bustMax = TRUE;
      newBlue = maxRgb;
    }
  
  if (bustMax)
    {
      const int offset = maxRgb * (factor - 1);
      newRed = (origColor->red + offset > maxRgb) ? maxRgb : origColor->red + offset;
      newGreen = (origColor->green + offset > maxRgb) ? maxRgb : origColor->green + offset;
      newBlue = (origColor->blue + offset > maxRgb) ? maxRgb : origColor->blue + offset;
    }
  
  result->red = newRed;
  result->green = newGreen;
  result->blue = newBlue;
  
  gboolean failures[1];
  gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), result, 1, TRUE, TRUE, failures);
  
  if (numFailures > 0)
    {
      printf("Warning: error calculating highlight color; using original color instead (RGB=%d %d %d).", origColor->red, origColor->green, origColor->blue);
      result->pixel = origColor->pixel;
      result->red = origColor->red;
      result->green = origColor->green;
      result->blue = origColor->blue;
    }
}


/* Reduce the brightness of a color by the given factor (e.g. 0.5 for half brightness) */
static void reduceColorBrightness(GdkColor *origColor, const double factor, GdkColor *result)
{
  result->pixel = 0;
  result->red = origColor->red * factor;
  result->green = origColor->green * factor;
  result->blue = origColor->blue * factor;
  
  gboolean failures[1];
  gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), result, 1, TRUE, TRUE, failures);
  
  if (numFailures > 0)
    {
      printf("Warning: error calculating highlight color; using original color instead (RGB=%d %d %d).", origColor->red, origColor->green, origColor->blue);
      result->pixel = origColor->pixel;
      result->red = origColor->red;
      result->green = origColor->green;
      result->blue = origColor->blue;
    }
}


/* Increase/decrease color brightness by given factor */
void adjustColorBrightness(GdkColor *origColor, const double factor, GdkColor *result)
{
  if (factor > 1)
    increaseColorBrightness(origColor, factor, result);
  else if (factor < 1)
    reduceColorBrightness(origColor, factor, result);
}


/* Utility to take a gdk color and return a slightly darker shade of it, for higlighting
 * things of that color when they are selected. */
void getSelectionColor(GdkColor *origColor, GdkColor *result)
{
  const double factor = 0.8;
  adjustColorBrightness(origColor, factor, result);
}


/* Utility to take a gdk color and return a much darker shade of it, for use as a drop-shadow */
void getDropShadowColor(GdkColor *origColor, GdkColor *result)
{
  const double factor = 0.3;
  adjustColorBrightness(origColor, factor, result);
}


gboolean mspIsExon(const MSP const *msp)
{
  return (msp && 
          (msp->type == BLXMSP_EXON_CDS ||
           msp->type == BLXMSP_EXON_UTR || 
           msp->type == BLXMSP_EXON_UNK));
}

/* Determine whether the given MSP is in a coding region or untranslated region. For 
 * exons, this is determined by the exon type. For introns, we have to look at the
 * adjacent exons to determine whether to show them as CDS or UTR - only show it as
 * CDS if there is a CDS exon on both sides of the intron. */
gboolean mspIsCds(const MSP const *msp, const BlxSequence *blxSeq)
{
  gboolean result = FALSE;
  
  if (mspIsExon(msp))
    {
      result = msp && msp->type == BLXMSP_EXON_CDS;
    }
  else if (mspIsIntron(msp))
    {
      /* Find the nearest exons before and after this MSP */
      const MSP *prevExon = NULL;
      const MSP *nextExon = NULL;
    
      GList *mspItem = blxSeq->mspList;
      for ( ; mspItem; mspItem = mspItem->next)
	{
	  const MSP *curMsp = (const MSP *)(mspItem->data);
	
	  if (mspIsExon(curMsp))
	    {
	      const int curOffset = curMsp->qstart - msp->qstart;
	    
	      if (curOffset < 0 && (!prevExon || curOffset > prevExon->qstart - msp->qstart))
		{
		  /* Current MSP is before our MSP and is the smallest offset so far */
		  prevExon = curMsp;
		}
	      else if (curOffset > 0 && (!nextExon || curOffset < nextExon->qstart - msp->qstart))
		{
		  /* Current MSP is after our MSP and is the smallest offset so far */
		  nextExon = curMsp;
		}
	    }
	}

      /* If there isn't an exon either side, we don't know what to display. Default to UTR. */
      result = prevExon || nextExon;
      
      /* For the exon(s) that we found, they must both be CDS if the intron is to be CDS. If
       * we only have an exon on one side, and it is CDS, then the intron is CDS. */
      result &= (!prevExon || prevExon->type == BLXMSP_EXON_CDS) &&
		(!nextExon || nextExon->type == BLXMSP_EXON_CDS);
    }
  
  return result;
}

gboolean mspIsIntron(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_INTRON);
}

gboolean mspIsSnp(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_SNP);
}

gboolean mspIsBlastMatch(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_MATCH);
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


/* Sorts the given values so that val1 is less than val2 if forwards is true,
 * or the reverse if forwards is false. */
void sortValues(int *val1, int *val2, gboolean forwards)
{
  if ((forwards && *val1 > *val2) || (!forwards && *val1 < *val2))
    {
      int temp = *val1;
      *val1 = *val2;
      *val2 = temp;
    }
}


/* Returns the upper and lower extents of the query and subject sequence ranges in 
 * the given MSP. Any of the return values can be NULL if they are not required. */
void getMspRangeExtents(const MSP *msp, int *qSeqMin, int *qSeqMax, int *sSeqMin, int *sSeqMax)
{
  if (qSeqMin)
    *qSeqMin = msp->qstart < msp->qend ? msp->qstart : msp->qend;
  
  if (qSeqMax)
    *qSeqMax = msp->qstart < msp->qend ? msp->qend : msp->qstart;
  
  if (sSeqMin)
    *sSeqMin = msp->sstart < msp->send ? msp->sstart : msp->send;
  
  if (sSeqMax)
    *sSeqMax = msp->sstart < msp->send ? msp->send : msp->sstart;
}


/* Returns the upper and lower extents of the query and subject sequence ranges in
 * the given gap range. Any of the return values can be passed as NULL if they are not required. */
void getCoordRangeExtents(CoordRange *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax)
{
  if (qRangeMin)
    *qRangeMin = min(range->qStart, range->qEnd);
  
  if (qRangeMax)
    *qRangeMax = max(range->qStart, range->qEnd);
  
  if (sRangeMin)
    *sRangeMin = min(range->sStart, range->sEnd);
  
  if (sRangeMax)
    *sRangeMax = max(range->sStart, range->sEnd);
}


char convertBaseToCorrectCase(const char charToConvert, const BlxSeqType seqType)
{
  char result = (seqType == BLXSEQ_PEPTIDE) ? toupper(charToConvert) : tolower(charToConvert);
  return result;
}


/* Returns the base at the given index in the reference sequence. If
 * 'complement' is true, the result is complemented. */
char getRefSeqBase(char *refSeq,
		   const int qIdx, 
		   const gboolean complement, 
		   const IntRange const *refSeqRange,
		   const BlxSeqType seqType)
{
  char result = ' ';
  
  if (qIdx >= refSeqRange->min && qIdx <= refSeqRange->max)
    {
      char base = refSeq[qIdx - refSeqRange->min];
      base = convertBaseToCorrectCase(base, seqType);
      
      if (!complement)
	{
	  result = base;
	}
      else
	{
	  switch (base)
	  {
	    case 'c':
	      result = 'g';
	      break;
	    case 'g':
	      result = 'c';
	      break;
	    case 'a':
	      result = 't';
	      break;
	    case 't':
	      result = 'a';
	      break;
	    default:
	      result = ' ';
	      break;
	  }
	}
    }
  
  return result;
}


/* Given an index into the displayed sequence, a reading frame, and the base number within that
 * reading frame, return the index into the DNA sequence that will give the equivalent DNA base.
 * If the display sequence is a peptide sequence, it will convert the coord to a DNA coord. If the
 * display is reversed, the display coord will be inverted. */
int convertDisplayIdxToDnaIdx(const int displayIdx, 
			      const BlxSeqType displaySeqType,
			      const int frame, 
			      const int baseNum, 
			      const int numFrames,
			      const gboolean displayRev,
			      const IntRange const *refSeqRange)
{
  int dnaIdx = displayIdx;
  
  if (displayIdx != UNSET_INT)
    {
      if (displaySeqType == BLXSEQ_PEPTIDE)
	{
	  /* Convert the input peptide coord to a dna coord */
	  dnaIdx = (displayIdx * numFrames) - numFrames + frame + baseNum - 1;

	  /* Convert from 1-based coord to real coords */
	  dnaIdx += refSeqRange->min - 1;
	}
      
      /* If the display is reversed, we need to invert the result. (For example, if the 
       * result is index '2' out of the range '12345', then we convert it to '4', which is the
       * equivalent position in the range '54321'. */
      if (displayRev)
	{
	  dnaIdx = refSeqRange->max - dnaIdx + refSeqRange->min;
	}
    }

  return dnaIdx;
}


/* Given an index into a dna sequence and the reading frame, find the index into 
 * the display sequence. Converts the index to a peptide index if displaying
 * peptide sequences (Also returns the base number of the DNA base within the
 * reading frame in this case (i.e. whether it's the 1st, 2nd or 3rd out of the 
 * triplet). */
int convertDnaIdxToDisplayIdx(const int dnaIdx, 
			      const BlxSeqType displaySeqType,
			      const int frame, 
			      const int numFrames, 
			      const gboolean displayRev,
			      const IntRange const *dnaIdxRange,
			      int *baseNum)
{
  int displayIdx = dnaIdx;
  
  if (dnaIdx != UNSET_INT)
    {
      /* If the display is reversed (i.e. showing increasing values from right-to-left),
       * invert the index (i.e. as if it were the normal left-to-right index for this
       * same position - for example, if the index is '4' out of the range '54321', convert
       * it to '2', which is the equivalent position in the range '12345'. */
      if (displayRev)
	{
	  displayIdx = dnaIdxRange->max - dnaIdx + dnaIdxRange->min;
	}
      
      if (displaySeqType == BLXSEQ_PEPTIDE)
	{
	  /* Must convert to 1-based coords to get correct reading frame (because that's
	   * the coordinate system used in blixem's input files). */
	  displayIdx -= dnaIdxRange->min - 1;

	  /* We're displaying peptides, so convert the dna coord to a peptide coord */
	  gdouble fraction = ((gdouble)(displayIdx - frame + 1) / (gdouble)numFrames) ;
	  displayIdx = ceil(fraction);
      
	  /* Find the base number of this DNA coord within the codon, if requested */
	  if (baseNum)
	    {
	      *baseNum = roundNearest((fraction - (int)fraction) * numFrames);
	  
	      if (*baseNum < 1)
		{
		  *baseNum += numFrames;
		}
	    }
	}
      else if (baseNum)
	{
	  /* For dna sequences, we only have one reading frame */
	  *baseNum = 1;
	}
    }
  
  return displayIdx;
}


/* Get the start coord in the given display range, reversed as necessary if 'reverse' is true.
 * Result is a coord in the DNA sequence - converts as necessary if the display range is in terms
 * of peptide coords */
int getStartDnaCoord(const IntRange const *displayRange, 
		     const int frame,
		     const BlxSeqType displaySeqType, 
		     const gboolean displayRev, 
		     const int numFrames,
		     const IntRange const *refSeqRange)
{
  int result = displayRange->min;
  
  /* Convert the display coord to coords into the ref seq, which is a DNA sequence. We want
   * the first base in the codon, if this is a peptide sequence. */
  const int baseNum = 1;
  result = convertDisplayIdxToDnaIdx(result, displaySeqType, frame, baseNum, numFrames, displayRev, refSeqRange);
  
  return result;
}


/* Get the end coord in the given display range, reversed as necessary if 'reverse' is true.
 * Result is a coord in the DNA sequence - converts as necessary if the display range is in terms
 * of peptide coords */
int getEndDnaCoord(const IntRange const *displayRange, 
		   const int frame,
		   const BlxSeqType displaySeqType, 
		   const gboolean displayRev, 
		   const int numFrames,
		   const IntRange const *refSeqRange)
{
  int result = displayRange->max;
  
  /* Convert the display coord to coords into the ref seq, which is a DNA sequence. We want
   * the last base in the codon, if this is a peptide sequence. */
  const int baseNum = numFrames;
  result = convertDisplayIdxToDnaIdx(result, displaySeqType, frame, baseNum, numFrames, displayRev, refSeqRange);
  
  return result;
}


/* Given a base index on the query sequence, find the corresonding base 
 * in subject sequence. The return value is always UNSET_INT if there is not
 * a corresponding base at this position. However, in this case the start/end
 * of the sequence (or nearest gap) is returned in the nearestIdx argument. If
 * the base exists then the return idx and nearestIdx will be the same. 
 * nearestIdx can be null if not required. */
int gapCoord(const MSP *msp, 
	     const int qIdx, 
	     const int numFrames, 
	     const BlxStrand strand, 
	     const gboolean displayRev)
{
  int result = UNSET_INT;
  
  if (mspIsBlastMatch(msp) || mspIsExon(msp))
    {
      int qSeqMin = UNSET_INT, qSeqMax = UNSET_INT, sSeqMin = UNSET_INT, sSeqMax = UNSET_INT;
      getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
      
      const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
      const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
      const gboolean sameDirection = (qForward == sForward);
      
      if (!msp->gaps || g_slist_length(msp->gaps) < 1)
	{
	  /* If strands are in the same direction, find the offset from qSeqMin and add it to sSeqMin.
	   * If strands are in opposite directions, find the offset from qSeqMin and subtract it from sSeqMax. */
	  int offset = (qIdx - qSeqMin)/numFrames ;
	  result = (sameDirection) ? sSeqMin + offset : sSeqMax - offset ;
	  
	  if (result < sSeqMin)
	    {
	      result = UNSET_INT;
	    }
	  else if (result > sSeqMax)
	    {
	      result = UNSET_INT;
	    }
	}
      else
	{
	  /* Look to see if x lies inside one of the "gaps" ranges. */
          GSList *rangeItem = msp->gaps;

	  for ( ; rangeItem ; rangeItem = rangeItem->next)
	    {
	      CoordRange *curRange = (CoordRange*)(rangeItem->data);
	      
	      int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
	      getCoordRangeExtents(curRange, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
	      
	      /* We've "found" the value if it's in or before this range. Note that the
	       * the range values are in decreasing order if the q strand is reversed. */
	      gboolean found = qForward ? qIdx <= qRangeMax : qIdx >= qRangeMin;
	      
	      if (found)
		{
		  gboolean inRange = (qForward ? qIdx >= qRangeMin : qIdx <= qRangeMax);
		  
		  if (inRange)
		    {
		      /* It's inside this range. Calculate the actual index. */
		      int offset = (qIdx - qRangeMin) / numFrames;
		      result = sameDirection ? sRangeMin + offset : sRangeMax - offset;
		    }

		  break;
		}
	    }
	}
    }
  
  return result;
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
  return msp->sSequence ? msp->sSequence->sequence : NULL;
}


/* Get the sequence name for an MSP. This removes the postfixed 'x' or 'i' from
 * exon and intron names; otherwise it just returns a copy of the name in the MSP.
 * The result should be free'd with g_free. */
char* mspGetSeqName(const MSP *msp)
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


/* Return a char representation of a strand, i.e. "+" for forward strand, "-" for
 * reverse, or "." for none. */
char getStrandAsChar(const BlxStrand strand)
{
  char result = '.';
  
  if (strand == BLXSTRAND_FORWARD)
    {
      result = '+';
    }
  else if (strand == BLXSTRAND_REVERSE)
    {
      result = '-';
    }
  
  return result;
}


/* Utility function to extract the variant name from a 
 * long sequence name (of the form LL:LLdddddd.d, where L means letter and d
 * means digit). The result is a pointer into the original string, so should 
 * not be free'd. */
const char* getSeqVariantName(const char *longName)
{
  /* Ignore the text before the colon */
  char *cutPoint = strchr(longName, ':');
  
  if (cutPoint)
    ++cutPoint;
  
  const char *result = cutPoint ? cutPoint : longName;
  
  return result;
}

/* Extract the short version of a sequence name. If the original is of the format
 * "Em:AB12345.1", this cuts off the bit before the ":" and the bit after the "."
 * The result is a new string that must be free'd with g_free. */
char* getSeqShortName(const char *longName)
{
  const int len = strlen(longName);

  char *result = g_malloc(sizeof(char) * len);
 
  /* This bool will be set to true when we have passed the ":". It means that we are
   * ready to start copying chars into the result. If there is no ":", set it to true
   * right at the beginning. */
  char *cutPoint = strchr(longName, ':');
  gboolean start = (cutPoint == NULL);
  
  int srcIdx = 0;
  int destIdx = 0;
  
  /* Copy contents after the : (or from beginning if : doesn't exist), and up to the . (if exists) */
  for ( ; srcIdx < len; ++srcIdx)
    {
      if (start)
	{
	  result[destIdx] = longName[srcIdx];
	  ++destIdx;
	}
      else if (longName[srcIdx] == ':')
	{
	  start = TRUE;
	}
      else if (longName[srcIdx] == '.')
	{
	  break;
	}
    }

  result[destIdx] = '\0';

  return result;
}


/* Return the full name of a BlxSequence (including prefix and variant) */
const char *sequenceGetFullName(const BlxSequence *seq)
{
  return seq->fullName;
}


/* Return the variant name of a BlxSequence (excludes prefix but includes variant) */
const char *sequenceGetVariantName(const BlxSequence *seq)
{
  return seq->variantName;
}


/* Return the display name of a BlxSequence (same as variant name for now) */
const char *sequenceGetDisplayName(const BlxSequence *seq)
{
  return seq->variantName;
}


/* Return the short name of a BlxSequence (excludes prefix and variant number) */
const char *sequenceGetShortName(const BlxSequence *seq)
{
  return seq->shortName;
}


/* Frees all memory used by a BlxSequence */
void destroyBlxSequence(BlxSequence *seq)
{
  if (seq)
    {
      g_free(seq->fullName);
      g_free(seq->variantName);
      g_free(seq->shortName);
      g_free(seq->sequence);
      
      g_free(seq);
    }
}


/* Utility to create a BlxSequence with the given name. */
BlxSequence* createEmptyBlxSequence(char *fullName)
{
  BlxSequence *seq = g_malloc(sizeof(BlxSequence));
  
  seq->fullName = g_strdup(fullName);
  
  /* The variant name: just cut off the prefix chars. We can use a pointer into
   * the original string. */
  seq->variantName = g_strdup(getSeqVariantName(fullName));
  
  /* The short name: cut off the prefix chars (before the ':') and the variant
   * number (after the '.'). Need to duplicate the string to change the end of it. */
  seq->shortName = g_strdup(seq->variantName);
  char *cutPoint = strchr(seq->shortName, '.');
  
  if (cutPoint)
    *cutPoint = '\0';

  seq->mspList = NULL;
  seq->sequence = NULL;
  
  return seq;
}


/* Round to nearest int (needed because there is no  round() function in ISO C90) */
int roundNearest(const double val)
{
  int truncated = (int)val;
  double remainder = val - truncated;
  int result = truncated;

  if (remainder > 0.5)
    {
      result = truncated + 1;
    }
  else if (remainder < -0.5)
    {
      result = truncated - 1;
    }
    
  return result;
}


/* Round to the nearest multiple of the given value */
int roundToValue(const int inputVal, const int roundTo)
{
  return roundNearest((double)inputVal / (double)roundTo) * roundTo;
}


/* Converts the given integer to a string. The result must be free'd with g_free */
char* convertIntToString(const int value)
{
  char result[numDigitsInInt(value) + 1];
  sprintf(result, "%d", value);
  return g_strdup(result);
}


/* Converts the given string to an integer. */
int convertStringToInt(const char *inputStr)
{
  return atoi(inputStr);
}


/***********************************************************
 *			   Dialogs			   *
 ***********************************************************/

/* Utility to create a scrollable text view from the given message text */
GtkWidget* createScrollableTextView(const char *messageText,
				    const gboolean wrapText,
				    PangoFontDescription *fontDesc,
				    int *height)
{
  /* Create a text buffer and copy the text in */
  GtkTextBuffer *textBuffer = gtk_text_buffer_new(NULL);
  gtk_text_buffer_set_text(textBuffer, messageText, -1);
  
  /* Create a text view to display the buffer */
  GtkWidget *textView = gtk_text_view_new();
  gtk_text_view_set_buffer(GTK_TEXT_VIEW(textView), textBuffer);
  gtk_text_view_set_editable(GTK_TEXT_VIEW(textView), FALSE);
  gtk_widget_modify_font(textView, fontDesc);
  
  if (wrapText)
    {
      gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(textView), TRUE);
    }
  
  /* Set the initial height based on the number of lines (but don't make it bigger than the parent window) */
  PangoContext *context = gtk_widget_get_pango_context(textView);
  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
  gint charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  
  /* Adjust height to include all the lines if possible (limit to the original height passed in thought) */
  const int calcHeight = (gtk_text_buffer_get_line_count(textBuffer) * charHeight) + 100; /* fudge to allow space for buttons */
  *height = min(*height, calcHeight);
  
  pango_font_metrics_unref(metrics);
  
  /* Put the text view in a scrollable window */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(scrollWin), textView);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  /* Return the outermost container */
  return scrollWin;
}


/* Utility to pop up a simple dialog with the given title and text, with just an "OK" button. */
void showMessageDialog(const char *title,  
		       const char *messageText,
		       GtkWidget *parent,
		       const int width,
		       const int maxHeight,
		       const gboolean wrapText,
		       PangoFontDescription *fontDesc)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  GTK_WINDOW(parent), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  int height = maxHeight;
  GtkWidget *child = createScrollableTextView(messageText, wrapText, fontDesc, &height);

  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);

  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
}


/* Get/set a CallbackData struct for this widget */
static CallbackData* widgetGetCallbackData(GtkWidget *widget)
{
  return widget ? (CallbackData*)(g_object_get_data(G_OBJECT(widget), "callbackData")) : NULL;
}


/* Set callback function and user-data for the given widget */
void widgetSetCallbackData(GtkWidget *widget, GtkCallback func, gpointer data)
{
  if (widget)
    { 
      CallbackData *callbackData = g_malloc(sizeof *callbackData);
      callbackData->func = func;
      callbackData->data = data;
      
      g_object_set_data(G_OBJECT(widget), "callbackData", callbackData);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyCustomWidget), NULL);
    }
}

/* Get the CallbackData struct for this widget, if it has one, and call it. */
void widgetCallCallback(GtkWidget *widget)
{
  CallbackData *callbackData = widgetGetCallbackData(widget);
  
  if (callbackData && callbackData->func)
    {
      callbackData->func(widget, callbackData->data);
    }
}

/* Call the stored CallbackData callbacks (if any exist) for this widget and
 * all of its children. */
void widgetCallAllCallbacks(GtkWidget *widget, gpointer data)
{
  if (widget && GTK_IS_WIDGET(widget))
    {
      widgetCallCallback(widget);
  
      if (GTK_IS_CONTAINER(widget))
	{
	  gtk_container_foreach(GTK_CONTAINER(widget), widgetCallAllCallbacks, NULL);
	}
    }
}

/* Generic callback to call all user-specified callbacks for all child widgets of
 * the given dialog if ACCEPT or APPLY responses received. Also closes the dialog 
 * if ACCEPT or REJECT responses received. */
void onResponseDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      widgetCallAllCallbacks(GTK_WIDGET(dialog), NULL);
      destroy = TRUE;
      break;
      
    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(GTK_WIDGET(dialog), NULL);
      destroy = FALSE;
      break;
      
    case GTK_RESPONSE_CANCEL:
    case GTK_RESPONSE_REJECT:
      destroy = TRUE;
      break;
      
    default:
      break;
  };
  
  if (destroy)
    {
      gtk_widget_destroy(GTK_WIDGET(dialog));
    }
}


/***********************************************************
 *			   Text parsing                    *
 ***********************************************************/

/* match to template with wildcards.   Authorized wildchars are * ? #
     ? represents any single char
     * represents any set of chars
   case-insensitive.   Example: *Nc*DE# fits abcaNchjDE23

   returns 0 if not found
           1 + pos of first sigificant match (i.e. not a *) if found
*/
int wildcardSearch(const char *textToSearch, const char *searchStr)
{
  char *textChar = (char*)textToSearch; /* to do: don't cast away const! */
  char *searchChar = (char*)searchStr;	/* to do: don't cast away const! */
  char *ts, *cs, *s = 0 ;
  int star=0;

  while (1)
    {
      switch(*searchChar)
	{
	case '\0':
	  {
	    if(!*textChar)
	      {
		return  ( s ? 1 + (s - textToSearch) : 1);
	      }
	    else if (!star)
	      {
		return 0;
	      }
	    else
	      {
		/* not success yet go back in template */
		searchChar = ts;
		textChar = cs + 1;
		
		if (ts == searchStr)
		  {
		    s = 0;
		  }
	      }
	    
	    break ;
	  }
	    
	case '?':
	  {
	    if (!*textChar)
	      {
		return 0;
	      }
	    
	    if(!s)
	      {
		s = textChar;
	      }
	    
	    searchChar++;
	    textChar++;
	    break;
	  }
	    
	case '*':
	  {
	    ts = searchChar;
	    
	    while( *searchChar == '?' || *searchChar == '*')
	      {
		searchChar++;
	      }
	    
	    if (!*searchChar)
	      {
		return s ? 1 + (s-textToSearch) : 1 ;
	      }
	    
	    while (freeupper(*textChar) != freeupper(*searchChar))
	      {
		if (*textChar)
		  textChar++;
		else
		  return 0;
	      }
	    
	    star = 1;
	    cs = textChar;
	    
	    if(!s)
	      {
		s = textChar;
	      }
	    
	    break;
	  }
	    
	default:
	  {
	    if (freeupper(*searchChar++) != freeupper(*textChar++))
	      {
		if(!star)
		  {
		    return 0;
		  }
		
		searchChar = ts;
		textChar = cs + 1;
		
		if(ts == searchStr)
		  {
		    s = 0;
		  }
	      }
	    else if(!s)
	      {
		s = textChar - 1 ;
	      }
	    
	    break;
	  }
	}
    }

  return 0;
}


/* Given a string, returns the string if is less than a fudge factor or returns
 * the string abbreviated in the form "xxxxx...yyyyy". The result should be 
 * free'd with g_free. */
gchar *abbreviateText(const char *inputStr, const int max_len)
{
  gchar *result = NULL;

  if (max_len > 0) 
    { 
      char abbrev[] = "<>";
      int inputStrLen = strlen(inputStr);
      
      if (inputStrLen <= max_len)
	{
	  result = g_strdup(inputStr);
	}
      else
	{
	  /* Get the first len/2 chars */
	  const int headChars = (max_len / 2) - 1 ;
	  char *headStr = headChars > 0 ? g_strndup(inputStr, headChars) : NULL;

	  const int tailChars = max_len - ((max_len / 2) + 1) ;
	  char *tailStr = g_strndup(inputStr + inputStrLen - tailChars, tailChars); /* trailing null fudged in here. */
	  
	  result = g_strconcat(headStr, abbrev, tailStr, NULL);
	  
	  g_free(headStr);
	  g_free(tailStr);
	}
    }
    
  return (result) ;
}


/* String comparison function with caseSensitive option. Returns true if the two 
 * strings are equal. The two strings must be null-terminated. */
gboolean stringsEqual(const char *str1, const char *str2, const gboolean caseSensitive)
{
  gboolean result = FALSE;
  
  if (caseSensitive)
    {
      result = !strcmp(str1, str2);
    }
  else
    {
      const int len1 = strlen(str1);
      result = (strlen(str2) == len1 && !g_ascii_strncasecmp(str1, str2, len1));
    }
  
  return result;
}


/* Parses a line of text containing a match description, which is of the form:
 * 
 *                       "\"name\" start end (length)"
 * or:			 "name start end (length)"
 *
 * Returns the number of valid fields that were found.
 * 
 * We could verify "name" more but it's probably not worth it - in fact, being
 * more flexible here allows us to use wildcards, because we can verify it 
 * with wildcard matching later.
 */
int parseMatchLine(const char *inputText,
		   char **matchNameOut,
		   int *matchStartOut, 
		   int *matchEndOut, 
		   int *matchLenOut)
{
  char sequence_name[1000] = {'\0'};
  int start = 0, end = 0, length = 0;
  char *format_str = "\"%[^\"]\"%d%d (%d)";

  int fields = sscanf(inputText, format_str, &sequence_name[0], &start, &end, &length);
  
  if (fields < 1)
    {
      /* Try again but without the quotes */
      format_str = "%[^\"]%d%d (%d)";
      fields = sscanf(inputText, format_str, &sequence_name[0], &start, &end, &length);
    }
  
  gboolean foundError = (fields < 1);
  int validFields = 0;

  if (fields > 0)
    {
      *matchNameOut = g_strdup(sequence_name);
      foundError = (sequence_name == NULL);
      
      if (!foundError)
	++validFields;
    }
    
  if (fields > 1)
    {
      *matchStartOut = start;
      foundError = (start < 1);
      
      if (!foundError)
	++validFields;
    }

  if (fields > 2)
    {
      *matchEndOut = end;
      foundError = (start >= end);
      
      if (!foundError)
	++validFields;
    }

  if (fields > 3)
    {
      *matchLenOut = length;
      foundError = (length < 1);
      
      if (!foundError)
	++validFields;
    }

  return fields;
}


/* Parse a list of matches using parseMatchLine(), and extract the names. Returns
 * a GList of sequence names. Both the list and all its entries should be free'd 
 * by the caller. */
GList* parseMatchList(const char *inputText)
{
  GList *matchList = NULL ;

  if (inputText)
    {
      /* Split the input text into separate strings based on the following delimiters: */
      char *delimiters = "\n\r,;";
      char **tokens = g_strsplit_set(inputText, delimiters, -1);   /* -1 means do all tokens. */
      char **token = tokens;
      
      if (token)
	{
	  char *match = *token;

	  while (token && match)
	    {
	      char *match_name ;
	      int start = 0, end = 0, length = 0 ;

	      if (parseMatchLine(match, &match_name, &start, &end, &length) > 0)
		{
		  matchList = g_list_append(matchList, match_name);
		}

	      token++ ;
	      match = *token ? *token : 0; /* token may be empty string if two delimiters next to each other */
	    }
	}
      
      g_strfreev(tokens);
    }
  
  return matchList ;
}


/* Request the text in the PRIMARY clipboard and pass it to the given callback function */
void requestPrimaryClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data)
{
#ifndef __CYGWIN__  /* Not applicable for Windows */

  GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_PRIMARY);
  gtk_clipboard_request_text(clipboard, callback, data);
  
#endif
}


/* Set the text in the PRIMARY clipboard */
void setPrimaryClipboardText(const char *text)
{
#ifndef __CYGWIN__   /* Not applicable for Windows */

  GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_PRIMARY);
  gtk_clipboard_set_text(clipboard, text, -1);
  
#endif
}


/* Request the text in the DEFAULT clipboard and pass it to the given callback function */
void requestDefaultClipboardText(GtkClipboardTextReceivedFunc callback, gpointer data)
{
  GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
  gtk_clipboard_request_text(clipboard, callback, data);
}


/* Set the text in the DEFAULT clipboard */
void setDefaultClipboardText(const char *text)
{
  GtkClipboard *clipboard = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);
  gtk_clipboard_set_text(clipboard, text, -1);
}


/* Returns the pointer to the BlxSequence that matches the given sequence name and
 * strand. Returns NULL if no match was found. */
BlxSequence *findBlxSequence(GList *seqList, const char *reqdName, const BlxStrand reqdStrand)
{
  BlxSequence *result = NULL;
  
  /* Loop through all sequences in the list */
  GList *listItem = seqList;
  
  for ( ; listItem; listItem = listItem->next)
    {
      BlxSequence *currentSeq = (BlxSequence*)(listItem->data);

      if (!strcmp(currentSeq->fullName, reqdName) && currentSeq->strand == reqdStrand)
	{
	  result = currentSeq;
	  break;
	}
    }
  
  return result;
}


/* Create and cache a blank pixmap drawable in the given widget */
GdkDrawable* createBlankPixmap(GtkWidget *widget)
{
  GdkWindow *window = GTK_IS_LAYOUT(widget) ? GTK_LAYOUT(widget)->bin_window : widget->window;
  
  GdkDrawable *drawable = gdk_pixmap_new(window, widget->allocation.width, widget->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(widget, drawable); /* deletes the old drawable, if there is one */
  
  /* Paint a blank rectangle for the background, the same color as the widget's background */
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *bgColor = &widget->style->bg[GTK_STATE_NORMAL];
  gdk_gc_set_foreground(gc, bgColor);
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);
  
  return drawable;
}


/* Utility to pop up a simple confirmation dialog box with the given title and text, 
 * with just an "OK" and "Cancel" button. Blocks until the user responds, and returns
 * their response ID. */
gint runConfirmationBox(GtkWidget *blxWindow, char *title, char *messageText)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  GTK_WINDOW(blxWindow), 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_CANCEL,
						  GTK_RESPONSE_REJECT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  /* Put message and icon into an hbox */
  GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 0);
  
  GtkWidget *image = gtk_image_new_from_stock(GTK_STOCK_DIALOG_QUESTION, GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start(GTK_BOX(hbox), image, TRUE, TRUE, 0);
  
  GtkWidget *label = gtk_label_new(messageText);
  gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
  gtk_widget_show_all(hbox);
  
  gint response = gtk_dialog_run(GTK_DIALOG(dialog));
  
  gtk_widget_destroy(dialog);
  
  return response;
}


/* Lookup the BlxColor with the given id in the given hash table and return a pointer to it */
BlxColor* getBlxColor(GList *colorList, const BlxColorId colorId)
{
  BlxColor *result = NULL;
  GList *item = colorList;
  
  for ( ; item; item = item->next)
    {
      BlxColor *blxColor = (BlxColor*)(item->data);

      if (blxColor->id == colorId)
	{
	  result = blxColor;
	  break;
	}
    }
  
  if (result && result->transparent)
    {
      /* return the background color instead */
      result = getBlxColor(colorList, BLXCOL_BACKGROUND);
    }
  else if (!result)
    {
      printf("Warning: color ID %d not found.\n", colorId);
    }
  
  return result;
}


/* Default handler for GLib log messages (e.g. from g_message() etc.) */
void defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data)
{
  printf("%s", (char*)message);
}


/* Default handler for critical GLib log messages (i.e. from g_error and g_critical.) */
void popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data)
{
  /* Display message to screen and in pop-up box */
  printf("***%s", (char*)message);
  messout((char*)message);
}


/* Tag the given prefix string onto the start of the given error's message */
void prefixError(GError *error, char *formatStr, ...)
{
//  g_prefix_error(error, prefixStr);  /* Only from GLib v2.16 */

  va_list argp;
  va_start(argp, formatStr);
  
  /* Print the prefix text and args into a new string. We're not sure how much space we need
   * for the args, so give a generous buffer but use vsnprintf to make sure we don't overrun.
   * (g_string_vprintf would avoid this problem but is only available from GLib version 2.14). */
  const int len = strlen(formatStr) + 200;
  char tmpStr[len];
  vsnprintf(tmpStr, len, formatStr, argp);

  va_end(argp);

  /* Append the error message */
  char *resultStr = g_malloc(len + strlen(error->message));
  snprintf(resultStr, len, "%s%s", tmpStr, error->message);
  
  /* Replace the error message text with the new string. */
  g_free(error->message);
  error->message = resultStr;
}


//
///* Tag the given postfix string onto the end of the given error's message */
//static void postfixError(GError *error, char *formatStr, ...)
//{
//  va_list argp;
//  va_start(argp, formatStr);
//  
//  /* Print the prefix text into a new string */
//  GString *tmp = g_string_new("");
//  g_string_vprintf(tmp, formatStr, argp);
//  
//  va_end(argp);
//  
//  /* Prepend the error text */
//  g_string_prepend(tmp, error->message);
//  
//  /* Replace the error message text with the new string. */
//  g_free(error->message);
//  error->message = tmp->str;
//  
//  /* Free the temporary GString object */
//  g_string_free(tmp, FALSE);
//}
//
//
///* Custom function to propagate a GError. Unlike g_propagate_error, this allows the destination
// * error to be non-null: if that is the case, the new error message is concatenated on the end
// * of the original error message. The src error is free'd and set to NULL. */
//void propagateError(GError **dest, GError **src)
//{
//  if (dest == NULL || src == NULL || *src == NULL)
//    return;
//  
//  if (*dest == NULL)
//    {
//      g_propagate_error(dest, *src);
//      *src = NULL;
//    }
//  else
//    {
//      postfixError(*error, (*src)->message);
//    }
//}


/* This utility reports the given error if it is non-null and then frees the GError and
 * sets the pointer to NULL. It reports it at the given log level.  */
void reportAndClearIfError(GError **error, GLogLevelFlags log_level)
{
  if (error && *error)
    {
      switch (log_level)
        {
          case G_LOG_LEVEL_ERROR:            g_error("%s", (*error)->message);            break;
          case G_LOG_LEVEL_CRITICAL:         g_critical("%s", (*error)->message);         break;
          case G_LOG_LEVEL_WARNING:          g_warning("%s", (*error)->message);          break;
          case G_LOG_LEVEL_MESSAGE:          g_message("%s", (*error)->message);          break;
          case G_LOG_LEVEL_INFO:             g_message("%s", (*error)->message);          break;
          case G_LOG_LEVEL_DEBUG:            g_debug("%s", (*error)->message);            break;
            
          default:
            break;
        };
      
      g_error_free(*error);
      *error = NULL;
    }
}


#ifdef DEBUG
/* Prints whitespace for log level (that is, how indented logging currently is, according to how
 * many DEBUG_ENTER and EXITs we've called), and then increase it by the increase amount (can be negative). */
void debugLogLevel(const int increaseAmt)
{
  static int count = 0;
  
  gboolean done = FALSE;
  
  /* If decreasing the count (i.e. printing an exit statement) decrease the log level
   * first. */
  if (increaseAmt < 0)
    {
      count += increaseAmt;
      if (count < 0) count = 0;
      done = TRUE;
    }
  
  /* Print whitespace */
  int i = 0;
  for ( ; i < count; ++i)
    {
      printf("  "); /* print 2 spaces per log level */
    }
  
  if (!done)
    {
      count += increaseAmt;
      if (count < 0) count = 0;
    }
}
#endif
