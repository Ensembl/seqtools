/*
 *  utilities.c
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */

#include "SeqTools/utilities.h"
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>


/* Test char to see if it's a iupac dna/peptide code. */
#define ISIUPACDNA(BASE) \
(((BASE) == 'a' || (BASE) == 'c'|| (BASE) == 'g' || (BASE) == 't'       \
  || (BASE) == 'u' || (BASE) == 'r'|| (BASE) == 'y' || (BASE) == 'm'    \
  || (BASE) == 'k' || (BASE) == 'w'|| (BASE) == 's' || (BASE) == 'b'    \
  || (BASE) == 'd' || (BASE) == 'h'|| (BASE) == 'v'                     \
  || (BASE) == 'n' || (BASE) == SEQUENCE_CHAR_GAP || (BASE) == SEQUENCE_CHAR_PAD))

#define ISIUPACPEPTIDE(PEPTIDE) \
(((PEPTIDE) == 'A' || (PEPTIDE) == 'B'|| (PEPTIDE) == 'C' || (PEPTIDE) == 'D'       \
  || (PEPTIDE) == 'E' || (PEPTIDE) == 'F'|| (PEPTIDE) == 'G' || (PEPTIDE) == 'H'    \
  || (PEPTIDE) == 'I' || (PEPTIDE) == 'K'|| (PEPTIDE) == 'L' || (PEPTIDE) == 'M'    \
  || (PEPTIDE) == 'N' || (PEPTIDE) == 'P'|| (PEPTIDE) == 'Q' || (PEPTIDE) == 'R'    \
  || (PEPTIDE) == 'S' || (PEPTIDE) == 'T'|| (PEPTIDE) == 'U' || (PEPTIDE) == 'V'    \
  || (PEPTIDE) == 'W' || (PEPTIDE) == 'X'|| (PEPTIDE) == 'Y' || (PEPTIDE) == 'Z'    \
  || (PEPTIDE) == SEQUENCE_CHAR_STOP || (PEPTIDE) == SEQUENCE_CHAR_GAP || (PEPTIDE) == SEQUENCE_CHAR_PAD))




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
void widgetClearCachedDrawable(GtkWidget *widget, gpointer data)
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
  
  gtk_misc_set_padding(GTK_MISC(label), DEFAULT_LABEL_X_PAD, 0);
  
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
  return (range != NULL && value >= range->min && value <= range->max);
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


/* Creates a GdkColor from the given color string in hex format, e.g. "#00ffbe". Returns false
 * and sets 'error' if there was a problem */
gboolean getColorFromString(const char *colorStr, GdkColor *color, GError **error)
{
  gboolean ok = gdk_color_parse(colorStr, color);

  if (ok)
    {
      gboolean failures[1];
      gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), color, 1, TRUE, TRUE, failures);
      
      if (numFailures > 0)
        {
          ok = FALSE;
        }
    }

  if (!ok)
    {
      g_free(color);
      color = NULL;
      g_set_error(error, BLX_ERROR, 1, "Error parsing color string '%s'\n", colorStr);
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


/* Utility to take a gdk color and return a greyscale version of it */
void convertToGrayscale(GdkColor *origColor, GdkColor *result)
{
  result->red = result->green = result->blue =
   ((11 * origColor->red) + (16 * origColor->green) + (5 * origColor->blue)) / 32;
   
  gboolean failures[1];
  gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), result, 1, TRUE, TRUE, failures);
  
  if (numFailures > 0)
    {
      printf("Warning: error calculating greyscale color (RGB=%d %d %d).", origColor->red, origColor->green, origColor->blue);
      
      /* Set it to black, for want of something better to do... */
      result->pixel = 0;
      result->red = 0;
      result->green = 0;
      result->blue = 0;
    }

}


/* Utility to take a gdk color and return a much darker shade of it, for use as a drop-shadow */
void getDropShadowColor(GdkColor *origColor, GdkColor *result)
{
  const double factor = 0.3;
  adjustColorBrightness(origColor, factor, result);
}


gboolean typeIsExon(const BlxMspType mspType)
{
  return (mspType == BLXMSP_CDS || mspType == BLXMSP_UTR || mspType == BLXMSP_EXON);
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
					 const gboolean fill)
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
      result = mspGetColor(prevExon, defaultColors, blxSeq, selected, usePrintColors, fill);
    }
  else if (nextIsUtr)
    {
      result = mspGetColor(nextExon, defaultColors, blxSeq, selected, usePrintColors, fill);
    }
  else if (prevExon)
    {
      /* Both exons (or the sole exon, if only one exists) are CDS, so the intron is CDS */
      result = mspGetColor(prevExon, defaultColors, blxSeq, selected, usePrintColors, fill);
    }
  else if (nextExon)
    {
      /* This is the only exon and it is CDS, so make the intron CDS */
      result = mspGetColor(nextExon, defaultColors, blxSeq, selected, usePrintColors, fill);
    }
  else
    {
      /* No exon exists adjacent to this intron: default to generic exon color for want of anything better to do. */
      if (fill)
        result = getGdkColor(BLXCOLOR_EXON_FILL, defaultColors, selected, usePrintColors);
      else
        result = getGdkColor(BLXCOLOR_EXON_LINE, defaultColors, selected, usePrintColors);
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

gboolean mspIsPolyATail(const MSP const *msp)
{
  return (msp && msp->type == BLXMSP_POLYA_TAIL);
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
  return mspIsExon(msp) || mspIsBlastMatch(msp) || mspIsPolyATail(msp);
}

/* Whether the MSP requires subject sequence strand to be set. Only matches 
 * require strand on the subject sequence, although exons may have them set. */
gboolean mspHasSStrand(const MSP const *msp)
{
  return mspIsBlastMatch(msp) || mspIsPolyATail(msp);
}

/* Whether the MSP requires the actual sequence data for the subject sequence. Only
 * matches require the sequence data. */
gboolean mspHasSSeq(const MSP  const *msp)
{
  return mspIsBlastMatch(msp) || mspIsPolyATail(msp);
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
	     const gboolean displayRev,
	     const gboolean showUnalignedSeq,
	     const gboolean limitUnalignedBases,
	     const int numUnalignedBases)
{
  int result = UNSET_INT;
  
  if (mspIsBlastMatch(msp) || mspIsExon(msp) || mspIsPolyATail(msp))
    {
      const gboolean qForward = (mspGetRefStrand(msp) == BLXSTRAND_FORWARD);
      const gboolean sForward = (mspGetMatchStrand(msp) == BLXSTRAND_FORWARD);
      const gboolean sameDirection = (qForward == sForward);
      const gboolean inGapsRange = (qIdx >= msp->qRange.min && qIdx <= msp->qRange.max);
      
      if (msp->gaps && g_slist_length(msp->gaps) >= 1 && inGapsRange)
	{
	  /* Gapped alignment. Look to see if x lies inside one of the "gaps" ranges. */
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
      else if (!inGapsRange && mspIsBlastMatch(msp))
	{
	  /* The q index is outside the alignment range but the option to show unaligned sequence 
	  * is enabled, so extend the range to include the allowed number of extra bases (or the 
	  * full range of the match sequence, if no limit is set), and check if the q index is 
	  * within that. */
	  int sMin = msp->sRange.min;
	  int sMax = msp->sRange.max;

	  if (!inGapsRange && showUnalignedSeq && mspGetMatchSeq(msp))
	    {
	      /* Get the full range of the match sequence */
	      sMin = 1;
	      sMax = mspGetMatchSeqLen(msp);
	      
	      if (limitUnalignedBases)
		{
		  /* Include up to 'numUnalignedBases' each side of the MSP range. (Make sure
		   * this doesn't take us outside the full range of the match sequence, though.) */
		  sMin = max(sMin, msp->sRange.min - numUnalignedBases);
		  sMax = min(sMax, msp->sRange.max + numUnalignedBases);
		}
	    }
	    
	  if (qIdx < msp->qRange.min)
	    {
	      /* Find the offset backwards from the low coord and subtract it from the low end
	       * of the s coord range (or add it to the high end, if the directions are opposite). */
	      const int offset = (msp->qRange.min - qIdx) / numFrames;
	      result = sameDirection ? msp->sRange.min - offset : msp->sRange.max + offset;
	    }
	  else
	    {
	      /* Find the offset forwards from the high coord and add it to the high end of the
	       * s coord range (or subtract it from the low end, if the directions are opposite). */
	      const int offset = (qIdx - msp->qRange.max) / numFrames;
	      result = sameDirection ? msp->sRange.max + offset : msp->sRange.min - offset;
	    }

	  /* If the result is still out of range, there's nothing to show for this qIdx. */
	  if (result < sMin || result > sMax)
	    {
	      result = UNSET_INT;
	    }
	}
      else
	{
	  /* If strands are in the same direction, find the offset from qRange.min and add it to 
	   * sRange.min. If strands are in opposite directions, find the offset from qRange.min and
	   * subtract it from sRange.max. Note that the offset could be negative if we're outside
	   * the alignment range. */
	  int offset = (qIdx - msp->qRange.min)/numFrames ;
	  result = (sameDirection) ? msp->sRange.min + offset : msp->sRange.max - offset ;
	  
	  if (result < msp->sRange.min || result > msp->sRange.max)
	    {
	      result = UNSET_INT;
	    }

	}
    }
  
  return result;
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


static void appendTextIfNonNull(GString *gstr, const char *separator, const char *text)
{
  if (text)
    {
      g_string_append_printf(gstr, "%s%s", separator, text);
    }
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


/* For the given BlxColor, return a pointer to the GdkColor that meets the given criteria */
static const GdkColor *blxColorGetColor(const BlxColor *blxColor, const gboolean selected, const gboolean usePrintColors)
{
  const GdkColor *result = NULL;
  
  if (usePrintColors)
    {
      if (selected)
	{
	  result = &blxColor->printSelected;
	}
      else
	{
	  result = &blxColor->print;
	}
    }
  else
    {
      if (selected)
	{
	  result = &blxColor->selected;
	}
      else
	{
	  result = &blxColor->normal;
	}
    }
    
  return result;
}


/* Return the color matching the given properties from the given style */
static const GdkColor *styleGetColor(const BlxStyle *style, const gboolean selected, const gboolean usePrintColors, const gboolean fill)
{
  const GdkColor *result = NULL;
  
  if (fill)
    {
      result = blxColorGetColor(&style->fillColor, selected, usePrintColors);
    }
  else
    {
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
			    const gboolean fill)
{
  const GdkColor *result = NULL;
  
  if (msp->style)
    {
      result = styleGetColor(msp->style, selected, usePrintColors, fill);
    }

  if (!result)
    {
      /* Use the default color for this MSP's type */
      switch (msp->type)
        {
          case BLXMSP_EXON:
            result = getGdkColor(fill ? BLXCOLOR_EXON_FILL : BLXCOLOR_EXON_LINE, defaultColors, selected, usePrintColors);
            break;
            
            case BLXMSP_CDS:
            result = getGdkColor(fill ? BLXCOLOR_CDS_FILL : BLXCOLOR_CDS_LINE, defaultColors, selected, usePrintColors);
            break;

          case BLXMSP_UTR:
            result = getGdkColor(fill ? BLXCOLOR_UTR_FILL : BLXCOLOR_UTR_LINE, defaultColors, selected, usePrintColors);
            break;
            
          case BLXMSP_INTRON:
             result = mspGetIntronColor(msp, defaultColors, blxSeq, selected, usePrintColors, fill);
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


/* Return summary info about a given MSP (e.g. for displaying in the status bar). The
 * result should be free'd with g_free. */
char* mspGetSummaryInfo(const MSP const *msp)
{
  char *result = NULL;
  
  if (msp)
    {
      GString *resultStr = g_string_new("");
      const char *separator = ";  ";
      
      g_string_append_printf(resultStr, "%s", mspGetSName(msp));
      appendTextIfNonNull(resultStr, separator, mspGetOrganism(msp));
      appendTextIfNonNull(resultStr, separator, mspGetGeneName(msp));
      appendTextIfNonNull(resultStr, separator, mspGetTissueType(msp));
      appendTextIfNonNull(resultStr, separator, mspGetStrain(msp));
      
      result = resultStr->str;
      g_string_free(resultStr, FALSE);
    }
  
  return result;
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
  char *result = NULL;

  const int len = strlen(longName);

  if (len > 0)
    {
      result = g_malloc(sizeof(char) * len);
      
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
    }
  
  return result;
}


/* Return the full name of a BlxSequence (including prefix and variant) */
const char *blxSequenceGetFullName(const BlxSequence *seq)
{
  return seq->fullName;
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
  if (!fullName && !idTag)
    {
      g_set_error(error, BLX_ERROR, 1, "Cannot create sequence: ID or name must be specified.");
      return NULL;
    }
  
  BlxSequence *seq = g_malloc(sizeof(BlxSequence));
  
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


/* Converts the given double to a string. The result must be free'd with g_free */
char* convertDoubleToString(const gdouble value)
{
  char result[numDigitsInInt((int)value) + 3]; /* the +3 includes decimal point, one decimal place, and terminating null */
  sprintf(result, "%1.1f", value);
  return g_strdup(result);
}


/* Converts the given string to an integer. */
int convertStringToInt(const char *inputStr)
{
  return atoi(inputStr);
}


/* Utility to return true if the given char is a whitespace char */
gboolean isWhitespaceChar(const char curChar)
{
  return (curChar == ' ' || curChar == '\t');
}


/***********************************************************
 *			   Dialogs			   *
 ***********************************************************/

/* Utility to set the height of a GtkTextView based on the number of lines of text it
 * contains (but not going above the given max height). Returns the height that was set. */
//static int textViewSetHeight(GtkWidget *textView, GtkTextBuffer *textBuffer, PangoFontDescription *fontDesc, int maxHeight, int *charHeightOut)
//{
//  PangoContext *context = gtk_widget_get_pango_context(textView);
//  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
//  gint charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
//  
//  if (charHeightOut)
//    {
//      *charHeightOut = charHeight;
//    }
//  
//  /* Adjust height to include all the lines if possible (limit to the original height passed in thought) */
//  const int calcHeight = (gtk_text_buffer_get_line_count(textBuffer) * charHeight);
//  int resultHeight = min(maxHeight, calcHeight);
//  
//  pango_font_metrics_unref(metrics);
//  
//  return resultHeight;
//}


/* Utility to create a scrollable text view from the given message text. If textBufferOut is
 * non-null its value is set to point to the text buffer. */
GtkWidget* createScrollableTextView(const char *messageText,
				    const gboolean wrapText,
				    PangoFontDescription *fontDesc,
                                    const gboolean useMarkup,
				    int *height,
                                    GtkTextView **textViewOut)
{
  /* Create a text buffer and copy the text in */
  GtkTextBuffer *textBuffer = gtk_text_buffer_new(NULL);
  
  GtkTextIter textIter;
  gtk_text_buffer_get_iter_at_line(textBuffer, &textIter, 0);
  
  if (useMarkup && messageText)
    {
      gtk_text_buffer_insert_markup(textBuffer, &textIter, messageText);
    }
  else if (messageText)
    {
      gtk_text_buffer_set_text(textBuffer, messageText, -1);
    }
  
  /* Create a text view to display the buffer */
  GtkWidget *textView = gtk_text_view_new();
  gtk_text_view_set_buffer(GTK_TEXT_VIEW(textView), textBuffer);
  gtk_text_view_set_editable(GTK_TEXT_VIEW(textView), FALSE);
  
  if (textViewOut)
    {
      *textViewOut = GTK_TEXT_VIEW(textView);
    }
  
  if (fontDesc)
    {
      gtk_widget_modify_font(textView, fontDesc);
    }
  
  if (wrapText)
    {
      gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(textView), TRUE);
    }
  
  /* Put the text view in a scrollable window */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(scrollWin), textView);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  /* Set the height to fit the number of lines of text, unless this would exceed the passed-in height */
//  int charHeight = UNSET_INT;
//  *height = textViewSetHeight(textView, textBuffer, fontDesc, *height, &charHeight);
  
  /* Return the outermost container */
  return scrollWin;
}


/* Utility to pop up a simple dialog with the given title and text, with just an "OK" 
 * button. Returns the dialog, and optionally sets a pointer to the text buffer in the 
 * textBuffer return arg. */
GtkWidget* showMessageDialog(const char *title,  
                             const char *messageText,
                             GtkWidget *parent,
                             const int width,
                             const int maxHeight,
                             const gboolean wrapText,
                             const gboolean useMarkup,
                             PangoFontDescription *fontDesc,
                             GtkTextView **textView)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  GTK_WINDOW(parent), 
						  GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);

  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);

  int height = maxHeight;
  GtkWidget *child = createScrollableTextView(messageText, wrapText, fontDesc, useMarkup, &height, textView);

  gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), child, TRUE, TRUE, 0);

  g_signal_connect(dialog, "response", G_CALLBACK(gtk_widget_destroy), NULL);
  gtk_widget_show_all(dialog);
  
  return dialog;
}


/* Get/set a CallbackData struct for this widget */
static CallbackData* widgetGetCallbackData(GtkWidget *widget)
{
  return widget ? (CallbackData*)(g_object_get_data(G_OBJECT(widget), "callbackData")) : NULL;
}


/* Set callback function and user-data for the given widget */
void widgetSetCallbackData(GtkWidget *widget, BlxResponseCallback func, gpointer data)
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
static void widgetCallCallback(GtkWidget *widget, const gint responseId)
{
  CallbackData *callbackData = widgetGetCallbackData(widget);
  
  if (callbackData && callbackData->func)
    {
      callbackData->func(widget, responseId, callbackData->data);
    }
}

/* Call the stored CallbackData callbacks (if any exist) for this widget and
 * all of its children. The response ID is passed as the user data */
void widgetCallAllCallbacks(GtkWidget *widget, gpointer data)
{
  if (widget && GTK_IS_WIDGET(widget))
    {
      gint responseId = GPOINTER_TO_INT(data);
      
      widgetCallCallback(widget, responseId);
  
      if (GTK_IS_CONTAINER(widget))
	{
	  gtk_container_foreach(GTK_CONTAINER(widget), widgetCallAllCallbacks, data);
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
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      destroy = TRUE;
      break;
      
    case GTK_RESPONSE_APPLY:
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
      destroy = FALSE;
      break;

    case BLX_RESPONSE_BACK:  /* fall through */
    case BLX_RESPONSE_FORWARD:
      widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
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
      /* If it's a persistent dialog, just hide it, otherwise destroy it */
      const gboolean isPersistent = GPOINTER_TO_INT(data);
      
      if (isPersistent)
        {
          gtk_widget_hide_all(GTK_WIDGET(dialog));
        }
      else
        {
          gtk_widget_destroy(GTK_WIDGET(dialog));
        }
    }
}


/* Utility to remove all of the child widgets from a dialog's content area. */
void dialogClearContentArea(GtkDialog *dialog)
{
  GtkWidget *contentArea = dialog->vbox;
  GtkWidget *actionArea = dialog->action_area;
  
  GList *child = gtk_container_get_children(GTK_CONTAINER(contentArea));
  
  for ( ; child; child = child->next)
    {
      GtkWidget *childWidget = GTK_WIDGET(child->data);
      
      if (childWidget != actionArea)
        {
          gtk_widget_destroy(childWidget);
          child->data = NULL;
        }
    }
  
//  gtk_container_foreach(GTK_CONTAINER(contentArea), destroyWidget, NULL);
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
  if (!textToSearch || !searchStr)
    {
      return 0;
    }
  
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
	    
	    while (toupper(*textChar) != toupper(*searchChar))
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
	    if (toupper(*searchChar++) != toupper(*textChar++))
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
  
  if (!str1 && !str2)
    {
      result = TRUE;
    }
  else if (!str1 || !str2)
    {
      result = FALSE;
    }
  else if (caseSensitive)
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


/* If error is non-null, tag the given prefix string onto the start of the given 
 * error's message. If error is null, do nothing. */
void prefixError(GError *error, char *formatStr, ...)
{
//  g_prefix_error(error, prefixStr);  /* Only from GLib v2.16 */

  if (error)
    {
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
}


/* Tag the given postfix string onto the end of the given error's message */
void postfixError(GError *error, char *formatStr, ...)
{
  va_list argp;
  va_start(argp, formatStr);
  
  
  /* Print the postfix text and args into a new string. We're not sure how much space we need
   * for the args, so give a generous buffer but use vsnprintf to make sure we don't overrun.
   * (g_string_vprintf would avoid this problem but is only available from GLib version 2.14). */
  const int len = strlen(formatStr) + 200;
  char tmpStr[len];
  vsnprintf(tmpStr, len, formatStr, argp);
  
  va_end(argp);
  
  /* Prepend the error message */
  char *resultStr = g_malloc(len + strlen(error->message));
  snprintf(resultStr, len, "%s%s", error->message, tmpStr);
  
  /* Replace the error message text with the new string. */
  g_free(error->message);
  error->message = resultStr;
}


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


/* Set the values in an IntRange. */
void intrangeSetValues(IntRange *range, const int val1, const int val2)
{
  range->min = val1 < val2 ? val1 : val2;
  range->max = val1 < val2 ? val2 : val1;
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


/* Draw the highlight box (for highlighting the current detail-view display range on 
 * big-picture widgets such as the grids and exon views). */
void drawHighlightBox(GdkDrawable *drawable,
                      const GdkRectangle const *rect,
                      const gint minWidth, 
                      GdkColor *color)
{
  GdkGC *gc = gdk_gc_new(drawable);

  gdk_gc_set_foreground(gc, color);
//  gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_SOLID, GDK_CAP_BUTT, GDK_JOIN_MITER);

  gdk_gc_set_function(gc, GDK_AND);
  
  const int width = max(minWidth, rect->width);
  
  gdk_draw_rectangle(drawable, gc, TRUE, 
		     rect->x, rect->y, width, rect->height);
                     
  g_object_unref(gc);
}


/* Create a given string from the given format and args. Like printf but it creates a string
 * of the correct length for you. The result should be free'd with g_free. */
char* blxprintf(char *formatStr, ...)
{
  va_list argp;
  va_start(argp, formatStr);

  /* Print the format text and args into a new string. We're not sure how much space we need
   * for the args, so give a generous buffer but use vsnprintf to make sure we don't overrun.
   * (g_string_vprintf would avoid this problem but is only available from GLib version 2.14). */
  const int len = strlen(formatStr) + 200;
  char *resultStr = g_malloc(len);
  vsnprintf(resultStr, len, formatStr, argp);
  
  va_end(argp);

  return resultStr;
}



/* Returns true if the given char is a valid IUPAC character of the given type */
gboolean isValidIupacChar(const char inputChar, const BlxSeqType seqType)
{
  gboolean isValidChar = FALSE;

  if (seqType == BLXSEQ_DNA)
    {
      /* Convert to correct case and then check if this letter is a valid base */
      char cp = tolower(inputChar);
      isValidChar = ISIUPACDNA(cp);
    }
  else if (seqType == BLXSEQ_PEPTIDE)
    {
      /* Convert to correct case and then check if this letter is a valid peptide */
      char cp = toupper(inputChar);
      isValidChar = ISIUPACPEPTIDE(cp);
    }

  return isValidChar;
}


/* Set the shadow style for the given status bar. Give the shadow style as described
 * in the GTK documentation but as a string, e.g. "GTK_SHADOW_IN". */
void setStatusBarShadowStyle(GtkWidget *statusBar, const char *shadowStyle)
{
  const char name[] = "StatusbarName";
  gtk_widget_set_name(statusBar, name);

  char parseString[500];
  sprintf(parseString, "style \"detailViewStatusbar\"\n"
          "{\n"
          "GtkStatusbar::shadow-type	      = %s\n"
          "}"
          "widget \"*%s*\" style \"detailViewStatusbar\"", shadowStyle, name);
  gtk_rc_parse_string(parseString);
}


/***********************************************************
 * Customisation to GtkTextView to allow pango markup.     
 * 
 * Taken from https://bugzilla.gnome.org/show_bug.cgi?id=59390
 * since this code hasn't been incorporated into GTK yet...
 * 
 ***********************************************************/

static void gtk_text_buffer_real_insert_markup         (GtkTextBuffer     *buffer,
                                                        GtkTextIter       *textiter,
                                                        const gchar       *markup,
                                                        GtkTextTag        *extratag);

static void
gtk_text_buffer_real_insert_markup (GtkTextBuffer *buffer,
                                    GtkTextIter   *textiter,
                                    const gchar   *markup,
                                    GtkTextTag    *extratag)
{
 PangoAttrIterator  *paiter;
  PangoAttrList      *attrlist;
  GtkTextMark        *mark;
  GError             *error = NULL;
  gchar              *text;

  g_return_if_fail (GTK_IS_TEXT_BUFFER (buffer));
  g_return_if_fail (textiter != NULL);
  g_return_if_fail (markup != NULL);
  g_return_if_fail (gtk_text_iter_get_buffer (textiter) == buffer);

  if (*markup == '\000')
    return;

  if (!pango_parse_markup(markup, -1, 0, &attrlist, &text, NULL, &error))
    {
      g_warning("Invalid markup string: %s", error->message);
      g_error_free(error);
      return;
    }

  if (attrlist == NULL)
    {
      gtk_text_buffer_insert(buffer, textiter, text, -1);
      g_free(text);
      return;
    }

  /* create mark with right gravity */
  mark = gtk_text_buffer_create_mark(buffer, NULL, textiter, FALSE);

  paiter = pango_attr_list_get_iterator(attrlist);

  do
    {
      PangoAttribute *attr;
      GtkTextTag     *tag;
      gint            start, end;

      pango_attr_iterator_range(paiter, &start, &end);

      if (end == G_MAXINT)  /* last chunk */
        end = start-1; /* resulting in -1 to be passed to _insert */

      tag = gtk_text_tag_new(NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_LANGUAGE)))
        g_object_set(tag, "language", pango_language_to_string(((PangoAttrLanguage*)attr)->value), NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_FAMILY)))
        g_object_set(tag, "family", ((PangoAttrString*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_STYLE)))
        g_object_set(tag, "style", ((PangoAttrInt*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_WEIGHT)))
        g_object_set(tag, "weight", ((PangoAttrInt*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_VARIANT)))
        g_object_set(tag, "variant", ((PangoAttrInt*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_STRETCH)))
        g_object_set(tag, "stretch", ((PangoAttrInt*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_SIZE)))
        g_object_set(tag, "size", ((PangoAttrInt*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_FONT_DESC)))
        g_object_set(tag, "font-desc", ((PangoAttrFontDesc*)attr)->desc, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_FOREGROUND)))
        {
          GdkColor col = { 0,
                           ((PangoAttrColor*)attr)->color.red,
                           ((PangoAttrColor*)attr)->color.green,
                           ((PangoAttrColor*)attr)->color.blue
                         };

          g_object_set(tag, "foreground-gdk", &col, NULL);
        }

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_BACKGROUND)))
        {
          GdkColor col = { 0,
                           ((PangoAttrColor*)attr)->color.red,
                           ((PangoAttrColor*)attr)->color.green,
                           ((PangoAttrColor*)attr)->color.blue
                         };

          g_object_set(tag, "background-gdk", &col, NULL);
        }

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_UNDERLINE)))
        g_object_set(tag, "underline", ((PangoAttrInt*)attr)->value, NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_STRIKETHROUGH)))
        g_object_set(tag, "strikethrough", (gboolean)(((PangoAttrInt*)attr)->value != 0), NULL);

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_RISE)))
        g_object_set(tag, "rise", ((PangoAttrInt*)attr)->value, NULL);

      /* PANGO_ATTR_SHAPE cannot be defined via markup text */

      if ((attr = pango_attr_iterator_get(paiter, PANGO_ATTR_SCALE)))
        g_object_set(tag, "scale", ((PangoAttrFloat*)attr)->value, NULL);

      gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), tag);

      if (extratag)
        {
          gtk_text_buffer_insert_with_tags(buffer, textiter, text+start, end - start, tag, extratag, NULL);
        }
      else
        {
          gtk_text_buffer_insert_with_tags(buffer, textiter, text+start, end - start, tag, NULL);
        }

      /* mark had right gravity, so it should be
       *  at the end of the inserted text now */
      gtk_text_buffer_get_iter_at_mark(buffer, textiter, mark);
    }
  while (pango_attr_iterator_next(paiter));

  gtk_text_buffer_delete_mark(buffer, mark);
  pango_attr_iterator_destroy(paiter);
  pango_attr_list_unref(attrlist);
  g_free(text);
}


/**
 * gtk_text_buffer_insert_markup:
 * @buffer: a #GtkTextBuffer
 * @iter: a position in the buffer
 * @markup: nul-terminated UTF-8 format text with pango markup to insert
 *
 * Inserts the text in @markup at position @iter. @markup
 * will be inserted in its entirety and must be nul-terminated
 * and valid UTF-8. Emits the "insert_text" signal, possibly multiple
 * times; insertion actually occurs in the default handler for
 * the signal. @iter will point to the end of the inserted text
 * on return
 *
 * Since: 2.6
 **/

void
gtk_text_buffer_insert_markup (GtkTextBuffer *buffer,
                               GtkTextIter   *iter,
                               const gchar   *markup)
{
  gtk_text_buffer_real_insert_markup (buffer, iter, markup, NULL);
}


/**
 * gtk_text_buffer_insert_markup_with_tag:
 * @buffer: a #GtkTextBuffer
 * @iter: a position in the buffer
 * @markup: nul-terminated UTF-8 format text in pango markup format to insert
 * @tag: additional text tag to apply to the whole text
 *
 * Just like <literal>gtk_text_buffer_insert_markup</literal>, only
 * that an additional tag can be specified that is applied to the
 * whole text to be inserted. This is useful to pass formatting
 * options to the text buffer that cannot be specified with
 * pango markup (e.g. text justification or wrap mode).
 * @markup will be inserted in its entirety and must be
 * nul-terminated and valid UTF-8 format
 *
 * Since: 2.6
 **/

void
gtk_text_buffer_insert_markup_with_tag (GtkTextBuffer *buffer,
                                        GtkTextIter   *iter,
                                        const gchar   *markup,
                                        GtkTextTag    *tag)
{
  gtk_text_buffer_real_insert_markup (buffer, iter, markup, tag);
}


/**
 * gtk_text_buffer_set_markup:
 * @buffer: a #GtkTextBuffer
 * @markup: nul-terminated UTF-8 text with pango markup to insert
 *
 * Deletes current contents of @buffer, and inserts the text
 * in @markup instead, which may contain pango markup.
 * @markup will be inserted in its entirety and must be
 * nul-terminated and valid UTF-8.
 *
 * Since: 2.6
 **/

void
gtk_text_buffer_set_markup_with_tag (GtkTextBuffer *buffer,
                                     const gchar   *markup,
                                     GtkTextTag    *tag)
{
  GtkTextIter  start, end;

  g_return_if_fail (GTK_IS_TEXT_BUFFER (buffer));
  g_return_if_fail (markup != NULL);

  if (*markup == '\000')
    return;

  gtk_text_buffer_get_bounds (buffer, &start, &end);

  gtk_text_buffer_delete (buffer, &start, &end);

  gtk_text_buffer_get_iter_at_offset (buffer, &start, 0);
  gtk_text_buffer_insert_markup_with_tag (buffer, &start, markup, tag);
}

/**
 * gtk_text_buffer_set_markup:
 * @buffer: a #GtkTextBuffer
 * @markup: nul-terminated UTF-8 text with pango markup to insert
 *
 * Deletes current contents of @buffer, and inserts the text
 * in @markup instead, which may contain pango markup and must
 * be valid UTF-8. @markup will be inserted in its entirety and
 * must be nul-terminated.
 *
 * Since: 2.6
 **/

void
gtk_text_buffer_set_markup (GtkTextBuffer *buffer,
                            const gchar   *markup)
{
  gtk_text_buffer_set_markup_with_tag(buffer, markup, NULL);
}

