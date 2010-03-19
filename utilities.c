/*
 *  utilities.c
 *  acedb
 *
 *  Created by Gemma Barson on 05/01/2010.
 *
 */

#include "SeqTools/utilities.h"

static CallbackData* widgetGetCallbackData(GtkWidget *widget);


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
      messerror("onExposeLabel: invalid widget type. Expected label [widget=%x]", label);
    }
  
  /* Only widgets that have a pixmap set will be shown in print output */
  GtkWidget *parent = gtk_widget_get_parent(label);
  GdkDrawable *drawable =  gdk_pixmap_new(parent->window, label->allocation.width, label->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(label, drawable);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GtkStyle *style = gtk_widget_get_style(label);
  GdkColor *bgColour = &style->bg[GTK_STATE_NORMAL];
  gdk_gc_set_foreground(gc, bgColour);
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
  return range->min + roundNearest(((double)(range->max - range->min)) / 2.0);
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
  
  return count;
}


/* Create a GdkColor from one of our internally-defined colours, e.g. GDK_YELLOW.
 * At the moment this is just a wrapper around gdk_color_parse. */
GdkColor getGdkColor(const char *colour)
{
  GdkColor result;
  
  gboolean ok = gdk_color_parse(colour, &result);

  if (ok)
    {
      gboolean failures[1];
      gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), &result, 1, TRUE, TRUE, failures);
      
      if (numFailures > 0)
	ok = FALSE;
    }

  if (!ok)
    {
      messerror("Error parsing colour string '%s'", colour);
    }
  
  return result;
}


gboolean mspIsExon(const MSP const *msp)
{
  return (msp && msp->score == -1);
}

gboolean mspIsIntron(const MSP const *msp)
{
  return (msp && msp->score == -2);
}

gboolean mspIsBlastMatch(const MSP const *msp)
{
  return (msp && msp->score > 0);
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
void getSMapMapRangeExtents(SMapMap *range, int *qRangeMin, int *qRangeMax, int *sRangeMin, int *sRangeMax)
{
  if (qRangeMin)
    *qRangeMin = range->r1 < range->r2 ? range->r1 : range->r2;
  
  if (qRangeMax)
    *qRangeMax = range->r1 < range->r2 ? range->r2 : range->r1;
  
  if (sRangeMin)
    *sRangeMin = range->s1 < range->s2 ? range->s1 : range->s2;
  
  if (sRangeMax)
    *sRangeMax = range->s1 < range->s2 ? range->s2 : range->s1;
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
			      const gboolean rightToLeft,
			      const IntRange const *refSeqRange)
{
  int dnaIdx = displayIdx;
  
  if (displaySeqType == BLXSEQ_PEPTIDE)
    {
      /* Convert the input peptide coord to a dna coord */
      dnaIdx = (displayIdx * numFrames) - numFrames + frame + baseNum - 1;
    }
  
  /* If the display is reversed, we need to invert the result. (For example, if the 
   * result is index '2' out of the range '12345', then we convert it to '4', which is the
   * equivalent position in the range '54321'. */
  if (rightToLeft)
    {
      dnaIdx = refSeqRange->max - dnaIdx + refSeqRange->min;
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
			      const gboolean rightToLeft,
			      const IntRange const *dnaIdxRange,
			      int *baseNum)
{
  int displayIdx = dnaIdx;
  
  /* If the display is reversed (i.e. showing increasing values from right-to-left),
   * invert the index (i.e. as if it were the normal left-to-right index for this
   * same position - for example, if the index is '4' out of the range '54321', convert
   * it to '2', which is the equivalent position in the range '12345'. */
  if (rightToLeft)
    {
      displayIdx = dnaIdxRange->max - dnaIdx + dnaIdxRange->min;
    }
  
  if (displaySeqType == BLXSEQ_PEPTIDE)
    {
      /* We're displaying peptides, so convert the dna coord to a peptide coord */
      gdouble fraction = ((gdouble)(displayIdx - frame + 1) / (gdouble)numFrames) ;
      displayIdx = ceil(fraction);
  
      /* Find the base number of this DNA coord within the codon, if requested */
      if (baseNum)
	{
	  *baseNum = roundNearest((fraction - (int)fraction) * numFrames);
      
	  if (*baseNum == 0)
	    {
	      *baseNum = 3;
	    }
	}
    }
  else if (baseNum)
    {
      /* For dna sequences, we only have one reading frame */
      *baseNum = 1;
    }
  
  return displayIdx;
}


/* Get the start coord in the given display range, reversed as necessary if 'reverse' is true.
 * Result is a coord in the DNA sequence - converts as necessary if the display range is in terms
 * of peptide coords */
int getStartDnaCoord(const IntRange const *displayRange, 
		     const int frame,
		     const BlxSeqType displaySeqType, 
		     const gboolean rightToLeft, 
		     const int numFrames,
		     const IntRange const *refSeqRange)
{
  int result = displayRange->min;
  
  /* Convert the display coord to coords into the ref seq, which is a DNA sequence. We want
   * the first base in the codon, if this is a peptide sequence. */
  const int baseNum = 1;
  result = convertDisplayIdxToDnaIdx(result, displaySeqType, frame, baseNum, numFrames, rightToLeft, refSeqRange);
  
  return result;
}


/* Get the end coord in the given display range, reversed as necessary if 'reverse' is true.
 * Result is a coord in the DNA sequence - converts as necessary if the display range is in terms
 * of peptide coords */
int getEndDnaCoord(const IntRange const *displayRange, 
		   const int frame,
		   const BlxSeqType displaySeqType, 
		   const gboolean rightToLeft, 
		   const int numFrames,
		   const IntRange const *refSeqRange)
{
  int result = displayRange->max;
  
  /* Convert the display coord to coords into the ref seq, which is a DNA sequence. We want
   * the last base in the codon, if this is a peptide sequence. */
  const int baseNum = numFrames;
  result = convertDisplayIdxToDnaIdx(result, displaySeqType, frame, baseNum, numFrames, rightToLeft, refSeqRange);
  
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
	     const Strand strand, 
	     const gboolean rightToLeft, 
	     int *nearestIdx)
{
  int result = UNSET_INT;
  
  int qSeqMin = UNSET_INT, qSeqMax = UNSET_INT, sSeqMin = UNSET_INT, sSeqMax = UNSET_INT;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
  
  BOOL qForward = ((strchr(msp->qframe, '+'))) ? TRUE : FALSE ;
  BOOL sForward = ((strchr(msp->sframe, '+'))) ? TRUE : FALSE ;
  BOOL sameDirection = (qForward == sForward);
  
  Array gaps = msp->gaps ;
  if (!gaps || arrayMax(gaps) < 1)
    {
      /* If strands are in the same direction, find the offset from qSeqMin and add it to sSeqMin.
       * If strands are in opposite directions, find the offset from qSeqMin and subtract it from sSeqMax. */
      int offset = (qIdx - qSeqMin)/numFrames ;
      result = (sameDirection) ? sSeqMin + offset : sSeqMax - offset ;
      
      if (result < sSeqMin)
	{
	  result = UNSET_INT;
	  
	  if (nearestIdx)
	    {
	      *nearestIdx = sSeqMin;
	    }
	}
      else if (result > sSeqMax)
	{
	  result = UNSET_INT;
	  
	  if (nearestIdx)
	    {
	      *nearestIdx = sSeqMax;
	    }
	}
      else if (nearestIdx)
	{
	  *nearestIdx = result;
	}
    }
  else
    {
      gboolean gapsReversedWrtDisplay = (rightToLeft == qForward);
      
      /* Look to see if x lies inside one of the reference sequence ranges. */
      int i = 0;
      for ( ; i < arrayMax(gaps) ; ++i)
	{
	  SMapMap *curRange = arrp(gaps, i, SMapMap) ;
	  
	  int qRangeMin, qRangeMax, sRangeMin, sRangeMax;
	  getSMapMapRangeExtents(curRange, &qRangeMin, &qRangeMax, &sRangeMin, &sRangeMax);
	  
	  /* We've "found" the value if it's in or before this range. Note that the
	   * the range values are in decreasing order if the q strand is reversed. */
	  BOOL found = qForward ? qIdx <= qRangeMax : qIdx >= qRangeMin;
	  	  if (found)
	    {
	      BOOL inRange = qForward ? qIdx >= qRangeMin : qIdx <= qRangeMax;
	      
	      if (inRange)
		{
		  /* It's inside this range. Calculate the actual index. */
		  int offset = (qIdx - qRangeMin) / numFrames;
		  result = sameDirection ? sRangeMin + offset : sRangeMax - offset;
		  
		  if (nearestIdx)
		    {
		      *nearestIdx = result;
		    }
		}
	      else if (nearestIdx && (i == 0 || gapsReversedWrtDisplay))
		{
		  /* Remember the start of the current range. This will be the value we
		   * return if the index is before of the first range, or if we're in a gap
		   * and gaps are ordered in the opposite direction to the display (i.e. gap
		   * coords go from low to high and display shows high to low or v.v.) */
		  *nearestIdx = curRange->s1;
		}
	      
	      break;
	    }
	  else if (nearestIdx && !gapsReversedWrtDisplay)
	    {
	      /* Remember the end of the current range (which is the result we will return
	       * if qIdx lies after this range but before the next - unless gaps are reversed
	       * with respect to the display). */
	      *nearestIdx = curRange->s2;
	    }
	}
    }
  
  return result;
}


/* Return the reading frame of the ref sequence that the given MSP is a match against */
int mspGetRefFrame(const MSP const *msp, const BlxSeqType seqType)
{
  if (seqType == BLXSEQ_DNA)
    {
      /* Ignore the frame in  the msp. For DNA matches we only have one frame on each strand. */
      return 1;
    }
  else
    {
      return atoi(&msp->qframe[2]);
    }
}


/* Return the strand of the ref sequence that the given MSP is a match against */
Strand mspGetRefStrand(const MSP const *msp)
{
  return msp->qframe[1] == '+' ? FORWARD_STRAND : REVERSE_STRAND;
}


/* Return the strand of the subject sequence that the given MSP is a match on */
Strand mspGetSubjectStrand(const MSP const *msp)
{
  return msp->sframe[1] == '+' ? FORWARD_STRAND : REVERSE_STRAND;
}


/* Add the given MSP to the given hash table */
void addMspToHashTable(GHashTable *hashTable, MSP *msp, char *hashKey)
{
  /* See if this MSP's sequence already exists in the hash table */
  SubjectSequence *subjectSeq = (SubjectSequence*)g_hash_table_lookup(hashTable, hashKey);
  
  if (subjectSeq)
    {
      /* Append the MSP to the existing sequence's list */
      subjectSeq->mspList = g_list_append(subjectSeq->mspList, msp);
    }
  else
    {
      /* Create a new SubjectSequence struct containing this MSP, and add it to the hash table */
      subjectSeq = g_malloc(sizeof(SubjectSequence));
      subjectSeq->seqName = msp->sname;
      subjectSeq->mspList = g_list_append(NULL, msp);
      
      g_hash_table_insert(hashTable, hashKey, subjectSeq);  
    }
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


/* Utility to pop up a simple dialog with the given title and text, with just an "OK" button. */
void showMessageDialog(const char *title,  
		       const char *messageText,
		       GtkWidget *parent,
		       const int initWidth,
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
  
  /* Create the text buffer and copy the text in */
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
  
  
  /* Put the text view in a scrollable window, and add it to the dialog */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(scrollWin), textView);
  
  GtkWidget *contentArea = GTK_DIALOG(dialog)->vbox;
  gtk_box_pack_start(GTK_BOX(contentArea), scrollWin, TRUE, TRUE, 0);
  
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  
  /* Set the initial height based on the number of lines (but don't make it bigger than the parent window) */
  PangoContext *context = gtk_widget_get_pango_context(textView);
  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
  gint charHeight = (pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / PANGO_SCALE;
  
  int initHeight = (gtk_text_buffer_get_line_count(textBuffer) * charHeight) + 100; /* fudge to allow space for buttons */
  initHeight = min(maxHeight, initHeight);
  gtk_window_set_default_size(GTK_WINDOW(dialog), initWidth, initHeight);
  
  pango_font_metrics_unref(metrics);

  /* Ensure dialog is destroyed when user responds */
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


