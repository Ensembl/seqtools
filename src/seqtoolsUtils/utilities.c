/*  File: utilities.c
 *  Author: Gemma Barson, 2010-01-05
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
 * Description: See utilities.h
 *----------------------------------------------------------------------------
 */

#include <seqtoolsUtils/utilities.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/utsname.h>


#define SEQTOOLS_TOOLBAR_NAME	"SeqtoolsToolbarName"
#define DEFAULT_PFETCH_WINDOW_WIDTH_CHARS	      85
#define DEFAULT_PFETCH_WINDOW_HEIGHT_FRACTION         0.7


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



/* Globals */
gboolean g_printCachedOnly = FALSE;
PrintScaleType g_printScaleType = PRINT_FIT_BOTH;


/* A struct used to record warning/error messages */
typedef struct _BlxMessage
  {
    char *text;
    time_t timestamp;
    GLogLevelFlags logLevel;
  } BlxMessage;



static CallbackData*   widgetGetCallbackData(GtkWidget *widget);
static void            displayMessageAsPopup(const gchar *message, GLogLevelFlags log_level, const char *titlePrefix, GtkWindow *parent, GtkStatusbar *statusBar);
static void            displayMessageAsList(GSList *messageList, const char *titlePrefix, const gboolean bringToFront);
static BlxMessage*     createBlxMessage(const char *text, const GLogLevelFlags logLevel);
static void            destroyBlxMessage(BlxMessage **blxMessage);
static GSList**        getMessageList();
void                   destroyMessageList();
static char*           blxMessageGetDisplayText(const BlxMessage *msg, const gboolean incTimestamp);
static void            addMessagesToBuffer(GSList *messageList, GtkTextBuffer *textBuffer, GtkTextTag *normalTag, GtkTextTag *highlightTag);
static gboolean        getUseScrolledMessages();
static void            setUseScrolledMessages(const gboolean newValue);
static gboolean        onSetUseScrolledMessages(GtkWidget *button, const gint responseId, gpointer data);
static gboolean        onSetUsePopupMessages(GtkWidget *button, const gint responseId, gpointer data);
static char*           getDialogIcon(GLogLevelFlags log_level);
static void            printMessageToStatusbar(const gchar *message, gpointer data);


/***********************************************************
 *		         Widgets			   * 
 ***********************************************************/


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


/* Recursively call the given function on this widget and any child widgets */
void callFuncOnAllChildWidgets(GtkWidget *widget, gpointer data)
{
  GtkCallback func = (GtkCallback)data;
  func(widget, NULL);
  
  if (GTK_IS_CONTAINER(widget))
    gtk_container_foreach(GTK_CONTAINER(widget), callFuncOnAllChildWidgets, data);
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
      g_signal_connect(G_OBJECT(label), "expose-event", G_CALLBACK(onExposePrintable), NULL);
    }

  return label;
}



/* Expose-event handler for labels/entrys that are required to be shown during printing. */
gboolean onExposePrintable(GtkWidget *widget, GdkEventExpose *event, gpointer callbackData)
{
  if (!widget || !(GTK_IS_LABEL(widget) || GTK_IS_ENTRY(widget)))
    {
      g_critical("Invalid widget type for printable label/text-entry [%p]\n", widget);
    }

  /* Get the window. (Labels don't have their own so get the parent's window) */
  GtkWidget *parent = gtk_widget_get_parent(widget);
  
  /* Only widgets that have a pixmap set will be shown in print output */
  GdkDrawable *drawable =  gdk_pixmap_new(parent->window, widget->allocation.width, widget->allocation.height, -1);
  gdk_drawable_set_colormap(drawable, gdk_colormap_get_system());
  widgetSetDrawable(widget, drawable);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GtkStyle *style = gtk_widget_get_style(widget);
  GdkColor *bgColor = &style->bg[GTK_STATE_NORMAL];
  gdk_gc_set_foreground(gc, bgColor);
  gdk_draw_rectangle(drawable, gc, TRUE, 0, 0, widget->allocation.width, widget->allocation.height);

  PangoLayout *layout = NULL;
  
  if (GTK_IS_LABEL(widget))
    {
      layout = gtk_label_get_layout(GTK_LABEL(widget));
    }
  else
    {
      layout = gtk_entry_get_layout(GTK_ENTRY(widget));
    }
  
  if (layout)
    {
      gtk_paint_layout(widget->style, drawable, GTK_STATE_NORMAL, TRUE, NULL, widget, NULL, 0, 0, layout);
    }
  
  return FALSE; /* let the default handler continue */
}


/***********************************************************
 *		         Ranges/values			   * 
 ***********************************************************/

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

/* Utility to set the given range to the given length, centred on the given coord */
void centreRangeOnCoord(IntRange *range, const int coord, const int length)
{
  range->min = coord - (length / 2);
  range->max = range->min + length;
}


/* Simple utility to determine whether the given value is within the given range.
 * Returns false if the given value is an unset int or the given range is null */
gboolean valueWithinRange(const int value, const IntRange const *range)
{
  return (range != NULL && value >= range->min && value <= range->max);
}


/* Return true if two IntRanges overlap */
gboolean rangesOverlap(const IntRange const *range1, const IntRange const *range2)
{
  return (range1->min <= range2->max && range1->max >= range2->min);
}

/* Return true if two IntRanges are equal */
gboolean rangesEqual(const IntRange const *range1, const IntRange const *range2)
{
  return (range1->min == range2->min && range1->max == range2->max);
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


/* Utility to bounds-limit the first range to within the second. If maintainLen is true, maintains
 * the length of the range if possible by shifting the range. */
void boundsLimitRange(IntRange *range, const IntRange const *limit, const gboolean maintainLen)
{
  const int len = getRangeLength(range);
  
  if (range->min < limit->min)
    {
      range->min = limit->min;
      
      if (maintainLen)
        {
          range->max = range->min + len;
        }
    }
  
  if (range->max > limit->max)
    {
      range->max = limit->max;
      
      if (maintainLen)
        {
          range->min = range->max - len;

          /* If limit is shorter than range, we'll have gone lower than the min again */
          boundsLimitValue(&range->min, limit);
       }
    }
  
  
  boundsLimitValue(&range->max, limit);
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


/***********************************************************
 *		         Misc                              * 
 ***********************************************************/

/* Determine (or give our best guess) the sequence type of a sequence, based on the characters it
 * contains. Returns BLXSEQ_DNA for nucleotide sequences or BLXSEQ_PEPTIDE for peptide sequences. */
BlxSeqType determineSeqType(char *seq)
{
    char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
    char *primenuc    = "ACGTUN";
    char *protonly    = "EFIPQZ";

    /* Simplified version of Sean Eddy's */

    int  pos;
    char c;
    int  po = 0;			/* count of protein-only */
    int  nt = 0;			/* count of t's */
    int  nu = 0;			/* count of u's */
    int  na = 0;			/* count of nucleotides */
    int  aa = 0;			/* count of amino acids */
    int  no = 0;			/* count of others */
  
    /* Look at the first 300 characters
     */
    for (pos = 0; seq[pos] && pos < 300; pos++)
    {
	c = toupper(seq[pos]);

	if (strchr(protonly, c)) 
	    po++;
	else if (strchr(primenuc, c)) {
	    na++;
	    if (c == 'T') nt++;
	    else if (c == 'U') nu++;
	}
	else if (strchr(aminos, c)) 
	    aa++;
	else if (isalpha(c)) 
	    no++;
    }

    if (po > 0) return BLXSEQ_PEPTIDE;
    else if (na > aa) return BLXSEQ_DNA;
    else return BLXSEQ_PEPTIDE;
} /* determineSeqType */


/* Use this routine to add -install in programs that use Dotter */
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


char* getSystemErrorText()
{
  char *result = NULL;

#ifdef SUN
  /* horrible hack for Sunos/Macs(?) which are not standard C compliant */
  result = blxprintf(SYSERR_FORMAT, errno, sys_errlist[errno]) ;
#elif defined(MACINTOSH)
  result = blxprintf(SYSERR_FORMAT, errno) ;
#else
  result = blxprintf(SYSERR_FORMAT, errno, strerror(errno)) ;
#endif

  return result;
}
  
  
/* Get the BlxStyle with the given name. Returns null and sets the error if 
 * we expected to find the name but didn't. The special-case input string "."
 * means "value not set", so in this case we return NULL and don't set the error. */
BlxStyle* getBlxStyle(const char *styleName, GSList *styles, GError **error)
{
  BlxStyle *result = NULL;

  if (styleName && strcmp(styleName, ".") != 0)
    {
      GSList *item = styles;
      
      for ( ; item; item = item->next)
	{
	  BlxStyle *style = (BlxStyle*)(item->data);
	  if (style && stringsEqual(style->styleName, styleName, FALSE))
	    {
	      result = style;
	      break;
	    }
	}
	
      if (!result)
	{
          /* Only report the same missing style once */
          static GSList *seenStyles = NULL;
          GSList *item = seenStyles;
          gboolean seen = FALSE;
          for ( ; item && !seen; item = item->next)
            {
              seen = stringsEqual((char*)(item->data), styleName, FALSE);
            }
          
          if (!seen)
            {
              seenStyles = g_slist_append(seenStyles, g_strdup(styleName));
              g_set_error(error, SEQTOOLS_ERROR, SEQTOOLS_ERROR_NO_STYLE, "Requested style '%s' does not exist.\n", styleName);
            }
	}
    }

  return result;
}

/***********************************************************
 *		         Colors				   * 
 ***********************************************************/

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
      g_set_error(error, SEQTOOLS_ERROR, 1, "Error parsing color string '%s'\n", colorStr);
    }
    
  return ok;
}


/* Increase the brightness of a color by the given factor */
static void increaseColorBrightness(const GdkColor const *origColor, const double factor, GdkColor *result)
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
      g_warning("Error calculating highlight color; using original color instead (RGB=%d %d %d).", origColor->red, origColor->green, origColor->blue);
      result->pixel = origColor->pixel;
      result->red = origColor->red;
      result->green = origColor->green;
      result->blue = origColor->blue;
    }
}


/* Reduce the brightness of a color by the given factor (e.g. 0.5 for half brightness) */
static void reduceColorBrightness(const GdkColor const *origColor, const double factor, GdkColor *result)
{
  result->pixel = 0;
  result->red = origColor->red * factor;
  result->green = origColor->green * factor;
  result->blue = origColor->blue * factor;
  
  gboolean failures[1];
  gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), result, 1, TRUE, TRUE, failures);
  
  if (numFailures > 0)
    {
      g_warning("Error calculating highlight color; using original color instead (RGB=%d %d %d).", origColor->red, origColor->green, origColor->blue);
      result->pixel = origColor->pixel;
      result->red = origColor->red;
      result->green = origColor->green;
      result->blue = origColor->blue;
    }
}


/* Increase/decrease color brightness by given factor */
void adjustColorBrightness(const GdkColor const *origColor, const double factor, GdkColor *result)
{
  if (factor > 1)
    increaseColorBrightness(origColor, factor, result);
  else if (factor < 1)
    reduceColorBrightness(origColor, factor, result);
}


/* Utility to take a gdk color and return a slightly darker shade of it, for higlighting
 * things of that color when they are selected. */
void getSelectionColor(const GdkColor const *origColor, GdkColor *result)
{
  const double factor = 0.8;
  adjustColorBrightness(origColor, factor, result);
}


/* Utility to take a gdk color and return a greyscale version of it */
void convertToGrayscale(const GdkColor const *origColor, GdkColor *result)
{
  result->red = result->green = result->blue =
   ((11 * origColor->red) + (16 * origColor->green) + (5 * origColor->blue)) / 32;
   
  gboolean failures[1];
  gint numFailures = gdk_colormap_alloc_colors(gdk_colormap_get_system(), result, 1, TRUE, TRUE, failures);
  
  if (numFailures > 0)
    {
      g_warning("Error calculating greyscale color (RGB=%d %d %d).", origColor->red, origColor->green, origColor->blue);
      
      /* Set it to black, for want of something better to do... */
      result->pixel = 0;
      result->red = 0;
      result->green = 0;
      result->blue = 0;
    }

}


/* Utility to take a gdk color and return a much darker shade of it, for use as a drop-shadow */
void getDropShadowColor(const GdkColor const *origColor, GdkColor *result)
{
  const double factor = 0.3;
  adjustColorBrightness(origColor, factor, result);
}


/* Return a pointer to the BlxColor with the given id */
BlxColor* getBlxColor(GArray *defaultColors, const int colorId)
{
  BlxColor *result = &g_array_index(defaultColors, BlxColor, colorId);

  if (!result)
    {
      g_warning("Color ID %d not found in the default colors array.\n", colorId);
    }
  
  return result;
}


/* Lookup the BlxColor with the given id in the given hash table and return a 
 * pointer to the gdk color with the given properties */
GdkColor* getGdkColor(const int colorId, GArray *defaultColors, const gboolean selected, const gboolean usePrintColors)
{
  GdkColor *result = NULL;
  
  BlxColor *blxColor = getBlxColor(defaultColors, colorId);
  
  if (blxColor)
    {
      if (usePrintColors)
	{
	  result = selected ? &blxColor->printSelected : &blxColor->print;
	}
      else
	{
	  result = selected ? &blxColor->selected : &blxColor->normal;
	}
    }
  else
    {
      g_warning("Error: requested invalid color ID %d", colorId);
    }
  
  return result;
}


/* This basically does what gdk_color_to_string does, but that function isn't available in
 * older versions of GDK... */
char* convertColorToString(GdkColor *color)
{
  //  char *result = gdk_color_to_string(&widget->style->bg[GTK_STATE_NORMAL]);
  //  return result;
  
  /* Need to make sure the color is allocated (otherwise the 'pixel' field may be zero) */
  gboolean failures[1];
  gdk_colormap_alloc_colors(gdk_colormap_get_system(), color, 1, TRUE, TRUE, failures);
  
  const int hexLen = 8; /* to fit a string of the format '#ffffff', plus the terminating character */
  
  char *result = g_malloc(sizeof(char) * hexLen);
  sprintf(result, "#%x", color->pixel);
  
  return result;
}


/* Free all memory used by a BlxColor */
void destroyBlxColor(BlxColor *blxColor)
{
  if (blxColor)
    {
      g_free(blxColor->name);
      g_free(blxColor->desc);
    }
}

/* Create a BlxColor. The result should be destroyed with destroyBlxColor. Returns
 * NULL if there was any problem. If "selected" state colors are not passed, this function
 * calculates a slightly darker shade of the normal color to use for when it is selected.
 * Inserts the new color into the given list */
void createBlxColor(GArray *defaultColors,
		    int colorId,
		    const char *name, 
		    const char *desc, 
		    const char *normalCol, 
		    const char *printCol,
		    const char *normalColSelected,
		    const char *printColSelected)
{
  BlxColor *result = &g_array_index(defaultColors, BlxColor, colorId);
  
  result->transparent = FALSE;
  GError *error = NULL;
  
  if (!normalCol)
    {
      result->transparent = TRUE;
    }
  else if (getColorFromString(normalCol, &result->normal, &error)) 
    {
      /* find a "selected" version of it, if not passed one */
      if (normalColSelected)
	{
	  if (!getColorFromString(normalColSelected, &result->selected, &error))
	    {
	      prefixError(error, "Error getting 'selected' color: using normal color instead. ");
	      reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
	      result->selected = result->normal;
	    }
	}
    else
	{
	  getSelectionColor(&result->normal, &result->selected); 
	}
    
    /* Parse the print color */
    if (getColorFromString(printCol, &result->print, &error))
      {
	/* find a "selected" version of it, if not passed one */
	if (printColSelected)
	  {
	    if (!getColorFromString(printColSelected, &result->printSelected, &error))
	      {
		prefixError(error, "Error getting print color for selected items: using normal print color instead. ");
		reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
		result->printSelected = result->print;
	      }
	  }
	else
	  {
	    getSelectionColor(&result->print, &result->printSelected); 
	  }
	}
      else
	{
	/* Error parsing the print color: use the normal color but give a warning */
	prefixError(error, "Error getting print colors: using normal colors instead. ");
	reportAndClearIfError(&error, G_LOG_LEVEL_MESSAGE);
	result->print = result->normal;
	result->printSelected = result->selected;
	}
    }
  else
    {
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      g_free(result);
      result = NULL;
    }
  
  if (result)
    {
      /* Set the other properties */
      result->name = g_strdup(name);
      result->desc = g_strdup(desc);
    }
}


/* Calculate how many pixels wide a base is */
gdouble pixelsPerBase(const gint displayWidth, const IntRange const *displayRange)
{
  gdouble displayLen = (gdouble)(displayRange->max - displayRange->min) + 1;
  return ((gdouble)displayWidth / displayLen);
}


/* Convert a base index to a coord within the given rectangle - the x coord if 
 * horizontal is true, the y coord otherwise. Returns the number of pixels 
 * to the coord from 0 at the left edge (or top edge, if vertical), including 
 * any pixels up to the left/top edge of the rectangle. 
 * 
 * 'displayRev' should be passed as true if the display is reversed (i.e. if 
 * values decrease from left-to-right/top-to-bottom). 
 * 
 * The result is clipped to lie within the rectangle if 'clip' is true */
gdouble convertBaseIdxToRectPos(const gint dnaIdx, 
 			        const GdkRectangle const *rect, 
 			        const IntRange const *dnaDispRange,
 			        const gboolean horizontal,
 			        const gboolean displayRev,
 			        const gboolean clip)
{
  int baseIdx = invertCoord(dnaIdx, dnaDispRange, displayRev);
  
  gdouble numBasesFromEdge = (gdouble)(baseIdx - dnaDispRange->min); /* 0-based index from edge */
  
  if (clip && numBasesFromEdge < 0)
    {
      numBasesFromEdge = 0;
    }
  
  const int rectLength = horizontal ? rect->width : rect->height;
  gdouble pixelsFromEdge = numBasesFromEdge * pixelsPerBase(rectLength, dnaDispRange);
  
  const int rectStart = horizontal ? rect->x : rect->y;
  gdouble result = (gdouble)rectStart + pixelsFromEdge;
  
  if (clip && result > rectStart + rectLength)
    {
      result = rectStart + rectLength;
    }
  
  return result;
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


/* Convert the given range of display coords to dna coords */
void convertDisplayRangeToDnaRange(const IntRange const * displayRange, 
                                   const BlxSeqType displaySeqType,
                                   const int numFrames,
                                   const gboolean displayRev,
                                   const IntRange const *refSeqRange,
                                   IntRange *result)
{
  const int q1 = convertDisplayIdxToDnaIdx(displayRange->min, displaySeqType, 1, 1, numFrames, displayRev, refSeqRange); /* 1st base in 1st reading frame */
  const int q2 = convertDisplayIdxToDnaIdx(displayRange->max, displaySeqType, numFrames, numFrames, numFrames, displayRev, refSeqRange); /* last base in last frame */
  intrangeSetValues(result, q1, q2);
}


/* Get the reading frame (i.e. modulus 3) for the coord. Inverts the frame if displayRev is true. */
static int getCoordReadingFrame(const int coord, const int numFrames, const gboolean displayRev)
{
  int result = coord % numFrames;
  
  if (result < 1)
    result += numFrames;
    
  if (displayRev)
    result = numFrames - result + 1;
    
  return result;
}


/* Given an index into the displayed sequence, a reading frame, and the base number within that
 * reading frame, return the index into the DNA sequence that will give the equivalent DNA base.
 * If the display sequence is a peptide sequence, it will convert the coord to a DNA coord. If the
 * display is reversed, the display coord will be inverted. */
int convertDisplayIdxToDnaIdx(const int displayIdx, 
			      const BlxSeqType srcSeqType,
			      const int frame, 
			      const int baseNum, 
			      const int numFrames,
			      const gboolean displayRev,
			      const IntRange const *refSeqRange)
{
  int dnaIdx = displayIdx;
  const int base = baseNum;
  
  if (srcSeqType == BLXSEQ_PEPTIDE)
    {
      /* Convert the input peptide coord to a dna coord */
      dnaIdx = (displayIdx * numFrames) - numFrames + frame + base - 1;
    }
  
  if (displayRev)
    {
      int origReadingFrame = getCoordReadingFrame(dnaIdx, numFrames, displayRev);
      
      /* If the display is reversed, we need to invert the result. For example, if the 
       * result is index '2' out of the range '12345', then we convert it to '4', which is the
       * equivalent position in the range '54321'. */
      dnaIdx = refSeqRange->max - dnaIdx + refSeqRange->min;

      if (srcSeqType == BLXSEQ_PEPTIDE)
        {
          /* Undo any offset that was added in convertDnaIdxToDisplayIdx to make the reading frame
           * of the inverted coord the same as the base number of the forward coord. */
          int newReadingFrame = getCoordReadingFrame(dnaIdx, numFrames, FALSE);
          int offset = origReadingFrame - newReadingFrame;
          
          if (offset < 0) 
            offset += numFrames;

          dnaIdx += offset;
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
  
  const gboolean peptides = (displaySeqType == BLXSEQ_PEPTIDE);
  int base = peptides ? (dnaIdx - frame + 1) % numFrames : 1;

  if (base < 1)
    base += numFrames;
  
  /* If the display is reversed (i.e. showing increasing values from right-to-left),
   * invert the index (i.e. as if it were the normal left-to-right index for this
   * same position - for example, if the index is '4' out of the range '54321', convert
   * it to '2', which is the equivalent position in the range '12345'. */
  if (displayRev)
    {
      /* Invert the coord and base */
      displayIdx = dnaIdxRange->max - dnaIdx + dnaIdxRange->min;
      base = numFrames - base + 1;
  
      if (peptides)
    {
          /* When we do a mod-3 of the inverted coord, we want it to give the same reading frame
           * as the mod-3 of the original coord */
          const int origReadingFrame = getCoordReadingFrame(dnaIdx, numFrames, displayRev);
          const int newReadingFrame = getCoordReadingFrame(displayIdx, numFrames, FALSE);
          int offset = origReadingFrame - newReadingFrame;
          if (offset < 0) 
            offset += numFrames;
          displayIdx += offset;
        }
    }

  /* Convert from nucleotide to peptide coords */
  if (peptides)
    displayIdx = ceil((gdouble)(displayIdx - frame + 1) / (gdouble)numFrames);
  
  /* Set the output arg, if requested */
      if (baseNum)
    *baseNum = base;
      
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


/* For the given BlxColor, return a pointer to the GdkColor that meets the given criteria */
const GdkColor *blxColorGetColor(const BlxColor *blxColor, const gboolean selected, const gboolean usePrintColors)
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


/* Converts the given double to a string. The result must be free'd with g_free. Numdp specifies the
 * number of decimal places to show */
char* convertDoubleToString(const gdouble value, const int numDp)
{
  char format[numDp + 5];
  sprintf(format, "%%1.%df", numDp);
  
  char result[numDigitsInInt((int)value) + numDp + 2]; /* the +2 includes the and terminating null */
  
  sprintf(result, format, value);
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


/* Get the child of the given widget that has the given name (which could be the given 
 * widget itself.) Assumes there is only one, so returns the first one found. */
GtkWidget* getNamedChildWidget(GtkWidget *widget, const gchar *searchName)
{
  GtkWidget *result = NULL;
  const gchar *widgetName = gtk_widget_get_name(widget);
 
  if (stringsEqual(widgetName, searchName, TRUE))
    {
      result = widget;
    }
  else if (GTK_IS_CONTAINER(widget))
    {
      /* recurse over children until we find the detail view (assumes there is only one!) */
      GList *child = gtk_container_get_children(GTK_CONTAINER(widget));
      
      for ( ; child && !result; child = child->next)
	{
	  GtkWidget *childWidget = GTK_WIDGET(child->data);
	  result = getNamedChildWidget(childWidget, searchName);
	}
    }
    
  return result;
}


/* Send a string to a file, "protecting" it by placing quotes around it and escaping quotes
 * inside it. */
void stringProtect(FILE *file, const char *string)
{
  const char *cp;
 
  fputc(' ', file);
  fputc('"', file);
  if (string)
    for(cp = string; *cp; ++cp)
      {
	if (*cp == '"' || *cp == '$')
	  fputc('$', file);
	fputc(*cp, file);
      }
  fputc('"', file);
  
}


/* Read a protected string */
char *stringUnprotect(char **textp, char *target)
{
  char *cp, *cpd;
  int count = 0;

 redo:
  cp = *textp;
  cpd = target;
  
  while (*cp)
    {
      if (*cp == '"')
	{
	  cp++;						    /* skip quote */
	  break ;
	}
      else
	cp++ ;
    }

  while (*cp != '"' && *cp)
    {
      if (*cp == '$')
	cp++;

      if (cpd)
	*cpd++ = *cp;
      else
	count++;

      cp++;
    }
  
  if (!target)
    {
      target = g_malloc(count+1);
      goto redo;
    }
  
  *cp = '\0' ;
  *textp = cp+1; /* skip quote */

  return target;
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
      CallbackData *callbackData = widgetGetCallbackData(widget);
      
      if (!callbackData)
        callbackData = g_malloc(sizeof *callbackData);
      
      callbackData->func = func;
      callbackData->data = data;
      
      g_object_set_data(G_OBJECT(widget), "callbackData", callbackData);
      g_signal_connect(G_OBJECT(widget), "destroy", G_CALLBACK(onDestroyCustomWidget), NULL);
    }
}

/* Get the CallbackData struct for this widget, if it has one, and call it. Returns false if
 * there was an error */
static gboolean widgetCallCallback(GtkWidget *widget, const gint responseId)
{
  gboolean result = TRUE;
  CallbackData *callbackData = widgetGetCallbackData(widget);
  
  if (callbackData && callbackData->func)
    {
      result = callbackData->func(widget, responseId, callbackData->data);
    }
  
  return result;
}

/* Call the stored CallbackData callbacks (if any exist) for this widget and
 * all of its children. The response ID is passed as the user data. Returns false if
 * any callback returns an error. */
gboolean widgetCallAllCallbacks(GtkWidget *widget, gpointer data)
{
  gboolean result = TRUE;
  
  if (widget && GTK_IS_WIDGET(widget))
    {
      gint responseId = GPOINTER_TO_INT(data);
      
      if (!widgetCallCallback(widget, responseId))
        {
          result = FALSE;
        }
  
      if (GTK_IS_CONTAINER(widget))
	{
          GList *childItem = gtk_container_get_children(GTK_CONTAINER(widget));
          
          for ( ; childItem; childItem = childItem->next)
            {
              GtkWidget *child = GTK_WIDGET(childItem->data);
              
              if (!widgetCallAllCallbacks(child, data))
                {
                  result = FALSE;
                }
            }
	}
    }
  
  return result;
}


/* Generic callback to call all user-specified callbacks for all child widgets of
 * the given dialog if ACCEPT or APPLY responses received. Also closes the dialog 
 * if ACCEPT or REJECT responses received.
 * The user data must be a boolean (which has been converted to a pointer using
 * GINT_TO_POINTER), that is true if this dialog is persistent or false 
 * otherwise. If it is persistent, the dialog will be hidden rather than being 
 * destroyed. */
void onResponseDialog(GtkDialog *dialog, gint responseId, gpointer data)
{
  gboolean destroy = TRUE;
  
  switch (responseId)
  {
    case GTK_RESPONSE_ACCEPT:
      /* Destroy if successful */
      destroy = widgetCallAllCallbacks(GTK_WIDGET(dialog), GINT_TO_POINTER(responseId));
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


/* Gets the dialog widget for the given dialog id. Returns null if the widget has not
 * been created yet. */
GtkWidget* getPersistentDialog(GtkWidget* dialogList[], const int dialogId)
{
  GtkWidget *result = NULL;
  
  if (dialogList[dialogId])
    {
      result = dialogList[dialogId];
    }
  
  return result;
}

/* Add a newly-created dialog to the list of persistent dialogs. The dialog should not
 * exist in the list yet. */
void addPersistentDialog(GtkWidget* dialogList[], const int dialogId, GtkWidget *widget)
{
  if (dialogId == 0)
    {
      g_warning("Code error: cannot add a dialog with ID %d. Dialog will not be persistent.\n", dialogId);
    }
  else
    {
      if (dialogList[dialogId])
        {
          g_warning("Creating a dialog that already exists. Old dialog will be destroyed. Dialog ID=%d.\n", dialogId);
          gtk_widget_destroy(dialogList[dialogId]);
          dialogList[dialogId] = NULL;
        }
      
      dialogList[dialogId] = widget;
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


/* Create and cache a blank pixmap drawable in the given widget. The size of the
 * resultant pixmap is the same as the widget's allocation size. */
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
  
  g_object_unref(gc);
  
  return drawable;
}


/* Like createBlankPixmap, but allows the caller to specify a different size to
 * the widget allocation size */
GdkDrawable* createBlankSizedPixmap(GtkWidget *widget, GdkDrawable *window, const int width, const int height)
{
  GdkDrawable *bitmap = gdk_pixmap_new(window, width, height, -1);
  gdk_drawable_set_colormap(bitmap, gdk_colormap_get_system());
  widgetSetDrawable(widget, bitmap); /* deletes the old drawable, if there is one */
  
  /* Paint a blank rectangle for the background, the same color as the widget's background */
  GdkGC *gc = gdk_gc_new(bitmap);
  
  GdkColor *bgColor = &widget->style->bg[GTK_STATE_NORMAL];
  gdk_gc_set_foreground(gc, bgColor);
  gdk_draw_rectangle(bitmap, gc, TRUE, 0, 0, width, height);
  
  g_object_unref(gc);
  
  return bitmap;
}


/* Utility to pop up a simple confirmation dialog box with the given title and text, 
 * with just an "OK" and "Cancel" button. Blocks until the user responds, and returns
 * their response ID. */
gint runConfirmationBox(GtkWidget *widget, char *title, char *messageText)
{
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  GTK_WINDOW(widget), 
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
      if (log_level & G_LOG_LEVEL_ERROR)
        g_error("%s", (*error)->message);
      else if (log_level & G_LOG_LEVEL_CRITICAL)
        g_critical("%s", (*error)->message);
      else if (log_level & G_LOG_LEVEL_WARNING)
        g_warning("%s", (*error)->message);
      else if (log_level & G_LOG_LEVEL_DEBUG)
        g_debug("%s", (*error)->message);
      else
        g_message("%s", (*error)->message); /* message or info */
            
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
  const int width = max(minWidth, rect->width);

  cairo_t *cr = gdk_cairo_create(drawable);
  
  gdk_cairo_set_source_color(cr, color);
  cairo_rectangle(cr, rect->x, rect->y, width, rect->height);
  cairo_clip(cr);
  cairo_paint_with_alpha(cr, 0.2);
  
  cairo_destroy(cr);
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



/* Copy a segment of the given sequence into a new string. The result must be
 * free'd with g_free by the caller. The given indices must be 0-based. */
static gchar *copySeqSegment(const char const *inputSeq, const int idx1, const int idx2)
{

  const int minIdx = min(idx1, idx2);
  const int maxIdx = max(idx1, idx2);
  
  const int segmentLen = maxIdx - minIdx + 1;

  gchar *segment = g_malloc(sizeof(gchar) * segmentLen + 1);

  
  strncpy(segment, inputSeq + minIdx, segmentLen);
  segment[segmentLen] = '\0';

  return segment;
}



/* Copy a segment of the given sequence (which is always the DNA sequence and
 * always the forward strand).
 *
 * The result is complemented if the reverse strand is requested, but only if
 * the allowComplement flag allows it.
 * 
 * The result is translated to a peptide sequence if the destination seq type 
 * is peptide.
 *
 * The result is reversed if the reverseResult flag is true (regardless of the
 * strand requested - this is because the caller often wants the result in the
 * opposite direction to that indicated by the strand, because the display may be
 * reversed, so we leave it up to the caller to decide).
 */
gchar *getSequenceSegment(const char const *dnaSequence,
                          IntRange *qRangeIn,                  /* the range to extract, in nucleotide coords */
			  const BlxStrand strand,
			  const BlxSeqType srcSeqType,        /* whether input sequence is nucleotide or peptide */
			  const BlxSeqType destSeqType,       /* whether result sequence should be nucleotide or peptide */
			  const int frame,
			  const int numFrames, 
			  const IntRange const *refSeqRange,
			  const BlxBlastMode blastMode,
			  char **geneticCode,
			  const gboolean displayRev,
			  const gboolean reverseResult,
			  const gboolean allowComplement,
			  GError **error)
{
  gchar *result = NULL;
  GError *tmpError = NULL;
  
  if (destSeqType == BLXSEQ_DNA && srcSeqType == BLXSEQ_PEPTIDE)
    {
      /* We shouldn't try to convert from peptide to nucleotide. If we get here it's a coding error */
      g_set_error(error, SEQTOOLS_ERROR, SEQTOOLS_ERROR_SEQ_SEGMENT, "Error: requested conversion of peptide sequence to nucleotide sequence.\n");
      return result;
    }
  
  IntRange qRange = {qRangeIn->min, qRangeIn->max};
  
  if (qRange.min < refSeqRange->min || qRange.max > refSeqRange->max)
    {
      /* We might request up to 3 bases beyond the end of the range if we want to 
       * show a partial triplet at the start/end. (It's a bit tricky for the caller to
       * specify the exact bases they want here so we allow it but just limit it to 
       * the actual range so that they can't index beyond the end of the range.) Any
       * further out than one triplet is probably indicative of an error, so give a warning. */
      if (qRange.min < refSeqRange->min - (numFrames + 1) || qRange.max > refSeqRange->max + (numFrames + 1))
	{
          g_set_error(&tmpError, SEQTOOLS_ERROR, SEQTOOLS_ERROR_SEQ_SEGMENT, "Requested query sequence %d - %d out of available range: %d - %d.\n", qRange.min, qRange.max, refSeqRange->min, refSeqRange->max);
	}
      
      if (qRange.min < refSeqRange->min)
        {
          qRange.min = refSeqRange->min;
        }
      
      if (qRange.max > refSeqRange->max)
        {
          qRange.max = refSeqRange->max;
        }
    }
  
  /* Get 0-based indices into the sequence */
  const int idx1 = qRange.min - refSeqRange->min;
  const int idx2 = qRange.max - refSeqRange->min;

  /* Copy the required segment from the ref seq. Must pass 0-based indices into the sequence */
  result = copySeqSegment(dnaSequence, idx1, idx2);

  /* If a nucleotide seq, complement if this is the reverse strand (and if allowed). NB the old code didn't used to do
   * this for tblastn or blastp modes so I've kept that behaviour, but I'm not sure why it is that way. */
  if (srcSeqType == BLXSEQ_DNA && strand == BLXSTRAND_REVERSE && allowComplement && blastMode != BLXMODE_TBLASTN && blastMode != BLXMODE_BLASTP)
    {
      blxComplement(result);
      
      if (!result)
        {
          g_set_error(error, SEQTOOLS_ERROR, SEQTOOLS_ERROR_SEQ_SEGMENT, "Error getting sequence segment: Failed to complement the reference sequence for the range %d - %d.\n", qRange.min, qRange.max);
          g_free(result);
          return NULL;
        }
    }
  
  /* Reverse if requested. */
  if (reverseResult)
    {
      g_strreverse(result);
    }

  /* Translate if we're going from a nucleotide seq to a peptide seq */
  if (srcSeqType == BLXSEQ_DNA && destSeqType == BLXSEQ_PEPTIDE)
    {
      char *tmp = blxTranslate(result, geneticCode);
	  
      g_free(result); /* delete the original because it's no longer required */
      result = tmp;
	  
      if (!result)
        {
          g_set_error(error, SEQTOOLS_ERROR, SEQTOOLS_ERROR_SEQ_SEGMENT,
                      "Error getting the sequence segment: Failed to translate the DNA sequence for the reference sequence range %d - %d\n", 
                      qRange.min, qRange.max) ;
          
          return NULL;
        }
    }
  
  if (!result)
    {
      g_set_error(error, SEQTOOLS_ERROR, SEQTOOLS_ERROR_SEQ_SEGMENT, "Failed to find sequence segment for the range %d - %d\n", qRange.min, qRange.max);
    }
  
  if (tmpError && *error != NULL)
    {
      prefixError(*error, tmpError->message);
    }
  else if (tmpError)
    {
      g_propagate_error(error, tmpError);
    }

  return result;
}


/* Invert the given coord's position within the given range. Only invert if the bool is true; otherwise
 * returns the original coord */
int invertCoord(const int coord, const IntRange const *range, const gboolean invert)
{
  int result = invert ? range->max - coord + range->min : coord;
  return result;
}



/* Tries to return a fixed font from the list given in pref_families, returns
 * TRUE if it succeeded in finding a matching font, FALSE otherwise.
 * The list of preferred fonts is treated with most preferred first and least
 * preferred last.  The function will attempt to return the most preferred font
 * it finds.
 *
 * @param widget         Needed to get a context, ideally should be the widget you want to
 *                       either set a font in or find out about a font for.
 * @param pref_families  List of font families (as text strings).
 * @param points         Size of font in points.
 * @param weight         Weight of font (e.g. PANGO_WEIGHT_NORMAL)
 * @param font_out       If non-NULL, the font is returned.
 * @param desc_out       If non-NULL, the font description is returned.
 * @return               TRUE if font found, FALSE otherwise.
 */
const char* findFixedWidthFontFamily(GtkWidget *widget, GList *pref_families)
{
  /* Find all the font families available */
  PangoContext *context = gtk_widget_get_pango_context(widget) ;
  PangoFontFamily **families;
  gint n_families;
  pango_context_list_families(context, &families, &n_families) ;
  
  /* Loop through the available font families looking for one in our preferred list */
  gboolean found_most_preferred = FALSE;
  gint most_preferred = g_list_length(pref_families);
  PangoFontFamily *match_family = NULL;

  gint family;
  for (family = 0 ; (family < n_families && !found_most_preferred) ; family++)
    {
      const gchar *name = pango_font_family_get_name(families[family]) ;
      
      /* Look for this font family in our list of preferred families */
      GList *pref = g_list_first(pref_families) ;
      gint current = 1;
      
      while(pref)
	{
	  char *pref_font = (char *)pref->data ;
	  
	  if (g_ascii_strncasecmp(name, pref_font, strlen(pref_font)) == 0
#if GLIB_MAJOR_VERSION >= 1 && GLIB_MINOR_VERSION >= 4
	      && pango_font_family_is_monospace(families[family])
#endif
	      )
	    {
	      /* We prefer ones nearer the start of the list */
              if(current <= most_preferred)
                {
		  most_preferred = current;
		  match_family = families[family];

                  if(most_preferred == 1)
		    {
		      found_most_preferred = TRUE;
		    }
                }

	      break;
	    }
	  
	  pref = g_list_next(pref);
	  ++current;
	}
    }

  const char *result = NULL;
  if (match_family)
    {
      result = pango_font_family_get_name(match_family);
      g_debug("Using fixed-width font '%s'\n", result);
    }
  else
    {
      g_critical("Could not find a fixed-width font. Alignments may not be displayed correctly.\n");
    }
  
  return result;
}


/* We need a fixed-width font for displaying alignments. Find one from a 
 * list of possibilities. */
const char* findFixedWidthFont(GtkWidget *widget)
{
  GList *fixed_font_list = NULL ;

//  fixed_font_list = g_list_append(fixed_font_list, "inconsolata");
//  fixed_font_list = g_list_append(fixed_font_list, "consolas");
//  fixed_font_list = g_list_append(fixed_font_list, "droid sans mono");
//  fixed_font_list = g_list_append(fixed_font_list, "proggy");
//  fixed_font_list = g_list_append(fixed_font_list, "profont");
  fixed_font_list = g_list_append(fixed_font_list, "Andale mono");
  fixed_font_list = g_list_append(fixed_font_list, "Lucida sans typewriter");
  fixed_font_list = g_list_append(fixed_font_list, "deja vu sans mono");
  fixed_font_list = g_list_append(fixed_font_list, "Bitstream vera sans mono");
  fixed_font_list = g_list_append(fixed_font_list, "monaco");
  fixed_font_list = g_list_append(fixed_font_list, "Lucida console");
  fixed_font_list = g_list_append(fixed_font_list, "Courier 10 pitch");
  fixed_font_list = g_list_append(fixed_font_list, "Courier new");
  fixed_font_list = g_list_append(fixed_font_list, "Courier");
  fixed_font_list = g_list_append(fixed_font_list, "Monospace");
  fixed_font_list = g_list_append(fixed_font_list, "fixed");
  
  const char *fontFamily = findFixedWidthFontFamily(widget, fixed_font_list);
  g_list_free(fixed_font_list);
  
  return fontFamily;
}


/* Utility to get the character width and height of a given pango font */
void getFontCharSize(GtkWidget *widget, PangoFontDescription *fontDesc, gdouble *width, gdouble *height)
{
  PangoContext *context = gtk_widget_get_pango_context(widget);
  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
  
  if (height)
    {
      *height = (gdouble)(pango_font_metrics_get_ascent (metrics) + pango_font_metrics_get_descent (metrics)) / (gdouble)PANGO_SCALE;
    }
  
  if (width)
    {
      *width = (gdouble)pango_font_metrics_get_approximate_digit_width(metrics) / (gdouble)PANGO_SCALE;
    }
  
  pango_font_metrics_unref(metrics);
}


/***********************************************************
 *                      Memory
***********************************************************/

/* Create a handle. A handle is just a GSList so should be initialised to NULL */
BlxHandle handleCreate()
{
  BlxHandle handle = g_malloc(sizeof(BlxHandleStruct));
  handle->memoryList = NULL;
  return handle;
}

/* Utility to allocate memory of the given size with g_malloc. Returns the allocated memory and
 * adds a pointer to it to the given handle (list), so that all memory allocated to this list can be destroyed
 * with handleDestroy */
gpointer handleAlloc(BlxHandle *handle, size_t numBytes)
{
  gpointer result = g_malloc(numBytes);
  
  /* prepend so that we delete data in reverse order */
  (*handle)->memoryList = g_slist_prepend((*handle)->memoryList, result); 

  return result;
}

/* Utilitiy to call g_free on all elements in the given handle. Frees the handle too and sets it to null. */
void handleDestroy(BlxHandle *handle)
{
  DEBUG_ENTER("handleDestroy");

  GSList *listItem = (*handle)->memoryList;
  for ( ; listItem; listItem = listItem->next)
    {
      g_free(listItem->data);
      listItem->data = NULL;
    }
  
  g_slist_free((*handle)->memoryList);
  g_free(*handle);
  *handle = NULL;
  
  DEBUG_EXIT("handleDestroy returning ");
}


/***********************************************************
 *                      Toolbars
***********************************************************/

static void buttonAttach(GtkHandleBox *handlebox, GtkWidget *toolbar, gpointer data)
{
  gtk_widget_set_usize(toolbar, 1, -2);
}


static void buttonDetach(GtkHandleBox *handlebox, GtkWidget *toolbar, gpointer data)
{
  gtk_widget_set_usize(toolbar, -1, -2);
}


/* This creates a handle box where we can attach/detach a toolbar */
GtkWidget *createToolbarHandle()
{
  GtkWidget *handleBox = gtk_handle_box_new();

  /* this stops the toolbar forcing the size of the parent window */
  g_signal_connect(GTK_OBJECT(handleBox), "child-attached", G_CALLBACK(buttonAttach), NULL);
  g_signal_connect(GTK_OBJECT(handleBox), "child-detached", G_CALLBACK(buttonDetach), NULL);

  return handleBox;
}


/* Makes the given widget into a toolbar item on the given toolbar.
 * 'position' is the position; can be 0 to add to the start or -1 to add to the end. */
GtkToolItem* addToolbarWidget(GtkToolbar *toolbar, GtkWidget *widget, const int position)
{
  GtkToolItem *toolItem = gtk_tool_item_new();
  gtk_container_add(GTK_CONTAINER(toolItem), widget);
  gtk_toolbar_insert(toolbar, toolItem, position);
  
  gtk_tool_item_set_visible_horizontal(toolItem, TRUE);
  gtk_tool_item_set_visible_vertical(toolItem, TRUE);
  gtk_tool_item_set_is_important(toolItem, TRUE);
  
  return toolItem;
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


/***********************************************************
 *                     Message handlers
 ***********************************************************/

/* GLib doesn't have a convenience function to log a message with the 
 * G_LOG_LEVEL_INFO flag, so this provides us one. In seqtools, this is usd
 * to provide overall task progress messages (whereas we use g_message to 
 * provide more detailed sub-task progress messages). */
void g_message_info(char *formatStr, ...)
{
  va_list argp;
  va_start(argp, formatStr);
  
  g_logv(G_LOG_DOMAIN, G_LOG_LEVEL_INFO, formatStr, argp);
  
  va_end(argp);
}


/* Default handler for GLib log messages (e.g. from g_message() etc.) */
void defaultMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data)
{
  const char *warningText = "Warning: ";
  
  const int len = strlen(warningText) + 1;
  char prefixText[len];
  prefixText[len - 1] = '\0';
  
  if (log_level & G_LOG_LEVEL_WARNING)
    strcpy(prefixText, warningText);
  else
    prefixText[0] = '\0';
  
  /* print to console (indented by the same amount as debug output, if relevant) */
#ifdef DEBUG
  DEBUG_OUT("%s%s", prefixText, (char*)message);
#else
  if (log_level & G_LOG_LEVEL_INFO)
    fprintf(stderr, "%s%s", prefixText, (char*)message); /* info messages to to stderr */
  else
    printf("%s%s", prefixText, (char*)message);
#endif
  
  /* also display the message in the status bar (unless it's debug output) */
  if (!(log_level & G_LOG_LEVEL_DEBUG))
    {
      BlxMessageData *msgData = data ? (BlxMessageData*)data : NULL;
      GtkStatusbar *statusBar = msgData ? msgData->statusBar : NULL;
      printMessageToStatusbar(message, statusBar);
    }
}


/* Default handler for critical GLib log messages (i.e. from g_error and g_critical.) */
void popupMessageHandler(const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer data)
{
  const char *criticalText = "** Warning: ";
  const char *errorText = "** Error: ";

  const int len = max(strlen(criticalText), strlen(errorText)) + 1;
  char prefixText[len];
  prefixText[len - 1] = '\0';
  
  if (log_level & G_LOG_LEVEL_ERROR)
    strcpy(prefixText, errorText);
  else if (log_level & G_LOG_LEVEL_CRITICAL)
    strcpy(prefixText, criticalText);
  else
    prefixText[0] = '\0';
  
  /* print to console (indented by the same amount as debug output, if relevant) */
#ifdef DEBUG
  DEBUG_OUT("%s%s", prefixText, (char*)message);
#else
  if (log_level & G_LOG_LEVEL_INFO)
    fprintf(stderr, "%s%s", prefixText, (char*)message); /* info messages to to stderr */
  else
    printf("%s%s", prefixText, (char*)message);
#endif

  BlxMessageData *msgData = data ? (BlxMessageData*)data : NULL;
  GtkWindow *parent = msgData ? msgData->parent : NULL;
  GtkStatusbar *statusBar = msgData ? msgData->statusBar : NULL;
  char *titlePrefix = msgData ? msgData->titlePrefix : NULL;
  
  /* also display in the status bar */
  printMessageToStatusbar(message, statusBar);

  /* Record each message in a list */
  GSList **messageList = getMessageList();
  
  BlxMessage *blxMessage = createBlxMessage(message, log_level);
  *messageList = g_slist_append(*messageList, blxMessage);
  
  /* Display as popup or in a scrolled list, whichever method is active */
  if (getUseScrolledMessages())
    {
      displayMessageAsList(*messageList, titlePrefix ? titlePrefix : "", FALSE);
    }
  else
    {
      displayMessageAsPopup(message, log_level, titlePrefix ? titlePrefix : "", parent, statusBar);
    }
    
  /* Exit elegantly if it's a fatal error */
  if (log_level & G_LOG_LEVEL_ERROR)
    exit(EXIT_FAILURE);
}


/* Print the given message to the given statusbar (passed as the user data May be null if the
 * statusbar widget does not exist). */
static void printMessageToStatusbar(const gchar *message, gpointer data)
{
  if (data && GTK_IS_STATUSBAR(data))
    {
      GtkStatusbar *statusBar = GTK_STATUSBAR(data);

      /* Remove any newline chars */
      char *displayText = g_strdup(message);
      char *cp = strchr(displayText, '\n');
      
      while (cp)
        {
          *cp = ' ';
          cp = strchr(displayText, '\n');
        }
          
      guint contextId = gtk_statusbar_get_context_id(GTK_STATUSBAR(statusBar), "defaultMessages");
      gtk_statusbar_pop(GTK_STATUSBAR(statusBar), contextId);
      gtk_statusbar_push(GTK_STATUSBAR(statusBar), contextId, displayText);
        
      g_free(displayText);
    }
}


/* Display a warning/error message as a popup */
static void displayMessageAsPopup(const gchar *message, GLogLevelFlags log_level, const char *titlePrefix, GtkWindow *parent, GtkStatusbar *statusBar)
{
  char *title = blxprintf("%s%s", titlePrefix, "Error");
  
  GtkWidget *dialog = gtk_dialog_new_with_buttons(title, 
						  parent, 
						  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
						  GTK_STOCK_OK,
						  GTK_RESPONSE_ACCEPT,
						  NULL);
  
  g_free(title);
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  
  GtkWidget *vbox = gtk_vbox_new(FALSE, 12);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), vbox, TRUE, TRUE, 0);
  
  GtkWidget *hbox = gtk_hbox_new(FALSE, 12);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
  
  GtkWidget *image = gtk_image_new_from_stock(getDialogIcon(log_level), GTK_ICON_SIZE_DIALOG);
  gtk_box_pack_start(GTK_BOX(hbox), image, TRUE, TRUE, 0);
  
  GtkWidget *label = gtk_label_new(message);
  gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
  
  GtkWidget *button = gtk_check_button_new_with_mnemonic("Switch to _scrolled message window");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), getUseScrolledMessages());
  widgetSetCallbackData(button, onSetUseScrolledMessages, NULL);
  gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
  
  g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), NULL);
  
  gtk_widget_show_all(vbox);
  
  /* Block until the user clears the dialog */
  gtk_dialog_run(GTK_DIALOG(dialog));
}


/* Display a warning/error message in the scrolled message list */
static void displayMessageAsList(GSList *messageList, const char *titlePrefix, const gboolean bringToFront)
{
  static GtkWidget *dialog = NULL;
  static GtkWidget *button = NULL;
  static GtkTextView *textView = NULL;
  static GtkTextBuffer *textBuffer = NULL;
  static GtkTextTag *normalTag = NULL;       /* for normal text */
  static GtkTextTag *highlightTag = NULL;    /* for highlighting text */
  
  if (!dialog)
    {
      char *title = blxprintf("%s%s", titlePrefix, "Errors");
      
      dialog = gtk_dialog_new_with_buttons(title, 
                                           NULL, 
                                           GTK_DIALOG_DESTROY_WITH_PARENT,
                                           GTK_STOCK_OK,
                                           GTK_RESPONSE_ACCEPT,
                                           NULL);
      g_free(title);
      
      /* Hide the widget instead of destroying it when it is closed */
      g_signal_connect(G_OBJECT(dialog), "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
      g_signal_connect(dialog, "response", G_CALLBACK(onResponseDialog), GINT_TO_POINTER(TRUE));
      
      GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
      gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), vbox, TRUE, TRUE, 0);
      
      /* Create the text view */
      int width = 600;
      int height = 180;
      
      GtkWidget *child = createScrollableTextView(NULL, FALSE, NULL, TRUE, &height, &textView);
      gtk_box_pack_start(GTK_BOX(vbox), child, TRUE, TRUE, 0);
      gtk_window_set_default_size(GTK_WINDOW(dialog), width, height);
      textBuffer = gtk_text_view_get_buffer(textView);
      
      /* Always show the horizontal scrollbar, otherwise it can cover the last line in the list
       * when it magically appears */
      gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(child), GTK_POLICY_ALWAYS, GTK_POLICY_AUTOMATIC);
      
      /* Add some tags for normal/highlighted text */
      normalTag = gtk_text_buffer_create_tag(textBuffer, NULL, "foreground", BLX_BLACK, NULL);
      highlightTag = gtk_text_buffer_create_tag(textBuffer, NULL, "foreground", BLX_RED, NULL);
      
      /* Create a button to allow user to switch back to popup messages */
      button = gtk_check_button_new_with_mnemonic("Switch to _popup messages");
      widgetSetCallbackData(button, onSetUsePopupMessages, NULL);
      gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
    }
  else
    {
      /* Delete contents of text buffer */
      GtkTextIter startIter;
      GtkTextIter endIter;
      gtk_text_buffer_get_start_iter(textBuffer, &startIter);
      gtk_text_buffer_get_end_iter(textBuffer, &endIter);
      
      gtk_text_buffer_delete(textBuffer, &startIter, &endIter);
    }
  
  gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_ACCEPT);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), !getUseScrolledMessages());
  
  addMessagesToBuffer(messageList, textBuffer, normalTag, highlightTag);
  
  /* Don't bring the dialog to the front if it was already open, unless specifically requested. */
  gtk_widget_show_all(dialog);
  
  if (bringToFront || !GTK_WIDGET_VISIBLE(dialog))
    {
      gtk_window_present(GTK_WINDOW(dialog));
    }
  
  /* Scroll to the last line in the text view. We must make sure the text view height has
   * been calculated first or gtk_text_view_scroll_to_iter might not do the correct thing. */
  while(gtk_events_pending()) 
    {
      gtk_main_iteration();
    }
  
  GtkTextIter textIter;
  gtk_text_buffer_get_end_iter(textBuffer, &textIter);
  gtk_text_view_scroll_to_iter(textView, &textIter, 0.0, FALSE, 0.0, 0.0);
}


/* Create a BlxMessage. Result should be destroyed with destroyBlxMessage. Timestamps the
 * message with the current system time. */
static BlxMessage* createBlxMessage(const char *text, const GLogLevelFlags logLevel)
{
  BlxMessage *blxMessage = g_malloc(sizeof(BlxMessage));
  
  blxMessage->text = g_strdup(text);
  blxMessage->timestamp = time(NULL);
  blxMessage->logLevel = logLevel;
  
  return blxMessage;
}

/* free the memory used by a blxMessage */
static void destroyBlxMessage(BlxMessage **blxMessage)
{
  if (blxMessage && *blxMessage)
    {
      if ((*blxMessage)->text)
	{
	  g_free((*blxMessage)->text);
	}
      
      g_free(*blxMessage);
    }
}


/* Get the message list. Creates the list the first time it is requested. */
static GSList** getMessageList()
{
  static GSList *messageList = NULL;
  
  return &messageList;
}


/* Clear the message list, freeing any memory that it uses */
void destroyMessageList()
{
  GSList **messageList = getMessageList();
  
  if (messageList && *messageList)
    {
      GSList *msgItem = *messageList;
      for ( ; msgItem; msgItem = msgItem->next)
	{
	  BlxMessage *blxMessage = (BlxMessage*)(msgItem->data);
	  destroyBlxMessage(&blxMessage);
	}
      
      g_slist_free(*messageList);
      *messageList = NULL;
    }
}


/* Get the display text for a message, including the timestamp if 'incTimestamp' is true.
 * The result should be free'd with g_free */
static char* blxMessageGetDisplayText(const BlxMessage *msg, const gboolean incTimestamp)
{
  /* Use ctime to format the timestamp but cut off the newline char at the end. */
  char *timeText = g_strdup(ctime(&msg->timestamp));
  char *cutPoint = strchr(timeText, '\n');
  if (cutPoint)
    {
      *cutPoint = '\0';
    }
  
  char separatorText[] = " - "; /* separates timestamp from message body */
  
  char *result = g_malloc(strlen(msg->text) + strlen(separatorText) + strlen(timeText) + 1);
  
  sprintf(result, "%s%s%s", timeText, separatorText, msg->text);
  
  return result;
}


/* Add each BlxMessage in the given list to the given text buffer */
static void addMessagesToBuffer(GSList *messageList, 
                                GtkTextBuffer *textBuffer, 
                                GtkTextTag *normalTag, 
                                GtkTextTag *highlightTag)
{
  GSList *msgItem = messageList;  
  
  for ( ; msgItem; msgItem = msgItem->next)
    {
      const BlxMessage *msg = (const BlxMessage*)(msgItem->data);
      
      /* If it's the last in the list, highlight it. */
      GtkTextTag *tag = (msgItem->next == NULL ? highlightTag : normalTag);
      
      char *displayText = blxMessageGetDisplayText(msg, TRUE);
      
      GtkTextIter endIter;
      gtk_text_buffer_get_end_iter(textBuffer, &endIter);
      gtk_text_buffer_insert_with_tags(textBuffer, &endIter, displayText, -1, tag, NULL);
      
      g_free(displayText);
    }
}


/* Internal function to get a pointer to the useScrolledMessages flag. Creates the flag the first 
 * time it is requested.  */
static gboolean* getUseScrolledMessagesPtr()
{
  static gboolean useScrolledMessages = FALSE;
  return &useScrolledMessages;
}


/* Access function to set the 'use scrolled messages' flag */
static gboolean getUseScrolledMessages()
{
  return *getUseScrolledMessagesPtr();
}

/* Access function to set the 'use scrolled messages' flag */
static void setUseScrolledMessages(const gboolean newValue)
{
  gboolean *useScrolledMessages = getUseScrolledMessagesPtr();
  *useScrolledMessages = newValue;
}


/* Called when user requests to view messages as a list */
static gboolean onSetUseScrolledMessages(GtkWidget *button, const gint responseId, gpointer data)
{
  const gboolean useScrolledMessages = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  setUseScrolledMessages(useScrolledMessages);
  return TRUE;
}


/* Called when user requests to view messages as popups */
static gboolean onSetUsePopupMessages(GtkWidget *button, const gint responseId, gpointer data)
{
  const gboolean usePopupMessages = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
  setUseScrolledMessages(!usePopupMessages);
  return TRUE;
}


/* Utility to return the stock icon to use for a dialog for a given message severity */
static char* getDialogIcon(GLogLevelFlags log_level)
{
  char *result = NULL;
  
  if (log_level & G_LOG_LEVEL_ERROR)
        result = GTK_STOCK_DIALOG_ERROR;
  else if (log_level & G_LOG_LEVEL_CRITICAL || log_level & G_LOG_LEVEL_WARNING)
        result = GTK_STOCK_DIALOG_WARNING;
  else
        result = GTK_STOCK_DIALOG_INFO;
    
  return result;
}


/***********************************************************
 *		         Printing			   * 
 ***********************************************************/

/* Utility to get the drawable to print for the given widget */
static GdkDrawable* getPrintDrawable(GtkWidget *widget, const gboolean cachedOnly)
{
  GdkDrawable *drawable = widget->window; /* draw everything by default */
  
  if (cachedOnly)
    {
      if (widgetGetDrawable(widget) && (GTK_IS_LAYOUT(widget) || GTK_IS_DRAWING_AREA(widget)))
        {
          /* We've been asked to print a drawing area, so print the cached drawable directly */
          drawable = widgetGetDrawable(widget);
        }
      else if (GTK_IS_CONTAINER(widget))
        {
          /* Create a blank pixmap to draw the child widgets on to (clear any previous print first) */
          widgetClearCachedDrawable(widget, NULL);
          drawable = createBlankPixmap(widget);
          gtk_container_foreach(GTK_CONTAINER(widget), collatePixmaps, widget);
        }
    }
  
  return drawable;
}


/* Get the scale required to fit the image to the required dimension */
static double getPrintScale(GtkPrintContext *context, GdkDrawable *drawable, PrintScaleType scaleType)
{
  /* Get the page size */
  double ctxWidth = gtk_print_context_get_width(context);
  double ctxHeight = gtk_print_context_get_height(context);
  
  /* Get the image size */
  int imgWidth, imgHeight;
  gdk_drawable_get_size(drawable, &imgWidth, &imgHeight);
  
  double scale = 1;
  
  switch (scaleType)
  {
    case PRINT_FIT_WIDTH:
      scale = ctxWidth / (double)imgWidth;
      break;
      
    case PRINT_FIT_HEIGHT:
      scale = ctxHeight / (double)imgHeight;
      break;
      
    default:
      scale = min(ctxWidth / (double)imgWidth, ctxHeight / (double)imgHeight);
      break;
  };
  
  scale = min(scale, 1.0); /* don't scale larger */
  
  return scale;
}


/* Called after the user clicks ok in the print dialog. Calcualtes the number
 * of pages to print. */
void onBeginPrint(GtkPrintOperation *print, GtkPrintContext *context, gpointer data)
{
  GtkWidget *widget = GTK_WIDGET(data);
  
  GdkDrawable *drawable = getPrintDrawable(widget, g_printCachedOnly);
  double scale = getPrintScale(context, drawable, g_printScaleType);

  /* Get the scaled image size */
  int imgWidth, imgHeight;
  gdk_drawable_get_size(drawable, &imgWidth, &imgHeight);
  imgWidth = (int)((double)imgWidth * scale);
  imgHeight = (int)((double)imgHeight * scale);

  /* Get the page size */
  double ctxWidth = gtk_print_context_get_width(context);
  double ctxHeight = gtk_print_context_get_height(context);

  /* Calculate the number of pages. We always scale to fit in at least one direction
   * so we should only have multiple pages in one dimension. */
  int numPages = 1;
  
  switch (g_printScaleType)
  {
    case PRINT_FIT_WIDTH:
      numPages = ceil((double)imgHeight / ctxHeight);
      break;
      
    case PRINT_FIT_HEIGHT:
      numPages = ceil((double)imgWidth / ctxWidth);
      break;
      
    default:
      numPages = 1;
      break;
  };
  
  gtk_print_operation_set_n_pages(print, numPages);
}


/* This function is called recursively to combine all of the drawables for 
 * any child widgets onto a single drawable for their top-level parent widget,
 * which is passed as the user data. */
void collatePixmaps(GtkWidget *widget, gpointer data)
{
  GtkWidget *parent = GTK_WIDGET(data);
  GdkDrawable *drawable = widgetGetDrawable(widget);
  
  /* If this widget is visible and has a drawable set, draw it onto the window's drawable */
  if (drawable && GTK_WIDGET_VISIBLE(widget))
    {
      /* Get the left edge of the widget in terms of the blixem window's coordinates */
      int xSrc = 0, ySrc = 0;
      int xDest, yDest;
      gtk_widget_translate_coordinates(widget, parent, xSrc, ySrc, &xDest, &yDest);
      
      GdkGC *gc = gdk_gc_new(widget->window);
      gdk_draw_drawable(widgetGetDrawable(parent), gc, drawable, xSrc, ySrc, xDest, yDest, -1, -1); /* -1 means full width/height */
      g_object_unref(gc);
    }
  
  /* If this widget is a container, recurse over its children */
  if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), collatePixmaps, parent);
    }
}


/* Print handler - renders a specific page */
void onDrawPage(GtkPrintOperation *print, GtkPrintContext *context, gint pageNum, gpointer data)
{
  GtkWidget *widget = GTK_WIDGET(data);
  cairo_t *cr = gtk_print_context_get_cairo_context(context);

  GdkDrawable *drawable = getPrintDrawable(widget, g_printCachedOnly);

  double scale = getPrintScale(context, drawable, g_printScaleType);
  
  double ctxWidth = gtk_print_context_get_width(context);
  double ctxHeight = gtk_print_context_get_height(context);
  ctxWidth /= scale;
  ctxHeight /= scale;
  
  /* Get the coords on the image at which this page starts */
  int x = 0;
  int y = 0;

  switch (g_printScaleType)
  {
    case PRINT_FIT_WIDTH:
      y = pageNum * ctxHeight;
      break;
      
    case PRINT_FIT_HEIGHT:
      x = pageNum * ctxWidth;
      break;
      
    default:
      break;
  };
  
  /* Create a new pixmap for drawing just the section for this page */
  GdkDrawable *pagePixmap = gdk_pixmap_new(widget->window, ctxWidth, ctxHeight, -1);
  GdkGC *gc = gdk_gc_new(pagePixmap);
  GdkColor bgColor;
  gdk_color_parse("#ffffff", &bgColor);
  gboolean failures[1];
  gdk_colormap_alloc_colors(gdk_colormap_get_system(), &bgColor, 1, TRUE, TRUE, failures);
  gdk_gc_set_foreground(gc, &bgColor);
  gdk_draw_rectangle(pagePixmap, gc, TRUE, 0, 0, ctxWidth, ctxHeight);
  
  gdk_draw_drawable(pagePixmap, gc, drawable, x, y, 0, 0, ctxWidth, ctxHeight);

  /* Scale the image */
  cairo_scale(cr, scale, scale); 
  
  /* Paint the image */
  gdk_cairo_set_source_pixmap(cr, pagePixmap, 0, 0);
  cairo_paint(cr);
}


/* This function initiates a print operation. 
 *
 * The 'printCachedOnly' argument is an attempt to de-clutter the print
 * by only printing stuff we request by setting a cached drawable on it, i.e. 
 * so that we can exclude scrollbars etc. It works pretty well, except for the case 
 * where you want to include transient items that are not part of the cached
 * drawable in the print (such as the highlight box in blixem or the crosshair
 * in dotter). It's up to the application to set its cached drawables and make
 * sure they include everything necessary if they want to use this option, 
 * otherwise we just draw everything. */
void blxPrintWidget(GtkWidget *widget, 
                    GtkWidget *window,
                    GtkPrintSettings **printSettings, 
                    GtkPageSetup **pageSetup, 
                    const gboolean printCachedOnly,
                    const PrintScaleType scaleType)
{
  g_printCachedOnly = printCachedOnly;
  g_printScaleType = scaleType;
  
  /* Create a print operation, using the same settings as the last print, if there was one */
  GtkPrintOperation *print = gtk_print_operation_new();
  
  if (*printSettings != NULL)
    gtk_print_operation_set_print_settings(print, *printSettings);
  
  if (*pageSetup)
    gtk_print_operation_set_default_page_setup(print, *pageSetup);
  
  
  g_signal_connect (print, "begin_print", G_CALLBACK (onBeginPrint), widget);
  g_signal_connect(G_OBJECT(print), "draw-page", G_CALLBACK(onDrawPage), widget);
  
  /* Pop up the print dialog */
  GtkPrintOperationResult printResult = gtk_print_operation_run (print, 
								 GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
								 GTK_WINDOW(window),
								 NULL);
  
  /* If the user hit ok, remember the print settings for next time */
  if (printResult == GTK_PRINT_OPERATION_RESULT_APPLY)
    {
      if (*printSettings != NULL)
	{
	  g_object_unref(*printSettings);
	  *printSettings = NULL;
	}
      
      *printSettings = g_object_ref(gtk_print_operation_get_print_settings(print));
    }
  
  g_object_unref(print);
}


/* Recursively set the background color for the given widget and all its children */
void setWidgetBackgroundColor(GtkWidget *widget, gpointer data)
{
  GdkColor *color = (GdkColor*)data;
  GtkWidget *window = gtk_widget_get_toplevel(widget);
  
  window->style->bg[GTK_STATE_NORMAL] = *color;
  gtk_widget_modify_bg(widget, GTK_STATE_NORMAL, color);
  gtk_widget_modify_base(widget, GTK_STATE_NORMAL, color);
  
  if (GTK_IS_CONTAINER(widget))
    {
      gtk_container_foreach(GTK_CONTAINER(widget), setWidgetBackgroundColor, data);
    }
}


/* Find an executable and return its complete pathname.
 */
gboolean findCommand (char *command, char **resultOut)
{
  gboolean found = FALSE;
  
#if !defined(NO_POPEN)
  static char result[1025];
  char fileName[1025];
  FILE *file = NULL;
  
  char *pathEnv = g_malloc(strlen(getenv("PATH"))+1);
  strcpy(pathEnv, getenv("PATH"));
  
  /* Try each path in the enviroment var */
  char *path = strtok(pathEnv, ":");
  
  while (path) 
    {
      strcpy(fileName, path);
      strcat(fileName,"/");
      strcat(fileName, command);
      
      /* Check if the file exists in this path */
      file = fopen(fileName, "r");
      if (file)  //!access(fileName, F_OK) && !access(fileName, X_OK)) 
	{
	  found = TRUE;
	  fclose(file);
	  break;
	}
      
      path = strtok(0, ":");
    }
  
  g_free(pathEnv);
  
  if (found) 
    {
      strcpy(result, fileName);
      found = TRUE;
    }
  else 
    {
      char *msg = blxprintf("Can't find executable '%s' in path.\n", command);
      strcpy(result, msg);
      g_free(msg);
    }
  
  if (resultOut)
    {
      *resultOut = result;
    }
#else
  char *msg = blxprintf("Can't open executable '%s' - popen command is not defined.\n", command);
  strcpy(result, msg);
  g_free(msg);
#endif
  
  return found;
}


/* Emits the size-allocate signal on the given widget */
void forceResize(GtkWidget *widget)
{
  g_signal_emit_by_name(G_OBJECT(widget), "size-allocate", &widget->allocation);
}


/***********************************************************
 *		         Combo boxes			   * 
 ***********************************************************/

/* Callback for when the value in a 2-column combo box has changed. 
 * The value to update is an enum, a pointer to which is passed in the user
 * data.
 * This is called as a result of a response on a dialog, and the response
 * function of the dialog is responsible for updating the display. */
gboolean onComboChanged(GtkWidget *combo, const gint responseId, gpointer data)
{
  int *result = (int*)data;
  GtkTreeIter iter;
  
  if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(combo), &iter))
    {
      GtkTreeModel *model = gtk_combo_box_get_model(GTK_COMBO_BOX(combo));
      
      GValue val = {0};
      gtk_tree_model_get_value(model, &iter, COMBO_ENUM_COL, &val);
      
      *result = g_value_get_int(&val);
    }
  
  return TRUE;
}


/* Utility to add an item to our 2-column combo box */
void addComboItem(GtkComboBox *combo,
                         GtkTreeIter *parent, 
                         const int val,
                         const char *text,
                         const int initVal)
{
  GtkTreeStore *store = GTK_TREE_STORE(gtk_combo_box_get_model(combo));
  
  GtkTreeIter iter;
  gtk_tree_store_append(store, &iter, parent);
  
  gtk_tree_store_set(store, &iter, COMBO_ENUM_COL, val, COMBO_TEXT_COL, text, -1);
  
  if (val == initVal)
    {
      gtk_combo_box_set_active_iter(combo, &iter);
    }
}


/* Utility to create a 2-column combo box with column1 as an enum and column2 
 * as a text description */
GtkComboBox* createComboBox()
{
  /* Create a data-store for the combo box, and the combo itself */
  GtkTreeStore *store = gtk_tree_store_new(N_COMBO_COLUMNS, G_TYPE_INT, G_TYPE_STRING);
  GtkComboBox *combo = GTK_COMBO_BOX(gtk_combo_box_new_with_model(GTK_TREE_MODEL(store)));
  g_object_unref(store);
  
  /* Create a cell renderer to display the text column. */
  GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
  gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, FALSE);
  gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", COMBO_TEXT_COL, NULL);
  
  return combo;
}


/***********************************************************
 *		         Files    			   * 
 ***********************************************************/

/* Utility to ask the user for a file to save to. Returns the file name (or
 * NULL if the user cancels). The current file name, default file path, 
 * default file extension and dialog title can all be specified, or passed
 * as null to use defaults. */
const char* getSaveFileName(GtkWidget *widget, 
                            const char *currentName,   
                            const char *defaultPath,
                            const char *defaultExtension,
                            const char *title)
{
  const char *filename = NULL;
  
  GtkWindow *window = widget ? GTK_WINDOW(gtk_widget_get_toplevel(widget)) : NULL;
  
  GtkWidget *dialog = gtk_file_chooser_dialog_new (title ? title : "Save File",
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
  else if (defaultExtension)
    {
      gtk_file_chooser_set_current_name (GTK_FILE_CHOOSER (dialog), defaultExtension);
    }
  
  
  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
    {
      filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    }
  
  gtk_widget_destroy (dialog);
  return filename;
}

/* Utility to ask the user for a file to load. Returns the file name (or
 * NULL if the user cancels). The default file path and dialog title can be
 * specified, or passed as null to use defaults. */
const char* getLoadFileName(GtkWidget *widget, 
                            const char *defaultPath,
                            const char *title)
{
  const char *filename = NULL;
  
  GtkWindow *window = widget ? GTK_WINDOW(gtk_widget_get_toplevel(widget)) : NULL;
  
  GtkWidget *dialog = gtk_file_chooser_dialog_new (title ? title : "Open File",
                                                   window,
                                                   GTK_FILE_CHOOSER_ACTION_OPEN,
                                                   GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                                   GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                                   NULL);
  
  if (defaultPath)
    {
      gtk_file_chooser_set_current_folder (GTK_FILE_CHOOSER (dialog), defaultPath);
    }
  
  if (gtk_dialog_run (GTK_DIALOG (dialog)) == GTK_RESPONSE_ACCEPT)
    {
      filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (dialog));
    }
  
  gtk_widget_destroy (dialog);
  return filename;
}


/***********************************************************
 *		         Action groups			   * 
 ***********************************************************/

/* Utility to enable/disable an item in a menu. The action name must be the value of a valid action. */
void enableMenuAction(GtkActionGroup *action_group, const char *actionName, const gboolean enable)
{
  GtkAction *action = gtk_action_group_get_action(action_group, actionName);
  
  if (!action)
    g_warning("Error %s menu item: action '%s' not found.\n", (enable ? "enabling" : "disabling"), actionName);
  else
    gtk_action_set_sensitive(action, enable);
}


/* Utility to set the status of a toggle item in a menu. The action name must be the name 
 * of a valid toggle action. */
void setToggleMenuStatus(GtkActionGroup *action_group, const char *actionName, const gboolean active)
{
  GtkAction *action = gtk_action_group_get_action(action_group, actionName);
  
  if (!action)
    {
      g_warning("Error toggling menu item: action '%s' not found.\n", actionName);
    }
  else if (!GTK_IS_TOGGLE_ACTION(action))
    {
      g_warning("Error toggling menu item: action '%s' is not a valid toggle action.\n", actionName);
    }
  else
    {
      gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(action), active);
    }
}


/***********************************************************
 *		         External calls			   * 
 ***********************************************************/

/* call an external shell command and print output in a text_scroll window
 *
 * This is a replacement for the old graph based text window, it has the advantage
 * that it uses gtk directly and provides cut/paste/scrolling but...it has the
 * disadvantage that it will use more memory as it collects all the output into
 * one string and then this is _copied_ into the text widget.
 * 
 * If this proves to be a problem I expect there is a way to feed the text to the
 * text widget a line a time. */
GtkWidget* externalCommand (char *command, char *progName, GtkWidget *widget, GError **error)
{
  GtkWidget *resultWindow = NULL;

#if !defined(MACINTOSH)
  
  GString *resultText = getExternalCommandOutput(command, error);
  
  if (!error || *error == NULL)
    {
      char *title = blxprintf("%s - %s", progName, command);
      resultWindow = displayFetchResults(title, resultText->str, widget, NULL);
      
      g_free(title);
      g_string_free(resultText, TRUE);
    }
  
#endif
  
  return resultWindow;
}


/* Execute the given external command and return the output from the command as a GString. 
 * The result should be free'd with g_string_free. */
GString* getExternalCommandOutput(const char *command, GError **error)
{
  GString *resultText = g_string_new(NULL) ;
  
  char lineText[MAXLINE+1];
  
  FILE *pipe = popen (command, "r") ;
  
  if (pipe && !feof(pipe) && fgets (lineText, MAXLINE, pipe))
    {
      while (!feof (pipe))
        { 
          int len = strlen(lineText);
          
          if (len > 0)
            { 
              if (lineText[len-1] == '\n') 
                {
                  lineText[len-1] = '\0';
                }
              
              g_string_append_printf(resultText, "%s\n", lineText) ;
            }
          
          if (!fgets (lineText, MAXLINE, pipe))
            {
              break;
            }
        }
    }
  else
    {
      g_set_error(error, SEQTOOLS_ERROR, SEQTOOLS_ERROR_EXECUTING_CMD, "Error executing command: %s\n", command);
    }

  pclose (pipe);

  return resultText;
}


/* Display a message dialog showing the given display text. This utility functions sets
 * things like the font and default width based on properties of the main blixem window.
 * Returns a pointer to the dialog, and optionally sets a pointer to the text buffer in the
 * textBuffer return argument. */
GtkWidget* displayFetchResults(const char *title, 
                               const char *displayText, 
                               GtkWidget *widget, 
                               GtkTextBuffer **textBuffer)
{
  /* Use a fixed-width font */
  const char *fontFamily = findFixedWidthFont(widget);
  PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
  
  /* Set the initial width based on the default number of characters wide */
//  PangoContext *context = gtk_widget_get_pango_context(widget);
//  PangoFontMetrics *metrics = pango_context_get_metrics(context, fontDesc, pango_context_get_language(context));
//  getFontCharSize(widget, fontDesc, &charWidth, NULL);
  int charWidth = 8;

  const int initWidth = DEFAULT_PFETCH_WINDOW_WIDTH_CHARS * charWidth;
  const int maxHeight = widget->allocation.height * DEFAULT_PFETCH_WINDOW_HEIGHT_FRACTION;
  
  GtkTextView *textView = NULL;
  GtkWidget *result = showMessageDialog(title, displayText, NULL, initWidth, maxHeight, FALSE, FALSE, fontDesc, &textView);
  
  if (textBuffer && textView)
    {
      *textBuffer = gtk_text_view_get_buffer(textView);
    }
  
  /* Clean up */
//  pango_font_metrics_unref(metrics);
//  pango_font_description_free(fontDesc);
  
  return result;
}

