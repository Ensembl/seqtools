/*  File: belvuAlignment.c
 *  Author: Gemma Barson, 2011-04-12
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
 * Description: see belvuAlignment.h
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvuAlignment.h"
#include "belvuApp/belvuWindow.h"
#include <math.h>


#define DEFAULT_XPAD                            2
#define DEFAULT_YPAD                            2
#define DEFAULT_PADDING_CHARS                   1  /* number of char widths to use to pad between columns */
#define DEFAULT_NAME_COLUMN_PADDING_CHARS       2  /* number of char widths to use to pad after the name column */
#define WRAP_DISPLAY_PADDING_CHARS              4


/* Local function declarations */
static void               bg2fgColor(BelvuContext *bc, GdkColor *bgColor, GdkColor *result);
static gboolean           highlightAlignment(BelvuContext *bc, ALN *alnp);


/* Properties specific to the belvu alignment */
typedef struct _BelvuAlignmentProperties
{
  BelvuContext *bc;	            /* The belvu context */
  
  GtkWidget *headersArea;           /* Drawing widget for the header columns (i.e. name and coord columns) */
  GtkWidget *seqArea;               /* Drawing widget for the actual sequence */
  GdkRectangle headersRect;         /* Drawing area for the header columns (i.e. name and coord columns) */
  GdkRectangle seqRect;             /* Drawing area for the actual sequence */
  GtkAdjustment *hAdjustment;       /* Controls horizontal scroll bar */
  GtkAdjustment *vAdjustment;       /* Controls vertical scroll bar */
  
  int columnPadding;                /* Padding in between columns, in pixels */
  int nameColumnPadding;            /* Padding after name column, in pixels */
  char *title;                /* title to display at the top of the alignments */
  int wrapWidth;                    /* Number of characters after which to wrap (or UNSET_INT for no wrapping) */
  
  PangoFontDescription *fontDesc;   /* The fixed-width font to use for displaying the alignment */
  gdouble charWidth;                /* The width of each character in the display */
  gdouble charHeight;               /* The height of each character in the display */
} BelvuAlignmentProperties;


/***********************************************************
 *                         Properties                      *
 ***********************************************************/

static BelvuAlignmentProperties* belvuAlignmentGetProperties(GtkWidget *widget)
{
  return widget ? (BelvuAlignmentProperties*)(g_object_get_data(G_OBJECT(widget), "BelvuAlignmentProperties")) : NULL;
}

static void onDestroyBelvuAlignment(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  if (properties)
    {
      if (properties->title)
        {
          g_free(properties->title);
          properties->title = NULL;
        }
      
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuAlignment), "BelvuAlignmentProperties", NULL);
    }
}


/* Create the properties struct and initialise all values. */
static void belvuAlignmentCreateProperties(GtkWidget *belvuAlignment, 
                                           BelvuContext *bc,
                                           GtkWidget *headersArea,
                                           GtkWidget *seqArea,
                                           GtkAdjustment *hAdjustment,
                                           GtkAdjustment *vAdjustment,
                                           const char *title,
                                           const int wrapWidth)
{
  if (belvuAlignment)
    {
      BelvuAlignmentProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->headersArea = headersArea;
      properties->seqArea = seqArea;
      properties->columnPadding = 0; /* calculated in calculate-borders */
      properties->hAdjustment = hAdjustment;
      properties->vAdjustment = vAdjustment;
      properties->title = g_strdup(title);
      properties->wrapWidth = wrapWidth;
      
      /* Find a fixed-width font */
      const char *fontFamily = findFixedWidthFont(belvuAlignment);
      PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
      pango_font_description_set_size(fontDesc, pango_font_description_get_size(belvuAlignment->style->font_desc));
      gtk_widget_modify_font(seqArea, fontDesc);

      if (headersArea)
        gtk_widget_modify_font(headersArea, fontDesc);

      getFontCharSize(belvuAlignment, fontDesc, &properties->charWidth, &properties->charHeight);
      properties->fontDesc = fontDesc;

      g_object_set_data(G_OBJECT(belvuAlignment), "BelvuAlignmentProperties", properties);
      g_signal_connect(G_OBJECT(belvuAlignment), "destroy", G_CALLBACK (onDestroyBelvuAlignment), NULL);
    }
}


/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Clear cached drawables and redraw all */
void belvuAlignmentRedrawAll(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  widgetClearCachedDrawable(properties->seqArea, NULL);
  
  if (properties->headersArea)
    widgetClearCachedDrawable(properties->headersArea, NULL);
  
  gtk_widget_queue_draw(belvuAlignment);
}


/* Find the background color of the given residue index i in the given alignment */
static void findResidueBGcolor(BelvuContext *bc, ALN* alnp, int i, const gboolean isSelected, GdkColor *result) 
{
  int colorNum = 0;
  char *alnpSeq = alnGetSeq(alnp);
  
  if (!alnpSeq || i >= alnGetSeqLen(alnp))
    colorNum = WHITE;
  else if (bc->lowercaseOn && alnpSeq[i] >= 'a' && alnpSeq[i] <= 'z') 
    colorNum = ORANGE;
  else if (alnp->nocolor)
    colorNum = BOXCOLOR;
  else if (alnp->markup)
    colorNum = getMarkupColor(alnpSeq[i]);
  else if (colorByConservation(bc) || colorByResId(bc))
    colorNum = getConservColor(bc, alnpSeq[i], i);
  else
    colorNum = getColor(alnpSeq[i]);

  convertColorNumToGdkColor(colorNum, isSelected, result);
}


/* Draw a single line in the headers area */
static void drawSingleHeader(GtkWidget *widget,
                             GdkDrawable *drawable, 
                             BelvuAlignmentProperties *properties,
                             ALN *alnp, 
                             const int lineNum)
{
  int x = 0;
  const int y = properties->headersRect.y + (properties->charHeight * lineNum);
  
  const gboolean highlightAln = highlightAlignment(properties->bc, alnp);
  GdkGC *gc = gdk_gc_new(drawable);
  
  if (highlightAln)
    {
      /* Highlight the background */
      GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, properties->bc->defaultColors, highlightAln, FALSE);
      gdk_gc_set_foreground(gc, bgColor);
      gdk_draw_rectangle(drawable, gc, TRUE, properties->headersRect.x, y, properties->headersRect.width, properties->charHeight);
    }
  
  GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, highlightAln, FALSE);
  gdk_gc_set_foreground(gc, textColor);

  /* Draw the name */
  if (alnp->name)
    {
      x += properties->columnPadding;
      drawText(widget, drawable, gc, x, y, alnp->name, NULL, NULL);
      x += (properties->bc->maxNameLen * properties->charWidth) + properties->columnPadding;
    }
  
  /* Draw the coords (but only if this is not GC markup; note that we always
   * draw the column separator lines for the coords, though) */
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  x += properties->columnPadding;

  if (alnp->markup != GC)
    drawIntAsText(widget, drawable, gc, x, y, alnp->start, NULL, NULL);
  
  x += (properties->bc->maxStartLen * properties->charWidth) + properties->columnPadding;
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  x += properties->columnPadding;

  if (alnp->markup != GC)
    drawIntAsText(widget, drawable, gc, x, y, alnp->end, NULL, NULL);

  x += (properties->bc->maxEndLen * properties->charWidth) + properties->columnPadding;
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  
  /* Draw the score, if displaying scores (and if not a markup row) */
  if (properties->bc->displayScores)
    {
      x += properties->columnPadding;
      
      if (!alnp->markup)
        {
          char *tmpStr = blxprintf("%*.1f", alnp->score);
          drawText(widget, drawable, gc, x, y, tmpStr, NULL, NULL);
          g_free(tmpStr);
        }
      
      x += (properties->bc->maxScoreLen * properties->charWidth) + properties->columnPadding;
      gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
    }
  
  g_object_unref(gc);
}


/* Utility to return true if the given alignment should be highlighted */
static gboolean highlightAlignment(BelvuContext *bc, ALN *alnp)
{
  /* Return true if this alignment has the same name as the selected alignment */
  return (bc->selectedAln && stringsEqual(alnp->name, bc->selectedAln->name, TRUE));
}


/* Draw a single character in the given sequence */
static void drawSequenceChar(BelvuAlignmentProperties *properties, 
                             ALN *alnp,
                             const int colIdx,
                             const gboolean rowHighlighted,
                             GtkWidget *widget,
                             GdkDrawable *drawable,
                             GdkGC *gc,
                             GdkColor *defaultFgColor,
                             const int x, 
                             const int y)
{
  /* This base should be highlighted if the row is highlighted OR if the
   * column is highlighted (but not both) */
  const gboolean colHighlighted = (colIdx == properties->bc->highlightedCol - 1);
  const gboolean highlight = (rowHighlighted != colHighlighted);
  
  GdkColor bgColor;
  
  if (properties->bc->displayColors)
    {
      /* Draw the background */
      findResidueBGcolor(properties->bc, alnp, colIdx, highlight, &bgColor);
      gdk_gc_set_foreground(gc, &bgColor);
      gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
    }
  
  if (colIdx < alnGetSeqLen(alnp))
    {
      /* Draw the text. We get the text colour from the background color in 
       * color-by-conservation mode (if displaying colours is enabled) */
      if (colorByConservation(properties->bc) && properties->bc->displayColors)
        {
          GdkColor fgColor;
          bg2fgColor(properties->bc, &bgColor, &fgColor);
          gdk_gc_set_foreground(gc, &fgColor);
        }
      else
        {
          gdk_gc_set_foreground(gc, defaultFgColor);
        }
      
      char displayText[2];
      displayText[0] = alnGetSeq(alnp)[colIdx];
      displayText[1] = '\0';
      drawText(widget, drawable, gc, x, y, displayText, NULL, NULL);
    }
  
}


/* Draw a single line in the sequence area */
static void drawSingleSequence(GtkWidget *widget,
                                GdkDrawable *drawable, 
                                BelvuAlignmentProperties *properties,
                                ALN *alnp, 
                                const int lineNum)
{
  if (!alnGetSeq(alnp))
    return;
  
  const int startX = 0;
  int x = startX;
  const int y = properties->seqRect.y + (properties->charHeight * lineNum);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GtkAdjustment *hAdjustment = properties->hAdjustment;
  
  GdkColor *defaultFgColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, FALSE, FALSE);
  
  /* Loop through each column in the current display range and color
   * the text and background according to the relevant highlight colors */
  int colIdx = hAdjustment->value;
  int iMax = min(properties->bc->maxLen + 1, hAdjustment->value + hAdjustment->page_size);
  const gboolean rowHighlighted = highlightAlignment(properties->bc, alnp);
  
  for ( ; colIdx < iMax; ++colIdx)
    {
      drawSequenceChar(properties, alnp, colIdx, rowHighlighted, widget, drawable, gc, defaultFgColor, x, y);
      
      /* Increment the x position */
      x += properties->charWidth;
    }
  
  g_object_unref(gc);
}


static gboolean colorsEqual(GdkColor *color1, GdkColor *color2)
{
  return (color1->pixel == color2->pixel &&
          color1->red == color2->red &&
          color1->green == color2->green &&
          color1->blue == color2->blue);
}


/* Return the foreground color selected for a given
   background color in color-by-conservation mode 
 */
static void bg2fgColor(BelvuContext *bc, GdkColor *bgColor, GdkColor *result)
{
  int fgColorNum = BLACK; /* default colour for uncoloured or markup */
  
  /* See if the given bgcolor matches the max-conservation background colour */
  int testColor = *getConsColor(bc, CONS_LEVEL_MAX, FALSE);
  convertColorNumToGdkColor(testColor, FALSE, result);
  
  if (colorsEqual(bgColor, result))
    {
      fgColorNum = *getConsColor(bc, CONS_LEVEL_MAX, TRUE);
    }
  else
    {
      /* Try the mid-conservation colour */
      int testColor = *getConsColor(bc, CONS_LEVEL_MID, FALSE);
      convertColorNumToGdkColor(testColor, FALSE, result);
      
      if (colorsEqual(bgColor, result))
        {
          fgColorNum = *getConsColor(bc, CONS_LEVEL_MID, TRUE);
        }
      else
        {
          /* Lastly, try the low-conservation colour */
          int testColor = *getConsColor(bc, CONS_LEVEL_LOW, FALSE);
          convertColorNumToGdkColor(testColor, FALSE, result);
          
          if (colorsEqual(bgColor, result))
            fgColorNum = *getConsColor(bc, CONS_LEVEL_LOW, TRUE);
        }
    }

  /* Convert the foreground colour number into the GdkColor result. */
  convertColorNumToGdkColor(fgColorNum, FALSE, result);
}


static void drawWrappedSequences(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  /*  Column makeup:
   space, maxNameLen, 2xspace, maxEndLen, space, maxLen, maxEndLen
   i.e. 4 spaces
   */
  
  int i, oldpos, collapseRes = 0, collapsePos = 0, collapseOn = 0;
  int numSpaces = WRAP_DISPLAY_PADDING_CHARS;
  char collapseStr[10],
  ch[2] = " ";
  static int *pos=0;			/* Current residue position of sequence j */
  
  GdkColor bgColor;
  GdkColor *pBgColor = &bgColor;
  GdkGC *gcText = gdk_gc_new(drawable);
  GdkGC *gc = gdk_gc_new(drawable);
  
  BelvuContext *bc = properties->bc;
  
  
  /* Initialise the sequence and position array */
  char *seq = g_malloc(properties->wrapWidth + 1);
  
  if (!pos) 
    pos = (int *)g_malloc(bc->alignArr->len * sizeof(int *));
  
  int j = 0;
  for (j = 0; j < bc->alignArr->len; ++j) 
    pos[j] = g_array_index(bc->alignArr, ALN, j).start;
  
  
  int paragraph = 0;
  int totCollapsed = 0;
  int line = 1; /* leave a blank line (line=0) at the top for spacing */
  
  if (properties->title)
    {
      const int x = properties->columnPadding;
      const int y = line * properties->charHeight;
      
      drawText(widget, drawable, gcText, x, y, properties->title, NULL, NULL);
      line += 2; /* increment, and also add another blank line for spacing */
    }
  
  while (paragraph * properties->wrapWidth + totCollapsed < bc->maxLen)
    {
      for (j = 0; j < bc->alignArr->len; ++j)
	{
          ALN *alnp = &g_array_index(bc->alignArr, ALN, j);
          char *alnpSeq = alnGetSeq(alnp);
          
          if (alnp->hide) 
            continue;
          
          int alnstart = paragraph*properties->wrapWidth +totCollapsed; 
          int alnlen = ( (paragraph+1)*properties->wrapWidth +totCollapsed < bc->maxLen ? properties->wrapWidth : bc->maxLen - alnstart );
          int alnend = min(alnstart + alnlen, alnGetSeqLen(alnp));
          
          gboolean empty = TRUE;
          for (empty=1, i = alnstart; i < alnend; i++) 
            {
              if (alnpSeq && !isGap(alnpSeq[i]) && alnpSeq[i] != ' ') 
                {
                  empty = FALSE;
                  break;
                }
            }
          
          const int y = line * properties->charHeight;

          if (!empty) 
            {
              for (collapsePos = 0, oldpos = pos[j], i = alnstart; i < alnend; i++) 
                {	
                  const int xpos = bc->maxNameLen + bc->maxEndLen + numSpaces + i - alnstart - collapsePos;
                  const int x = xpos * properties->charWidth;
                  
                  if (alnpSeq[i] == '[') 
                    {
                      /* Experimental - collapse block: "[" -> Collapse start, "]" -> Collapse end, e.g.
                       1 S[..........]FKSJFE
                       2 L[RVLKLKYKNS]KDEJHF
                       3 P[RL......DK]FKEJFJ
                                 |
                                 V
                       1 S[  0]FKSJFE
                       2 L[ 10]KDEJHF
                       3 P[  4]FKEJFJ
                       
                       Minimum collapse = 5 chars. The system depends on absolute coherence in format; 
                       very strange results will be generated otherwise. Edges need to be avoided manually.
                       */
                      collapseOn = 1;
                      collapsePos += 1;
                      collapseRes = 0;
                      continue;
                    }
                  
                  if (alnpSeq[i] == ']') 
                    {
                      collapseOn = 0;
                      collapsePos -= 4;
                      alnend -= 3;
                      
                      sprintf(collapseStr, "[%3d]", collapseRes);
                      drawText(widget, drawable, gcText, x, y, collapseStr, NULL, NULL);

                      continue;
                    }
                  
                  if (collapseOn) 
                    {
                      collapsePos++;
                      if (!isGap(alnpSeq[i])) 
                        {
                          collapseRes++;
                          pos[j]++;
                        }
                      alnend++;
                      continue;
                    }
                  
                  
                  if (!isGap(alnpSeq[i]) && alnpSeq[i] != ' ') 
                    {
                      findResidueBGcolor(bc, alnp, i, FALSE, pBgColor);
                      gdk_gc_set_foreground(gc, pBgColor);
                      
                      gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
                      
                      pos[j]++;
                    }
                  
                  /* Foreground color */
                  GdkColor fgColor;
                  if (colorByConservation(bc))
                    bg2fgColor(bc, pBgColor, &fgColor);
                  else
                    convertColorNumToGdkColor(BLACK, FALSE, &fgColor);
                  
                  gdk_gc_set_foreground(gc, &fgColor);

                  *ch = alnpSeq[i];
                  drawText(widget, drawable, gc, x, y, ch, NULL, NULL);
                }
            }
          
          drawText(widget, drawable, gcText, properties->charWidth, y, alnp->name, NULL, NULL);
          
          if (!alnp->markup) 
            {
              char *tmpStr = blxprintf("%*d", bc->maxEndLen, oldpos);
              drawText(widget, drawable, gcText, (bc->maxNameLen + 3) * properties->charWidth, y, tmpStr, NULL, NULL);
              g_free(tmpStr);
              
              if (alnend == bc->maxLen) 
                {
                  char *tmpStr = blxprintf("%-d", pos[j] - 1);
                  drawText(widget, drawable, gcText, (bc->maxNameLen + bc->maxEndLen + alnlen + 5) * properties->charWidth, y, tmpStr, NULL, NULL);
                  g_free(tmpStr);
                }
            }
          
          line++;
        }
      
      paragraph++;
      line++;
      totCollapsed += collapsePos;
    }
  
  g_free(seq);
  g_object_unref(gc);
  g_object_unref(gcText);
}


/* Draw the alignment view */
static void drawBelvuColumns(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;
  
  /* Loop through each visible alignment */
  GtkAdjustment *vAdjustment = properties->vAdjustment;
  
  int i = vAdjustment->value;
  int lineNum = 0;
  const int iMax = min(bc->alignArr->len + 2, vAdjustment->value + vAdjustment->page_size);
  
  for ( ; i < iMax; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (!alnp->hide)
        {
          drawSingleHeader(widget, drawable, properties, alnp, lineNum);
          ++lineNum;
        }
    }
}


/* Draw the alignment view */
static void drawBelvuSequence(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;
  
  /* Loop through each visible alignment */
  GtkAdjustment *vAdjustment = properties->vAdjustment;
  
  int i = vAdjustment->value;
  int lineNum = 0;
  const int iMax = min(bc->alignArr->len + 2, vAdjustment->value + vAdjustment->page_size);
  
  for ( ; i < iMax; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (!alnp->hide)
        {
          drawSingleSequence(widget, drawable, properties, alnp, lineNum);
          ++lineNum;
        }
    }
}


/* Expose handler for the alignment section */
static gboolean onExposeBelvuColumns(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;

  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->headersRect.width, properties->headersRect.height);
          drawBelvuColumns(widget, bitmap, properties);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
        }
      else
	{
	  g_warning("Failed to draw Belvu alignment [%p] - could not create bitmap.\n", widget);
	}
    }

  return TRUE;
}

/* Expose handler for the alignment section */
static gboolean onExposeBelvuSequence(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, properties->seqRect.width, properties->seqRect.height);
          
          if (properties->wrapWidth == UNSET_INT)
            drawBelvuSequence(widget, bitmap, properties);
          else
            drawWrappedSequences(widget, bitmap, properties);
        }
      
      if (bitmap)
        {
          /* Push the bitmap onto the window */
          GdkGC *gc = gdk_gc_new(window);
          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
          g_object_unref(gc);
        }
      else
	{
	  g_warning("Failed to draw Belvu alignment [%p] - could not create bitmap.\n", widget);
	}
    }
  
  return TRUE;
}

/***********************************************************
 *                         Scrolling                       *
 ***********************************************************/

/* Called when the horizontal scroll position has changed */
static void onHScrollPosChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->wrapWidth == UNSET_INT)
    {
      /* We need to clear and re-draw the sequence area (unless we're displaying
       * wrapped sequences, in which case we've already drawn the full sequence
       * width in the drawing area */
      widgetClearCachedDrawable(properties->seqArea, NULL);
      gtk_widget_queue_draw(properties->seqArea);
    }
}

/* Called when the vertical scroll position has changed */
static void onVScrollPosChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->wrapWidth == UNSET_INT)
    {
      belvuAlignmentRedrawAll(belvuAlignment);
    }
}


/* Called when the horizontal scroll range has changed */
static void onHScrollRangeChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  /* Set the horizontal adjustment page size to be the number of characters
   * that will fit in the width */
  properties->hAdjustment->page_size = (int)(properties->seqArea->allocation.width / properties->charWidth);
  properties->hAdjustment->page_increment = properties->hAdjustment->page_size;

  /* Make sure the scroll position still lies within upper - page_size. */
  if (properties->hAdjustment->value > properties->hAdjustment->upper - properties->hAdjustment->page_size)
    properties->hAdjustment->value = properties->hAdjustment->upper - properties->hAdjustment->page_size;

  /* Make sure we're still within the lower limit */
  if (properties->hAdjustment->value < properties->hAdjustment->lower)
    properties->hAdjustment->value = properties->hAdjustment->lower;

  widgetClearCachedDrawable(properties->seqArea, NULL);
  gtk_widget_queue_draw(properties->seqArea);
}

/* Called when the vertical scroll range has changed */
static void onVScrollRangeChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  /* Set the horizontal adjustment page size to be the number of characters
   * that will fit in the width */
  properties->vAdjustment->page_size = (int)(properties->seqArea->allocation.height / properties->charHeight);
  properties->vAdjustment->page_increment = properties->vAdjustment->page_size;
  
  /* Make sure the scroll position still lies within upper - page_size. */
  if (properties->vAdjustment->value > properties->vAdjustment->upper - properties->vAdjustment->page_size)
    properties->vAdjustment->value = properties->vAdjustment->upper - properties->vAdjustment->page_size;
  
  /* Make sure we're still within the lower limit */
  if (properties->vAdjustment->value < properties->vAdjustment->lower)
    properties->vAdjustment->value = properties->vAdjustment->lower;

  belvuAlignmentRedrawAll(belvuAlignment);
}


/* Called when the vertical scroll range upper value has changed (i.e. the length
 * of the alignments array has changed). */
void updateOnVScrollSizeChaged(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  GtkAdjustment *vAdjustment = properties->vAdjustment;
  vAdjustment->upper = properties->bc->alignArr->len + 2;
  
  belvuAlignmentRedrawAll(belvuAlignment);
  gtk_adjustment_changed(vAdjustment);
}


/* This should be called when the alignment length has changed, i.e. when columns
 * have been deleted. */
void updateOnAlignmentLenChanged(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (properties->hAdjustment)
    {
      properties->hAdjustment->upper = properties->bc->maxLen + 1;
      gtk_adjustment_changed(properties->hAdjustment);
    }
}

/***********************************************************
 *                         Sizing                          *
 ***********************************************************/

/* Get the width of the drawing area that shows the sequences */
static int getAlignmentDisplayWidth(BelvuAlignmentProperties *properties)
{
  int result = 0;
  
  if (properties->wrapWidth == UNSET_INT)
    {
      /* Drawing can be very slow if we draw the full alignment, so we only
       * draw what's in the current display width. */
      result = properties->seqArea->allocation.width - properties->columnPadding;
    }
  else
    {
      /* The sequences are wrapped, so we only have a short width to deal with
       * and can therefore draw the full width. */
      result =
        properties->columnPadding + (properties->bc->maxNameLen * properties->charWidth) + 
        properties->nameColumnPadding + (properties->bc->maxEndLen * properties->charWidth) +
        properties->columnPadding + (properties->wrapWidth * properties->charWidth) + 
        properties->columnPadding + (properties->bc->maxEndLen * properties->charWidth) + 
        properties->columnPadding + (properties->bc->maxScoreLen * properties->charWidth);
    }
  
  return result;
}


static int getAlignmentDisplayHeight(BelvuAlignmentProperties *properties)
{
  int result = 0;
  
  if (properties->wrapWidth == UNSET_INT)
    {
      result = properties->seqArea->allocation.height - properties->charHeight - (2 * DEFAULT_YPAD);
    }  
  else
    {
      /* Divide the full sequence length by the display width to give the number
       * of paragraphs we require */
      const int displayLen = properties->wrapWidth;
      const int numParagraphs = ceil((double)properties->bc->maxLen / (double)displayLen);

      /* Multiply the standard height (plus one row for separation) by the number 
       * of paragraphs */
      result = (result + properties->charHeight) * numParagraphs;
      
      /* Add space for the title, if given (one line for the title and one as a spacer) */
      if (properties->title)
        result += properties->charHeight * 2;
    }
  
  return result;
}


static void calculateBelvuAlignmentBorders(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  /* Set the column padding. We use slightly different padding if the view is
   * wrapped because we don't draw column separators */
  if (properties->wrapWidth != UNSET_INT)
    {
      properties->columnPadding = (DEFAULT_PADDING_CHARS * properties->charWidth);
      properties->nameColumnPadding = (DEFAULT_NAME_COLUMN_PADDING_CHARS * properties->charWidth);
    }
  else
    {
      properties->columnPadding = (DEFAULT_PADDING_CHARS * properties->charWidth) / 2;
      properties->nameColumnPadding = 0;
    }
  
  /* Calculate the size of the drawing area for the sequences */
  properties->seqRect.x = DEFAULT_XPAD;
  properties->seqRect.y = DEFAULT_YPAD * 2 + properties->charHeight;
  properties->seqRect.width = getAlignmentDisplayWidth(properties);
  properties->seqRect.height = getAlignmentDisplayHeight(properties);
  
  if (properties->headersArea)
    {
      /* There is a separate drawing area for the columns that contain the row
       * headers; calculate its size. */
      properties->headersRect.x = DEFAULT_XPAD;
      properties->headersRect.y = DEFAULT_YPAD * 2 + properties->charHeight;
      properties->headersRect.height = properties->seqRect.height;
      
      properties->headersRect.width = (properties->bc->maxNameLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxStartLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxEndLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxScoreLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->columnPadding / 2);
      
      
      gtk_layout_set_size(GTK_LAYOUT(properties->headersArea), properties->headersRect.width, properties->headersRect.height);
  
      if (properties->headersRect.width != properties->headersArea->allocation.width)
        gtk_widget_set_size_request(properties->headersArea, properties->headersRect.width, -1);
    }
  
  gtk_layout_set_size(GTK_LAYOUT(properties->seqArea), properties->seqRect.width, properties->seqRect.height);

  if (properties->hAdjustment)
    gtk_adjustment_changed(properties->hAdjustment); /* signal that the scroll range has changed */

  if (properties->vAdjustment)
    gtk_adjustment_changed(properties->vAdjustment);
}


static void onSizeAllocateBelvuAlignment(GtkWidget *widget, GtkAllocation *allocation, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  calculateBelvuAlignmentBorders(belvuAlignment);
}


/***********************************************************
 *                Removing sequences                       *
 ***********************************************************/

void removeSelectedSequence(BelvuContext *bc, GtkWidget *belvuAlignment)
{ 
  if (!bc->selectedAln) 
    {
      g_critical("Please select a sequence to remove.\n");
      return;
    }
  
  const int idx = bc->selectedAln->nr - 1;
  g_array_remove_index(bc->alignArr, idx);
  arrayOrder(bc->alignArr);
  
  bc->saved = FALSE;
  
  g_message("Removed %s/%d-%d.  %d sequences left.\n\n", bc->selectedAln->name, bc->selectedAln->start, bc->selectedAln->end, bc->alignArr->len);
  
  bc->selectedAln = NULL;
  
  rmFinaliseGapRemoval(bc);
  updateOnVScrollSizeChaged(belvuAlignment);  
  onRowSelectionChanged(bc);
}


/* Get rid of seqs that are too gappy. */
void removeGappySeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  rmGappySeqs(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of partial seqs. */
void removePartialSeqs(BelvuContext *bc, GtkWidget *belvuAlignment)
{
  rmPartialSeqs(bc);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of redundant seqs (those that are more than the given % identical). */
void removeRedundantSeqs(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  mkNonRedundant(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of outlier seqs (those that are less than the given % identical to any other). */
void removeOutliers(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  rmOutliers(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);  
}

/* Get rid of seqs that have a score below the given value */
void removeByScore(BelvuContext *bc, GtkWidget *belvuAlignment, const double cutoff)
{
  rmScore(bc, cutoff);
  updateOnVScrollSizeChaged(belvuAlignment);
  centerHighlighted(bc, belvuAlignment);
}


/* Scroll the alignment list if necessary to make sure that the currently
 * highlighted sequence is visible */
void centerHighlighted(BelvuContext *bc, GtkWidget *belvuAlignment)
{
  if (bc->selectedAln) 
    {
      BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
      GtkAdjustment *vAdjustment = properties->vAdjustment;
      
      /* Create the range of y coords that we want to be in range */
      const int newIdx = bc->selectedAln->nr - 1;
      const int lower = properties->seqRect.y + (newIdx * properties->charHeight);
      const int upper = lower + properties->charHeight;
      
      gtk_adjustment_clamp_page(vAdjustment, lower, upper);
    }
}


/* Select the row at the given y coord */
static void selectRowAtCoord(BelvuAlignmentProperties *properties, const int y)
{
  BelvuContext *bc = properties->bc;
  
  const int rowIdx = (y - properties->seqRect.y) / properties->charHeight;
  
  /* Set the selected alignment. Note that some alignments are hidden so
   * we need to count how many visible alignments there are up to this row. */
  int i = 0;
  int count = -1;
  bc->selectedAln = NULL;
  
  for ( i = 0; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (!alnp->hide)
        {
          ++count;
          if (count == rowIdx)
            {
              bc->selectedAln = alnp;
              break;
            }
        }
    }
  
  /* Highlight any alignments that have the same name as the selected alignment */
  g_slist_free(bc->highlightedAlns); /* clear current list */
  bc->highlightedAlns = NULL;
  
  for (i = 0; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      
      if (highlightAlignment(bc, alnp))
        bc->highlightedAlns = g_slist_prepend(bc->highlightedAlns, alnp);
    }
}


/* Utility to find the column index (0-based from left edge) at the given
 * x coordinate in the sequence area of the alignment view. */
static int getColumnAtCoord(BelvuAlignmentProperties *properties, const int x)
{
  int result = (x - properties->seqRect.x) / properties->charWidth;
  return result;
}


/* Select the column at the given x coord */
static void selectColumnAtCoord(BelvuAlignmentProperties *properties, const int x, const gboolean highlightCol)
{
  BelvuContext *bc = properties->bc;

  const int colIdx = getColumnAtCoord(properties, x);
  
  if (colIdx >= 0 && colIdx < bc->maxLen)
    {
      bc->selectedCol = colIdx + properties->hAdjustment->value + 1;
      
      if (highlightCol)
        bc->highlightedCol = bc->selectedCol;
      else
        bc->highlightedCol = 0;
    }
}


/* Mouse button handler for the headers area of the alignment window */
static gboolean onButtonPressHeadersArea(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->type == GDK_BUTTON_PRESS && event->button == 1)  /* single click left button */
    {
      /* Select the clicked row */
      selectRowAtCoord(properties, event->y);
      onRowSelectionChanged(properties->bc);
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double click left button */
    {
      if (properties->bc->removingSeqs)
	{
	  /* Removed the clicked sequence (which will be the selected one) */
	  removeSelectedSequence(properties->bc, belvuAlignment);
	}
    }
  
  return handled;
}


/* Mouse button handler for the sequence area of the alignment window */
static gboolean onButtonPressSeqArea(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->type == GDK_BUTTON_PRESS &&
      (event->button == 1 || event->button == 2))  /* single click left or middle buttons */
    {
      /* If the middle button was pressed, highlight the selected column */
      const gboolean highlightCol = (event->button == 2);
      
      /* If the left button was pressed, select the clicked row */
      if (event->button == 1)
        {
          selectRowAtCoord(properties, event->y);
          onRowSelectionChanged(properties->bc);
        }
      
      /* Select the clicked column */
      selectColumnAtCoord(properties, event->x, highlightCol);
      onColSelectionChanged(properties->bc);
      
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double click left button */
    {
      if (properties->bc->removingSeqs)
	{
	  /* Removed the clicked sequence (which will be the selected one) */
	  removeSelectedSequence(properties->bc, belvuAlignment);
	}
    }
  
  return handled;
}


/* Mouse button handler release for the sequence area of the alignment window */
static gboolean onButtonReleaseSeqArea(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->button == 2)  /* released middle button */
    {
      /* Scroll to centre on the current column */
      int colIdx = getColumnAtCoord(properties, event->x);
      double newValue = colIdx + properties->hAdjustment->value - (properties->hAdjustment->page_size / 2.0);
      
      if (newValue > properties->hAdjustment->upper - properties->hAdjustment->page_size)
        newValue = properties->hAdjustment->upper - properties->hAdjustment->page_size;
        
      gtk_adjustment_set_value(properties->hAdjustment, newValue);
      
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double click left button */
    {
      if (properties->bc->removingSeqs)
	{
	  /* Removed the clicked sequence (which will be the selected one) */
	  removeSelectedSequence(properties->bc, belvuAlignment);
	}
    }
  
  return handled;
}


/* Mouse move handler for the sequence area of the alignment window */
static gboolean onMouseMoveSeqArea(GtkWidget *widget, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  if (event->state & GDK_BUTTON2_MASK) /* middle button pressed */
    {
      /* Update the selected/highlighted column */
      selectColumnAtCoord(properties, event->x, TRUE);
      onColSelectionChanged(properties->bc);
      
      handled = TRUE;
    }
  
  return handled;
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

/* Set style properties for the belvu alignment widgets */
static void setBelvuAlignmentStyle(BelvuContext *bc, GtkWidget *seqArea, GtkWidget *headersArea)
{
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);

  gtk_widget_modify_bg(seqArea, GTK_STATE_NORMAL, bgColor);

  if (headersArea)
    gtk_widget_modify_bg(headersArea, GTK_STATE_NORMAL, bgColor);
}


/* Create a widget that will draw the alignments. This type of widget is
 * used for both the standard and wrapped views - pass wrapWidth as UNSET_INT
 * for the standard view, or pass wrapWidth as the number of characters after
 * which to wrap for the wrapped view. */
GtkWidget* createBelvuAlignment(BelvuContext *bc, const char* title, const int wrapWidth)
{
  GtkWidget *belvuAlignment = NULL;      /* the outermost container */
  GtkWidget *seqArea = NULL;             /* the "sequence area", which will draw the actual peptide sequences */
  GtkWidget *headersArea = NULL;         /* optional "headers" column for drawing names, coords etc. */
  GtkAdjustment *hAdjustment = NULL;
  GtkAdjustment *vAdjustment = NULL;
  
  if (wrapWidth == UNSET_INT)
    {
      /* The standard (unwrapped) view. Two columns that share the same vertical
       * scrollbar but only the sequence area will have a horizontal scrollbar. */
      
      /* Place everything in a table */
      belvuAlignment = gtk_table_new(2, 3, FALSE);
      const int xpad = 0, ypad = 0;

      /* Create the scrollbars */
      hAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, 0, bc->maxLen + 1, 1, 0, 0));
      vAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, 0, bc->alignArr->len + 2, 1, 0, 0));
      GtkWidget *hScrollbar = gtk_hscrollbar_new(hAdjustment);
      GtkWidget *vScrollbar = gtk_vscrollbar_new(vAdjustment);
      
      /* Create the drawing areas */
      seqArea = gtk_layout_new(NULL, NULL);
      headersArea = gtk_layout_new(NULL, NULL);
      
      /* Place all the widgets in the table */
      gtk_table_attach(GTK_TABLE(belvuAlignment), headersArea,  0, 1, 0, 1, GTK_FILL, GTK_EXPAND | GTK_FILL, xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), seqArea,      1, 2, 0, 1, GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), hScrollbar,   1, 2, 1, 2, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);
      gtk_table_attach(GTK_TABLE(belvuAlignment), vScrollbar,   2, 3, 0, 1, GTK_SHRINK, GTK_EXPAND | GTK_FILL, xpad, ypad);

      /* Connect signals */
      gtk_widget_add_events(headersArea, GDK_BUTTON_PRESS_MASK);
      gtk_widget_add_events(seqArea, GDK_BUTTON_PRESS_MASK);
      gtk_widget_add_events(seqArea, GDK_BUTTON_RELEASE_MASK);
      gtk_widget_add_events(seqArea, GDK_POINTER_MOTION_MASK);

      g_signal_connect(G_OBJECT(headersArea), "expose-event", G_CALLBACK(onExposeBelvuColumns), belvuAlignment);  
      g_signal_connect(G_OBJECT(seqArea), "expose-event", G_CALLBACK(onExposeBelvuSequence), belvuAlignment);  
      g_signal_connect(G_OBJECT(seqArea), "size-allocate", G_CALLBACK(onSizeAllocateBelvuAlignment), belvuAlignment);
      
      g_signal_connect(G_OBJECT(headersArea), "button-press-event", G_CALLBACK(onButtonPressHeadersArea), belvuAlignment);
      g_signal_connect(G_OBJECT(seqArea), "button-press-event", G_CALLBACK(onButtonPressSeqArea), belvuAlignment);
      g_signal_connect(G_OBJECT(seqArea), "button-release-event", G_CALLBACK(onButtonReleaseSeqArea), belvuAlignment);
      g_signal_connect(G_OBJECT(seqArea), "motion-notify-event", G_CALLBACK(onMouseMoveSeqArea), belvuAlignment);

      g_signal_connect(G_OBJECT(hAdjustment), "changed", G_CALLBACK(onHScrollRangeChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(vAdjustment), "changed", G_CALLBACK(onVScrollRangeChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(hAdjustment), "value-changed", G_CALLBACK(onHScrollPosChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(vAdjustment), "value-changed", G_CALLBACK(onVScrollPosChangedBelvuAlignment), belvuAlignment);
    }
  else
    {
      /* For the wrapped view, we just have a single layout area. We'll need to
       * draw the entire contents at once so that it can be printed. We can
       * therefore just place it in a standard scrolled window. We also don't
       * need to worry about button-press events because this widget is not
       * interactive; it's just for printing. */
      seqArea = gtk_layout_new(hAdjustment, vAdjustment);

      belvuAlignment = gtk_scrolled_window_new(NULL, NULL);
      gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(belvuAlignment), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
      gtk_container_add(GTK_CONTAINER(belvuAlignment), seqArea);
      
      hAdjustment = gtk_layout_get_hadjustment(GTK_LAYOUT(seqArea));
      vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(seqArea));
      
      g_signal_connect(G_OBJECT(seqArea), "expose-event", G_CALLBACK(onExposeBelvuSequence), belvuAlignment);  
      g_signal_connect(G_OBJECT(seqArea), "size-allocate", G_CALLBACK(onSizeAllocateBelvuAlignment), belvuAlignment);
    }
  
  /* Set the style and properties */
  g_assert(belvuAlignment);
  setBelvuAlignmentStyle(bc, seqArea, headersArea);
  belvuAlignmentCreateProperties(belvuAlignment, bc, headersArea, seqArea, hAdjustment, vAdjustment, title, wrapWidth);

  return belvuAlignment;
}





