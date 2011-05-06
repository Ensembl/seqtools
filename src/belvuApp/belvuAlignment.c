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

#include <belvuApp/belvuAlignment.h>
#include <math.h>


#define DEFAULT_XPAD                            2
#define DEFAULT_YPAD                            2
#define DEFAULT_PADDING_CHARS                   1  /* number of char widths to use to pad between columns */
#define DEFAULT_NAME_COLUMN_PADDING_CHARS       2  /* number of char widths to use to pad after the name column */
#define WRAP_DISPLAY_PADDING_CHARS              4

/* Properties specific to the belvu alignment */
typedef struct _BelvuAlignmentProperties
{
  BelvuContext *bc;	            /* The belvu context */
  
  GtkWidget *headersArea;           /* Drawing widget for the header columns (i.e. name and coord columns) */
  GtkWidget *seqArea;               /* Drawing widget for the actual sequence */
  GdkRectangle headersRect;         /* Drawing area for the header columns (i.e. name and coord columns) */
  GdkRectangle seqRect;             /* Drawing area for the actual sequence */
  GtkAdjustment *hAdjustment;       /* Controls horizontal scroll bar */
  
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

static void drawText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const char *text)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, text);
  gdk_draw_layout(drawable, gc, x, y, layout);
  g_object_unref(layout);
}

static void drawIntAsText(GtkWidget *widget, GdkDrawable *drawable, GdkGC *gc, const int x, const int y, const int value)
{
  char *tmpStr = blxprintf("%d", value);
  drawText(widget, drawable, gc, x, y, tmpStr);
  g_free(tmpStr);
}


/* Convert one of the old acedb-style color numbers to a hex string */
static const char* convertColorNumToStr(const int colorNum)
{
  const char *result = NULL;
  
  switch (colorNum)
    {
      case WHITE: result = BLX_WHITE; break;
      case BLACK: result = BLX_BLACK; break;
      case LIGHTGRAY: result = BLX_LIGHT_GREY; break;
      case DARKGRAY: result = BLX_DARK_GREY; break;
      case RED: result = BLX_RED; break;
      case GREEN: result = BLX_GREEN; break;
      case BLUE: result = BLX_BLUE; break;
      case YELLOW: result = BLX_YELLOW; break;
      case CYAN: result = BLX_CYAN; break;
      case MAGENTA: result = BLX_MAGENTA; break;
      case LIGHTRED: result = BLX_LIGHT_RED; break;
      case LIGHTGREEN: result = BLX_LIGHT_GREEN; break;
      case LIGHTBLUE: result = BLX_SKY_BLUE; break;
      case DARKRED: result = BLX_DARK_RED; break;
      case DARKGREEN: result = BLX_DARK_GREEN; break;
      case DARKBLUE: result = BLX_DARK_BLUE; break;
      case PALERED: result = BLX_LIGHT_RED; break;
      case PALEGREEN: result = BLX_LIGHT_GREEN; break;
      case PALEBLUE: result = BLX_PALE_BLUE; break;
      case PALEYELLOW: result = BLX_PALE_YELLOW; break;
      case PALECYAN: result = BLX_LIGHT_CYAN; break;
      case PALEMAGENTA: result = BLX_PALE_MAGENTA; break;
      case BROWN: result = BLX_BROWN; break;
      case ORANGE: result = BLX_ORANGE; break;
      case PALEORANGE: result = BLX_PALE_ORANGE; break;
      case PURPLE: result = BLX_PURPLE; break;
      case VIOLET: result = BLX_VIOLET; break;
      case PALEVIOLET: result = BLX_PALE_VIOLET; break;
      case GRAY: result = BLX_GREY; break;
      case PALEGRAY: result = BLX_VERY_LIGHT_GREY; break;
      case CERISE: result = BLX_CERISE; break;
      case MIDBLUE: result = BLX_MID_BLUE; break;
      default: result = BLX_WHITE; break;
    };
  
  return result;
}


/* Convert an old-style ACEDB color number to a GdkColor */
static void convertColorNumToGdkColor(const int colorNum, GdkColor *result)
{
  const char *colorStr = convertColorNumToStr(colorNum);
  getColorFromString(colorStr, result, NULL);
}
  
  
/* Find the background color of the given residue index i in the given alignment */
static void findResidueBGcolor(BelvuContext *bc, ALN* alnp, int i, GdkColor *result) 
{
  int colorNum = 0;

  if (bc->lowercaseOn && alnp->seq[i] >= 'a' && alnp->seq[i] <= 'z') 
    colorNum = ORANGE;
  else if (alnp->nocolor)
    colorNum = BOXCOLOR;
  else if (alnp->markup)
    colorNum = getMarkupColor(alnp->seq[i]);
  else if (bc->color_by_conserv || bc->colorByResIdOn)
    colorNum = getConservColor(bc, alnp->seq[i], i);
  else
    colorNum = getColor(alnp->seq[i]);

  convertColorNumToGdkColor(colorNum, result);
  
  /* If this alignment is selected, use a slightly different shade to highlight it */
  if (alnp == bc->highlightedAln)
    getSelectionColor(result, result);
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
  
  const gboolean isSelected = (alnp == properties->bc->highlightedAln);
  GdkGC *gc = gdk_gc_new(drawable);
  
  if (isSelected)
    {
      /* Highlight the background */
      GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, properties->bc->defaultColors, isSelected, FALSE);
      gdk_gc_set_foreground(gc, bgColor);
      gdk_draw_rectangle(drawable, gc, TRUE, properties->headersRect.x, y, properties->headersRect.width, properties->charHeight);
    }
  
  GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, isSelected, FALSE);
  gdk_gc_set_foreground(gc, textColor);

  /* Draw the name */
  if (alnp->name)
    {
      x += properties->columnPadding;
      drawText(widget, drawable, gc, x, y, alnp->name);
      x += (properties->bc->maxNameLen * properties->charWidth) + properties->columnPadding;
    }
  
  /* Draw the coords */
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  x += properties->columnPadding;
  drawIntAsText(widget, drawable, gc, x, y, alnp->start);
  x += (properties->bc->maxStartLen * properties->charWidth) + properties->columnPadding;
  
  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);
  x += properties->columnPadding;
  drawIntAsText(widget, drawable, gc, x, y, alnp->end);
  x += (properties->bc->maxEndLen * properties->charWidth) + properties->columnPadding;

  gdk_draw_line(drawable, gc, x, y, x, y + properties->charHeight);

  g_object_unref(gc);
}


/* Draw a single line in the sequence area */
static void drawSingleSequence(GtkWidget *widget,
                                GdkDrawable *drawable, 
                                BelvuAlignmentProperties *properties,
                                ALN *alnp, 
                                const int lineNum)
{
  if (!alnp->seq)
    return;
  
  const int startX = 0;
  int x = startX;
  const int y = properties->seqRect.y + (properties->charHeight * lineNum);
  
  GdkGC *gc = gdk_gc_new(drawable);
  GtkAdjustment *hAdjustment = properties->hAdjustment;
  
  /* Loop through each character in the current display range and color
   * the background */
  int i = hAdjustment->value;
  const int displayLen = hAdjustment->page_size;
  const int iMax = hAdjustment->value + displayLen;
  
  for ( ; i < iMax; ++i)
    {
      GdkColor bgColor;
      findResidueBGcolor(properties->bc, alnp, i, &bgColor);
      gdk_gc_set_foreground(gc, &bgColor);
      
      gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
      x += properties->charWidth;
    }

  /* Draw the text for the current display range */
  GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, textColor);
  
  char *cp = alnp->seq + (int)hAdjustment->value;
  char *displayText = g_strndup(cp, displayLen);
  drawText(widget, drawable, gc, startX, y, displayText);
  
  g_free(displayText);
  g_object_unref(gc);
}


static gboolean colorsEqual(GdkColor *color1, GdkColor *color2)
{
  return (color1->pixel = color2->pixel);
}


/* Return the foreground color selected for a given
   background color in bc->color_by_conserv mode 
 */
static void bg2fgColor(BelvuContext *bc, GdkColor *bgColor, GdkColor *result)
{
  convertColorNumToGdkColor(bc->maxbgColor, result);
  if (colorsEqual(bgColor, result))
    return convertColorNumToGdkColor(bc->maxfgColor, result);
  
  convertColorNumToGdkColor(bc->midbgColor, result);
  if (colorsEqual(bgColor, result))
    return convertColorNumToGdkColor(bc->midfgColor, result);
  
  convertColorNumToGdkColor(bc->lowbgColor, result);
  if (colorsEqual(bgColor, result))
    return convertColorNumToGdkColor(bc->lowfgColor, result);
  
  /* Anything else is either uncolored or markup - make them black */
  convertColorNumToGdkColor(BLACK, result);
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
      
      drawText(widget, drawable, gcText, x, y, properties->title);
      line += 2; /* increment, and also add another blank line for spacing */
    }
  
  while (paragraph * properties->wrapWidth + totCollapsed < bc->maxLen)
    {
      for (j = 0; j < bc->alignArr->len; ++j)
	{
          ALN *alnp = &g_array_index(bc->alignArr, ALN, j);
          
          if (alnp->hide) 
            continue;
          
          int alnstart = paragraph*properties->wrapWidth +totCollapsed; 
          int alnlen = ( (paragraph+1)*properties->wrapWidth +totCollapsed < bc->maxLen ? properties->wrapWidth : bc->maxLen - alnstart );
          int alnend = alnstart + alnlen;
          
          gboolean empty = TRUE;
          for (empty=1, i = alnstart; i < alnend; i++) 
            {
              if (!isGap(alnp->seq[i]) && alnp->seq[i] != ' ') 
                {
                  empty = FALSE;
                  break;
                }
            }
          
          if (!empty) 
            {
              const int y = line * properties->charHeight;

              for (collapsePos = 0, oldpos = pos[j], i = alnstart; i < alnend; i++) 
                {	
                  const int xpos = bc->maxNameLen + bc->maxEndLen + numSpaces + i - alnstart - collapsePos;
                  const int x = xpos * properties->charWidth;
                  
                  if (alnp->seq[i] == '[') 
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
                  
                  if (alnp->seq[i] == ']') 
                    {
                      collapseOn = 0;
                      collapsePos -= 4;
                      alnend -= 3;
                      
                      sprintf(collapseStr, "[%3d]", collapseRes);
                      drawText(widget, drawable, gcText, x, y, collapseStr);

                      continue;
                    }
                  
                  if (collapseOn) 
                    {
                      collapsePos++;
                      if (!isGap(alnp->seq[i])) 
                        {
                          collapseRes++;
                          pos[j]++;
                        }
                      alnend++;
                      continue;
                    }
                  
                  
                  if (!isGap(alnp->seq[i]) && alnp->seq[i] != ' ') 
                    {
                      findResidueBGcolor(bc, alnp, i, pBgColor);
                      gdk_gc_set_foreground(gc, pBgColor);
                      
                      gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
                      
                      pos[j]++;
                    }
                  
                  /* Foreground color */
                  GdkColor fgColor;
                  if (bc->color_by_conserv)
                    bg2fgColor(bc, pBgColor, &fgColor);
                  else
                    convertColorNumToGdkColor(BLACK, &fgColor);
                  
                  gdk_gc_set_foreground(gc, &fgColor);

                  *ch = alnp->seq[i];
                  drawText(widget, drawable, gc, x, y, ch);
                }

              drawText(widget, drawable, gcText, properties->charWidth, y, alnp->name);
              
              if (!alnp->markup) 
                {
                  char *tmpStr = blxprintf("%*d", bc->maxEndLen, oldpos);
                  drawText(widget, drawable, gcText, (bc->maxNameLen + 3) * properties->charWidth, y, tmpStr);
                  g_free(tmpStr);
                  
                  if (alnend == bc->maxLen) 
                    {
                      char *tmpStr = blxprintf("%-d", pos[j] - 1);
                      drawText(widget, drawable, gcText, (bc->maxNameLen + bc->maxEndLen + alnlen + 5) * properties->charWidth, y, tmpStr);
                      g_free(tmpStr);
                    }
                }
              
              line++;
            }
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
  
  /* Loop through each alignment */
  int i = 0;
  for ( ; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      drawSingleHeader(widget, drawable, properties, alnp, i);
    }
}


/* Draw the alignment view */
static void drawBelvuSequence(GtkWidget *widget, GdkDrawable *drawable, BelvuAlignmentProperties *properties)
{
  BelvuContext *bc = properties->bc;
  
  /* Loop through each alignment */
  int i = 0;
  for ( ; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      drawSingleSequence(widget, drawable, properties, alnp, i);
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

/* Called when the horizontal scroll range has changed */
static void onHScrollRangeChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  /* Set the horizontal adjustment page size to be the number of characters
   * that will fit in the width */
  properties->hAdjustment->page_size = (int)(properties->seqArea->allocation.width / properties->charWidth);
  properties->hAdjustment->page_increment = properties->hAdjustment->page_size;
  
  widgetClearCachedDrawable(properties->seqArea, NULL);
  gtk_widget_queue_draw(properties->seqArea);
}

/* Called when the vertical scroll range upper value has changed (i.e. the length
 * of the alignments array has changed). */
void updateOnVScrollSizeChaged(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  GtkAdjustment *vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(properties->seqArea));
  vAdjustment->upper = properties->bc->alignArr->len;
  
  belvuAlignmentRedrawAll(belvuAlignment);
  gtk_adjustment_changed(vAdjustment);
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
        properties->columnPadding + (properties->bc->maxEndLen * properties->charWidth);
    }
  
  return result;
}


static int getAlignmentDisplayHeight(BelvuAlignmentProperties *properties)
{
  /* Calculate the height required to draw one row for each sequence */
  int result = properties->bc->alignArr->len * properties->charHeight + (2 * DEFAULT_YPAD);
  
  if (properties->wrapWidth != UNSET_INT)
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
  properties->seqRect.y = DEFAULT_YPAD;
  properties->seqRect.width = getAlignmentDisplayWidth(properties);
  properties->seqRect.height = getAlignmentDisplayHeight(properties);
  
  if (properties->headersArea)
    {
      /* There is a separate drawing area for the columns that contain the row
       * headers; calculate its size. */
      properties->headersRect.x = DEFAULT_XPAD;
      properties->headersRect.y = DEFAULT_YPAD;
      properties->headersRect.height = properties->seqRect.height;
      
      properties->headersRect.width = (properties->bc->maxNameLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxStartLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->bc->maxEndLen * properties->charWidth) + (2 * properties->columnPadding) +
                                      (properties->columnPadding / 2);
      
      
      gtk_layout_set_size(GTK_LAYOUT(properties->headersArea), properties->headersRect.width, properties->headersRect.height);
  
      if (properties->headersRect.width != properties->headersArea->allocation.width)
        gtk_widget_set_size_request(properties->headersArea, properties->headersRect.width, -1);
    }
  
  gtk_layout_set_size(GTK_LAYOUT(properties->seqArea), properties->seqRect.width, properties->seqRect.height);

  if (properties->hAdjustment)
    gtk_adjustment_changed(properties->hAdjustment); /* signal that the scroll range has changed */
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
  if (!bc->highlightedAln) 
    {
      g_critical("Please select a sequence to remove.\n");
      return;
    }
  
  const int idx = bc->highlightedAln->nr - 1;
  g_array_remove_index(bc->alignArr, idx);
  arrayOrder(bc->alignArr);
  
  bc->saved = FALSE;
  
  g_message("Removed %s/%d-%d.  %d sequences left.  Esc to cancel.\n\n", bc->highlightedAln->name, bc->highlightedAln->start, bc->highlightedAln->end, bc->alignArr->len);
  
  bc->highlightedAln = NULL;
  
  rmFinaliseGapRemoval(bc);
  
  updateOnVScrollSizeChaged(belvuAlignment);  
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


/* Center the vertical scrollbar on the currently highlighted sequence */
void centerHighlighted(BelvuContext *bc, GtkWidget *belvuAlignment)
{
  if (bc->highlightedAln) 
    {
      BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
      GtkAdjustment *vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(properties->seqArea));
      
      vAdjustment->value = bc->highlightedAln->nr * properties->charHeight;
      gtk_adjustment_changed(vAdjustment);
    }
}


/* Mouse button handler */
static gboolean onButtonPressBelvuAlignment(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;
  
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);

  if (event->type == GDK_BUTTON_PRESS && event->button == 1) /* left click */
    {
      /* Select the clicked sequence */
      const int idx = (event->y - properties->seqRect.y) / properties->charHeight;
    
      properties->bc->highlightedAln = &g_array_index(properties->bc->alignArr, ALN, idx);
    
      belvuAlignmentRedrawAll(belvuAlignment);
    
      handled = TRUE;
    }
  else if (event->type == GDK_2BUTTON_PRESS && event->button == 1) /* double-click */
    {
      if (properties->bc->removingSeqs)
	{
	  /* Removed the clicked sequence (which will be the selected one) */
	  removeSelectedSequence(properties->bc, belvuAlignment);
	}
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
  /* We'll put everything in a table */
  GtkWidget *belvuAlignment = gtk_table_new(2, 2, FALSE);
  const int xpad = 0, ypad = 0;
  
  /* Create the sequence drawing area and the columns drawing area (for the
   * name and coords etc). The sequence drawing area will have a horizontal 
   * scrollbar, and both sections will share the same vertical scrollbar. If the
   * display is wrapped, there is no headers area, so just use the default 
   * layout scrollbars for the sequence area.*/
  GtkAdjustment *hAdjustment = NULL;
  
  if (wrapWidth == UNSET_INT)
    {
      hAdjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0, 0, bc->maxLen + 1, 1, 0, 0));
      GtkWidget *hScrollbar = gtk_hscrollbar_new(hAdjustment);
      gtk_table_attach(GTK_TABLE(belvuAlignment), hScrollbar, 2, 3, 2, 3, GTK_EXPAND | GTK_FILL, GTK_SHRINK, xpad, ypad);

      g_signal_connect(G_OBJECT(hAdjustment), "value-changed", G_CALLBACK(onHScrollPosChangedBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(hAdjustment), "changed", G_CALLBACK(onHScrollRangeChangedBelvuAlignment), belvuAlignment);
    }
  
  /* Create the sequence area */
  GtkWidget *seqArea = gtk_layout_new(NULL, NULL);
  GtkWidget *seqScrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_container_add(GTK_CONTAINER(seqScrollWin), seqArea);
  gtk_table_attach(GTK_TABLE(belvuAlignment), seqScrollWin, 2, 3, 1, 2, GTK_EXPAND | GTK_FILL, GTK_EXPAND | GTK_FILL, xpad, ypad);
  
  g_signal_connect(G_OBJECT(seqArea), "expose-event", G_CALLBACK(onExposeBelvuSequence), belvuAlignment);  
  g_signal_connect(G_OBJECT(seqArea), "size-allocate", G_CALLBACK(onSizeAllocateBelvuAlignment), belvuAlignment);

  /* If we've already got a horizontal scrollbar, just set the scrollwin to use the vertical one. */
  if (hAdjustment)
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(seqScrollWin), GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
  else
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(seqScrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);

  /* Only create the header columns if the display is not wrapped */
  GtkWidget *headersArea = NULL;

  if (wrapWidth == UNSET_INT)
    {
      GtkAdjustment *vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(seqArea));
      headersArea = gtk_layout_new(NULL, vAdjustment);
      gtk_table_attach(GTK_TABLE(belvuAlignment), headersArea, 1, 2, 1, 2, GTK_FILL, GTK_EXPAND | GTK_FILL, xpad, ypad);
      g_signal_connect(G_OBJECT(headersArea), "expose-event", G_CALLBACK(onExposeBelvuColumns), belvuAlignment);  
    
      /* Also only connect button press handler if this is a standard (i.e. not wrapped) alignment window */
      gtk_widget_add_events(seqArea, GDK_BUTTON_PRESS_MASK);
      gtk_widget_add_events(headersArea, GDK_BUTTON_PRESS_MASK);
      g_signal_connect(G_OBJECT(seqArea), "button-press-event", G_CALLBACK(onButtonPressBelvuAlignment), belvuAlignment);
      g_signal_connect(G_OBJECT(headersArea), "button-press-event", G_CALLBACK(onButtonPressBelvuAlignment), belvuAlignment);
    }
  
  /* Set the style and properties */
  setBelvuAlignmentStyle(bc, seqArea, headersArea);
  belvuAlignmentCreateProperties(belvuAlignment, bc, headersArea, seqArea, hAdjustment, title, wrapWidth);

  return belvuAlignment;
}





