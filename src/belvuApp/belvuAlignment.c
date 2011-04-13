/*  File: belvuAlignment.c
 *  Author: Gemma Barson, 2011-04-12
 *  Copyright (c) 2009 - 2010 Genome Research Ltd
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


#define DEFAULT_XPAD                2
#define DEFAULT_YPAD                2
#define DEFAULT_PADDING_CHARS       1  /* number of char widths to use to pad between columns */


/* Properties specific to the belvu window */
typedef struct _BelvuAlignmentProperties
{
  BelvuContext *bc;	            /* The belvu context */
  
  GtkWidget *headersArea;           /* Drawing widget for the header columns (i.e. name and coord columns) */
  GtkWidget *seqArea;               /* Drawing widget for the actual sequence */
  GdkRectangle headersRect;         /* Drawing area for the header columns (i.e. name and coord columns) */
  GdkRectangle seqRect;             /* Drawing area for the actual sequence */
  
  int columnPadding;                /* Padding in between columns, in pixels */
  
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
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(belvuAlignment), "BelvuAlignmentProperties", NULL);
    }

  gtk_main_quit();  
}


/* Create the properties struct and initialise all values. */
static void belvuAlignmentCreateProperties(GtkWidget *belvuAlignment, 
                                           BelvuContext *bc,
                                           GtkWidget *headersArea,
                                           GtkWidget *seqArea)
{
  if (belvuAlignment)
    {
      BelvuAlignmentProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->headersArea = headersArea;
      properties->seqArea = seqArea;
      properties->columnPadding = 0; /* calculated in calculate-borders */
      
      /* Find a fixed-width font */
      const char *fontFamily = findFixedWidthFont(belvuAlignment);
      PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
      pango_font_description_set_size(fontDesc, pango_font_description_get_size(belvuAlignment->style->font_desc));
      gtk_widget_modify_font(headersArea, fontDesc);
      gtk_widget_modify_font(seqArea, fontDesc);
      
      getFontCharSize(belvuAlignment, fontDesc, &properties->charWidth, &properties->charHeight);
      properties->fontDesc = fontDesc;

      g_object_set_data(G_OBJECT(belvuAlignment), "BelvuAlignmentProperties", properties);
      g_signal_connect(G_OBJECT(belvuAlignment), "destroy", G_CALLBACK (onDestroyBelvuAlignment), NULL);
    }
}


/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

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
  
  const char * colorStr = convertColorNumToStr(colorNum);
  getColorFromString(colorStr, result, NULL);
}


/* Draw a single line in the alignment view */
static void drawSingleAlignment(GtkWidget *widget,
                                GdkDrawable *drawable, 
                                BelvuAlignmentProperties *properties,
                                ALN *alnp, 
                                const int lineNum,
                                const gboolean headers)
{
  int x = 0;
  const int y = properties->charHeight * lineNum;
  
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, properties->bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, textColor);
  
  if (headers)
    {
      /* Draw the name */
      if (alnp->name)
        {
          drawText(widget, drawable, gc, x, y, alnp->name);
          x += (properties->bc->maxNameLen * properties->charWidth) + properties->columnPadding;
        }
      
      /* Draw the coords */
      drawIntAsText(widget, drawable, gc, x, y, alnp->start);
      x += (properties->bc->maxStartLen * properties->charWidth) + properties->columnPadding;
      
      drawIntAsText(widget, drawable, gc, x, y, alnp->end);
      x += (properties->bc->maxEndLen * properties->charWidth) + properties->columnPadding;
    }
  else
    {
      /* Draw the sequence */
      if (alnp->seq)
        {
          const int startX = x;
          
          /* Loop through each char and color the background */
          int i = 0;
          for ( ; i < properties->bc->maxLen; ++i)
            {
              GdkColor bgColor;
              findResidueBGcolor(properties->bc, alnp, i, &bgColor);
              gdk_gc_set_foreground(gc, &bgColor);
              
              gdk_draw_rectangle(drawable, gc, TRUE, x, y, properties->charWidth, properties->charHeight);
              x += properties->charWidth;
            }

          gdk_gc_set_foreground(gc, textColor);
          drawText(widget, drawable, gc, startX, y, alnp->seq);
        }
    }
  
  g_object_unref(gc);
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
      drawSingleAlignment(widget, drawable, properties, alnp, i, TRUE);
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
      drawSingleAlignment(widget, drawable, properties, alnp, i, FALSE);
    }
}


static GdkDrawable* createBlankSizedPixmap(GtkWidget *widget, GdkDrawable *window, const int width, const int height)
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
          /* Push the bitmap onto the window, clipping it to the display area
           * (the clipping is necessary because the allocation may be larger
           * than the sequence area's allocation, because the sequence area may
           * have a horizontal scrollbar */
          GdkGC *gc = gdk_gc_new(window);

          GtkAdjustment *vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(properties->headersArea));
          const int height = (vAdjustment->value + properties->seqArea->allocation.height);
          
          GdkRectangle clipRect = {0, 0, widget->allocation.width, height};
          gdk_gc_set_clip_rectangle(gc, &clipRect);
          
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
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
          bitmap = createBlankSizedPixmap(widget, window, properties->seqRect.width, properties->seqRect.height);
          drawBelvuSequence(widget, bitmap, properties);
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

/* Called when the belvu alignment section has been scrolled */
static void onScrollChangedBelvuAlignment(GtkObject *object, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  gtk_widget_queue_draw(belvuAlignment);
}

/***********************************************************
 *                         Sizing                          *
 ***********************************************************/

static void calculateBelvuAlignmentBorders(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  properties->columnPadding = DEFAULT_PADDING_CHARS * properties->charWidth;
  
  /* Calculate the size of the drawing area for the header columns */
  properties->headersRect.x = DEFAULT_XPAD;
  properties->headersRect.y = DEFAULT_YPAD;
  properties->headersRect.width = (properties->bc->maxNameLen * properties->charWidth) + properties->columnPadding +
                                  (properties->bc->maxStartLen * properties->charWidth) + properties->columnPadding +
                                  (properties->bc->maxEndLen * properties->charWidth);
  properties->headersRect.height = properties->bc->alignArr->len * properties->charHeight + (2 * DEFAULT_YPAD);

  /* Calculate the size of the drawing area for the sequences */
  properties->seqRect.x = DEFAULT_XPAD;
  properties->seqRect.y = DEFAULT_YPAD;
  properties->seqRect.width = (properties->bc->maxLen * properties->charWidth);
  properties->seqRect.height = properties->headersRect.height;
  
  gtk_layout_set_size(GTK_LAYOUT(properties->headersArea), properties->headersRect.width, properties->headersRect.height);
  gtk_layout_set_size(GTK_LAYOUT(properties->seqArea), properties->seqRect.width, properties->seqRect.height);
  

  if (properties->headersRect.width != properties->headersArea->allocation.width)
    gtk_widget_set_size_request(properties->headersArea, properties->headersRect.width, -1);
}


static void onSizeAllocateBelvuAlignment(GtkWidget *widget, GtkAllocation *allocation, gpointer data)
{
  GtkWidget *belvuAlignment = GTK_WIDGET(data);
  calculateBelvuAlignmentBorders(belvuAlignment);
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

GtkWidget* createBelvuAlignment(BelvuContext *bc)
{
  /* Create the sequence drawing area and the columns drawing area (for the
   * name and coords etc). The sequence drawing area will have a horizontal 
   * scrollbar and they will both share the same vertical scrollbar */
  GtkWidget *seqArea = gtk_layout_new(NULL, NULL);
  GtkWidget *seqScrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(seqScrollWin), GTK_POLICY_ALWAYS, GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(seqScrollWin), seqArea);
  
  GtkAdjustment *vAdjustment = gtk_layout_get_vadjustment(GTK_LAYOUT(seqArea));
  GtkWidget *headersArea = gtk_layout_new(NULL, vAdjustment);

  /* Wrap everything in an hbox */
  GtkWidget *belvuAlignment = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(belvuAlignment), headersArea, FALSE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(belvuAlignment), seqScrollWin, TRUE, TRUE, 0);
  
  /* Set the properties and connect signals */
  belvuAlignmentCreateProperties(belvuAlignment, bc, headersArea, seqArea);

  g_signal_connect(G_OBJECT(seqArea), "expose-event", G_CALLBACK(onExposeBelvuSequence), belvuAlignment);  
  g_signal_connect(G_OBJECT(headersArea), "expose-event", G_CALLBACK(onExposeBelvuColumns), belvuAlignment);  
  g_signal_connect(G_OBJECT(seqArea), "size-allocate", G_CALLBACK(onSizeAllocateBelvuAlignment), belvuAlignment);

  g_signal_connect(G_OBJECT(vAdjustment), "changed", G_CALLBACK(onScrollChangedBelvuAlignment), belvuAlignment);
  g_signal_connect(G_OBJECT(vAdjustment), "value-changed", G_CALLBACK(onScrollChangedBelvuAlignment), belvuAlignment);
  
  return belvuAlignment;
}





