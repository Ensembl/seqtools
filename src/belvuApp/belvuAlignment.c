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

  GdkRectangle displayRect;         /* The area in which we do the drawing */
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
static void belvuAlignmentCreateProperties(GtkWidget *belvuAlignment, BelvuContext *bc)
{
  if (belvuAlignment)
    {
      BelvuAlignmentProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->columnPadding = 0; /* calculated in calculate-borders */
  
      /* Find a fixed-width font */
      const char *fontFamily = findFixedWidthFont(belvuAlignment);
      PangoFontDescription *fontDesc = pango_font_description_from_string(fontFamily);
      pango_font_description_set_size(fontDesc, pango_font_description_get_size(belvuAlignment->style->font_desc));

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


/* Draw a single line in the alignment view */
static void drawSingleAlignment(GtkWidget *widget,
                                GdkDrawable *drawable, 
                                GdkGC *gc,
                                BelvuAlignmentProperties *properties,
                                ALN *alnp, 
                                const int lineNum)
{
  int x = 0;
  const int y = properties->charHeight * lineNum;
  
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

  /* Draw the sequence */
  if (alnp->seq)
    {
      drawText(widget, drawable, gc, x, y, alnp->seq);
      x += (properties->bc->maxLen * properties->charWidth) + properties->columnPadding;
    }
}


/* Draw the alignment view */
static void drawBelvuAlignment(GtkWidget *widget, GdkDrawable *drawable)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(widget);
  BelvuContext *bc = properties->bc;
  
  GdkGC *gc = gdk_gc_new(drawable);
  GdkColor *textColor = getGdkColor(BELCOLOR_ALIGN_TEXT, bc->defaultColors, FALSE, FALSE);
  gdk_gc_set_foreground(gc, textColor);

  /* Loop through each alignment */
  int i = 0;
  for ( ; i < bc->alignArr->len; ++i)
    {
      ALN *alnp = &g_array_index(bc->alignArr, ALN, i);
      drawSingleAlignment(widget, drawable, gc, properties, alnp, i);
    }
  
  g_object_unref(gc);
}


/* Expose handler for the alignment section */
static gboolean onExposeBelvuAlignment(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  GdkDrawable *window = GTK_LAYOUT(widget)->bin_window;

  if (window)
    {
      drawBelvuAlignment(widget, window);
//      GdkDrawable *bitmap = widgetGetDrawable(widget);
//      
//      if (!bitmap)
//        {
//          /* There isn't a bitmap yet. Create it now. */
//	  bitmap = createBlankPixmap(widget);
//          drawBelvuAlignment(widget, bitmap);
//        }
//      
//      if (bitmap)
//        {
//          /* Push the bitmap onto the window */
//          GdkGC *gc = gdk_gc_new(window);
//          gdk_draw_drawable(window, gc, bitmap, 0, 0, 0, 0, -1, -1);
//          g_object_unref(gc);
//        }
//      else
//	{
//	  g_warning("Failed to draw Belvu alignment [%p] - could not create bitmap.\n", widget);
//	}
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
  
//  widgetClearCachedDrawable(belvuAlignment, NULL);
//  gtk_widget_queue_draw(belvuAlignment);
}

/***********************************************************
 *                         Sizing                          *
 ***********************************************************/

static int calculateBelvuAlignmentWidth(BelvuAlignmentProperties *properties)
{
  /* Calculate the total width of the layout */
  int result = (properties->bc->maxNameLen * properties->charWidth) + properties->columnPadding +
               (properties->bc->maxStartLen * properties->charWidth) + properties->columnPadding +
               (properties->bc->maxEndLen * properties->charWidth) + properties->columnPadding +
               (properties->bc->maxLen * properties->charWidth);
  
  return result;
}


static void calculateBelvuAlignmentBorders(GtkWidget *belvuAlignment)
{
  BelvuAlignmentProperties *properties = belvuAlignmentGetProperties(belvuAlignment);
  
  properties->columnPadding = DEFAULT_PADDING_CHARS * properties->charWidth;
  
  properties->displayRect.x = DEFAULT_XPAD;
  properties->displayRect.y = DEFAULT_YPAD;
  properties->displayRect.width = calculateBelvuAlignmentWidth(properties);
  properties->displayRect.height = properties->bc->alignArr->len * (properties->charHeight + DEFAULT_YPAD);
  
  gtk_layout_set_size(GTK_LAYOUT(belvuAlignment), properties->displayRect.width, properties->displayRect.height);
  
//  widgetClearCachedDrawable(belvuAlignment, NULL);
//  gtk_widget_queue_draw(belvuAlignment);
}


static void onSizeAllocateBelvuAlignment(GtkWidget *belvuAlignment, GtkAllocation *allocation, gpointer data)
{
  calculateBelvuAlignmentBorders(belvuAlignment);
}


/***********************************************************
 *                      Initialisation                     *
 ***********************************************************/

GtkWidget* createBelvuAlignment(BelvuContext *bc)
{
  /* Create the drawing area */
  GtkWidget *belvuAlignment = gtk_layout_new(NULL, NULL);

  /* Put it in a scrolled window */
  GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
  gtk_container_add(GTK_CONTAINER(scrollWin), belvuAlignment);
  
  /* Set the properties and connect signals */
  belvuAlignmentCreateProperties(belvuAlignment, bc);

  g_signal_connect(G_OBJECT(belvuAlignment), "expose-event", G_CALLBACK(onExposeBelvuAlignment), NULL);  
  g_signal_connect(G_OBJECT(belvuAlignment), "size-allocate", G_CALLBACK(onSizeAllocateBelvuAlignment), NULL);

  GtkAdjustment *adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrollWin));
  g_signal_connect(G_OBJECT(adjustment), "changed", G_CALLBACK(onScrollChangedBelvuAlignment), belvuAlignment);
  g_signal_connect(G_OBJECT(adjustment), "value-changed", G_CALLBACK(onScrollChangedBelvuAlignment), belvuAlignment);
  
  return scrollWin;
}





