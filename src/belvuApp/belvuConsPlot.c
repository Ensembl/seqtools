/*  File: belvuConsPlot.h
 *  Author: Gemma Barson, 2011-07-01
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
 * Description: See belvuConsPlot.h
 *----------------------------------------------------------------------------
 */

#include "belvuApp/belvuConsPlot.h"
#include "belvuApp/belvuWindow.h"
#include "seqtoolsUtils/utilities.h"
#include <gtk/gtk.h>

#define DEFAULT_CONS_PLOT_SCALE_HEIGHT      200   /* default height of the conservation plot scale */
#define CONS_PLOT_XPAD                      0     /* Default x padding for the conservation plot */
#define CONS_PLOT_YPAD                      4     /* Default y padding for the conservation plot */
#define TICKMARK_INTERVAL                   10    /* how many values between tickmark labels */
#define MAJOR_TICKMARK_HEIGHT               6     /* height of major tick marks in the sequence area header */
#define MINOR_TICKMARK_HEIGHT               3     /* height of minor tick marks in the sequence area header */



/***********************************************************
 *                         Properties                      *
 ***********************************************************/

/* Properties specific to the belvu window */
typedef struct _ConsPlotProperties
  {
    BelvuContext *bc;                   /* The belvu context */
    GtkActionGroup *actionGroup;
  } ConsPlotProperties;


/* Properties for generic windows */
static ConsPlotProperties* consPlotGetProperties(GtkWidget *consPlot)
{
  return consPlot ? (ConsPlotProperties*)(g_object_get_data(G_OBJECT(consPlot), "BelvuConsPlotProperties")) : NULL;
}


/* Does the job of destroying the window */
static void onDestroyConsPlot(GtkWidget *consPlot)
{
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  
  if (properties)
    {
      /* Free the properties struct itself */
      g_free(properties);
      properties = NULL;
      g_object_set_data(G_OBJECT(consPlot), "BelvuConsPlotProperties", NULL);
    }
}


/* Create the properties struct and initialise all values. */
static void consPlotCreateProperties(GtkWidget *consPlot, 
                                     BelvuContext *bc, 
                                     GtkActionGroup *actionGroup)
{
  if (consPlot)
    {
      ConsPlotProperties *properties = g_malloc(sizeof *properties);
      
      properties->bc = bc;
      properties->actionGroup = actionGroup;

      g_object_set_data(G_OBJECT(consPlot), "BelvuConsPlotProperties", properties);
      g_signal_connect(G_OBJECT(consPlot), "destroy", G_CALLBACK (onDestroyConsPlot), NULL);
    }
}


BelvuContext* consPlotGetContext(GtkWidget *consPlot)
{
  ConsPlotProperties *properties = consPlotGetProperties(consPlot);
  return (properties ? properties->bc : NULL);
}


/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

static void drawConsPlotHeader(GtkWidget *widget, GdkDrawable *drawable, BelvuContext *bc)
{
  IntRange range = {1, bc->maxLen};
  GdkRectangle rect = {CONS_PLOT_XPAD, 0, widget->allocation.width, widget->allocation.height};
  const int charWidth = 8;
  
  drawHScale(widget, 
             drawable,
             &range,
             &rect,
             charWidth,
             TICKMARK_INTERVAL / 2,
             TICKMARK_INTERVAL, 
             MINOR_TICKMARK_HEIGHT, 
             MAJOR_TICKMARK_HEIGHT);
}


static void drawConsPlot(GtkWidget *widget, GdkDrawable *drawable, BelvuContext *bc)
{
}


/***********************************************************
 *                        Events                           *
 ***********************************************************/

static gboolean onExposeConsPlotHeader(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  BelvuContext *bc = (BelvuContext*)data;
  GdkDrawable *window = widget->window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, widget->allocation.width, widget->allocation.height);
          drawConsPlotHeader(widget, bitmap, bc);
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
          g_warning("Failed to draw conservation plot [%p] - could not create bitmap.\n", widget);
        }
    }
  
  return TRUE;
}


static gboolean onExposeConsPlot(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  BelvuContext *bc = (BelvuContext*)data;
  
  GdkDrawable *window = widget->window;
  
  if (window)
    {
      GdkDrawable *bitmap = widgetGetDrawable(widget);
      
      if (!bitmap)
        {
          /* There isn't a bitmap yet. Create it now. */
          bitmap = createBlankSizedPixmap(widget, window, widget->allocation.width, widget->allocation.height);
          drawConsPlot(widget, bitmap, bc);
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
          g_warning("Failed to draw conservation plot [%p] - could not create bitmap.\n", widget);
        }
    }
  
  return TRUE;
}


/***********************************************************
 *                     Initialisation                      *
 ***********************************************************/

static void setConsPlotStyleProperties(GtkWidget *window, BelvuContext *bc)
{
  gtk_widget_set_name(window, BELVU_CONS_PLOT_WINDOW_NAME);
  
  /* Just hide the widget when it is closed rather than destroy it */
  g_signal_connect(window, "delete-event", G_CALLBACK(gtk_widget_hide_on_delete), NULL);
  
  /* Set default size based on scale height and alignment window width */
  GdkScreen *screen = gtk_widget_get_screen(window);
  const int screenWidth = gdk_screen_get_width(screen);
  const int screenHeight = gdk_screen_get_height(screen);
  
  const int width = screenWidth * DEFAULT_BELVU_WINDOW_WIDTH_FRACTION;
  const int height = DEFAULT_CONS_PLOT_SCALE_HEIGHT;
  
  gtk_window_set_default_size(GTK_WINDOW(window), width, height);
  
  /* Set the initial position */
  const int x = (screenWidth - width) / 4;
  const int y = (screenHeight - height) / 4;
  gtk_window_move(GTK_WINDOW(window), x, y);
}


void createConsPlot(BelvuContext *bc)
{
  /* Create the window */
  bc->consPlot = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  setConsPlotStyleProperties(bc->consPlot, bc);
  
  gtk_window_set_title(GTK_WINDOW(bc->consPlot), "Belvu - Conservation Profile");
  
  /* We must add all toplevel windows to the list of spawned windows */
  bc->spawnedWindows = g_slist_prepend(bc->spawnedWindows, bc->consPlot);
  
  /* Create the context menu and set a callback to show it */
  GtkActionGroup *actionGroup = NULL;
  GtkUIManager *uiManager = createUiManager(bc->consPlot, bc, &actionGroup);
  GtkWidget *contextmenu = createBelvuMenu(bc->consPlot, "/PlotContextMenu", uiManager);
  
  gtk_widget_add_events(bc->consPlot, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(bc->consPlot), "button-press-event", G_CALLBACK(onButtonPressBelvu), contextmenu);
  
  /* We'll place everything in a vbox */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(bc->consPlot), vbox);
  
  /* Add the drawing area */
  GtkWidget *drawing = gtk_drawing_area_new();
  gtk_box_pack_start(GTK_BOX(vbox), drawing, TRUE, TRUE, 0);
  g_signal_connect(G_OBJECT(drawing), "expose-event", G_CALLBACK(onExposeConsPlot), bc);
  
  /* Add the header */
  GtkWidget *header = gtk_drawing_area_new();
  gtk_box_pack_start(GTK_BOX(vbox), header, FALSE, TRUE, 0);
  g_signal_connect(G_OBJECT(header), "expose-event", G_CALLBACK(onExposeConsPlotHeader), bc);
  const int charHeight = 16;
  gtk_widget_set_size_request(header, -1, charHeight + CONS_PLOT_YPAD);
  
  /* Set default background color */
  GdkColor *bgColor = getGdkColor(BELCOLOR_BACKGROUND, bc->defaultColors, FALSE, FALSE);
  gtk_widget_modify_bg(drawing, GTK_STATE_NORMAL, bgColor);
  
  /* Set properties */
  consPlotCreateProperties(bc->consPlot, bc, actionGroup);
  
  gtk_widget_show_all(bc->consPlot);
  gtk_window_present(GTK_WINDOW(bc->consPlot));
}
