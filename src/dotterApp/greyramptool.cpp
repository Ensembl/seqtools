/*  File: greyramp.c
 *  Author: Gemma Barson, 2010-08-31
 *  Copyright [2018-2020] EMBL-European Bioinformatics Institute
 *  Copyright (c) 2006-2017 Genome Research Ltd
 * ---------------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
 * Description: Creates the greyramp tool window for a Dotter window.
 *              The greyramp tool is used to control how much noise is filtered
 *              out in the dotplot. The greyramp too is owned by the Dotter
 *              window and will be destroyed when that window is closed.
 *----------------------------------------------------------------------------
 */

#include <gtk/gtk.h>
#include <seqtoolsUtils/utilities.hpp>
#include <dotterApp/dotter_.hpp>

//#define MAXY                            20
//#define MINY                            140
//#define BORDER                          20
#define GREYRAMP_MIN                    0       /* min possible value to set the black/white point to */
#define GREYRAMP_MAX			255	/* max possible value to set the black/white point to */
#define GRADIENT_RECT_WIDTH             256     /* width of the gradient drawing area */
#define GRADIENT_RECT_HEIGHT            100     /* height of the gradient drawing area */
#define GRADIENT_RECT_HEIGHT_MIN        20      /* height of the gradient drawing area when "minimised" */
#define GRADIENT_RECT_MARKER_HEIGHT     10      /* height of the triangular markers on the gradient area */
#define GRADIENT_RECT_X_PADDING         10      /* x padding around the gradient area */
#define GRADIENT_RECT_Y_PADDING         10      /* y padding around the gradient area */
#define GRADIENT_RECT_X_PADDING_MIN     2       /* x padding around the gradient area when minimised */
#define GRADIENT_RECT_Y_PADDING_MIN     2       /* y padding around the gradient area when minimised */
#define GRADIENT_RECT_FRAME_PADDING     10      /* padding around the outside of the gradient widget */



typedef struct _CallbackItem
{
  GtkWidget *widget;
  GtkCallback func;
} CallbackItem;


class GreyrampProperties
{
public:
  GtkWidget *widget;
  GtkWidget *greyrampWindow;          /* the toplevel window the greyrampTool will be in IF
                                       * undocked from the main window */
  DotterWindowContext *dwc;
  GdkRectangle gradientRect;          /* the area where the gradient ectangle is drawn */

  GtkWidget *whiteSpinButton;         /* spin button to control the position of the white point in the gradient */
  GtkWidget *blackSpinButton;         /* spin button to control the position of the black point in the gradient */

  int blackPoint;                     /* current position for the black point */
  int whitePoint;                     /* current position for the white point */
  int lastBlackPoint;                 /* the previous value for the black-point spin button (so we can undo) */
  int lastWhitePoint;                 /* the previous value for the white-point spin button (so we can undo) */

  gboolean swapValues;                /* top and bottom values are swapped (i.e. so top is max and bottom is min) */

  gboolean draggingWhite;             /* set to true when dragging the white-point marker */
  gboolean draggingBlack;             /* set to true when dragging the black-point marker */
  gboolean draggingThreshold;         /* set to true when dragging the threshold marker */
  gint dragXPos;                      /* x position of where the drag started */

  GSList *callbackItems;              /* a list of callbacks to be called when the greyramp is updated */

  int greyrampMinDist;                /* minimum dist between black and white point */
};


/***********************************************************
 *                       Properties                        *
 ***********************************************************/

static GreyrampProperties* greyrampGetProperties(GtkWidget *greyramp)
{
  return greyramp ? (GreyrampProperties*)(g_object_get_data(G_OBJECT(greyramp), "GreyrampProperties")) : NULL;
}

static void onDestroyGreyramp(GtkWidget *greyramp)
{
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  if (properties)
    {
      if (properties->callbackItems)
        {
          GSList *item = properties->callbackItems;
          for ( ; item; item = item->next)
            {
              CallbackItem *callbackItem = (CallbackItem*)(item->data);
              delete callbackItem;
            }

          g_slist_free(properties->callbackItems);
        }

      delete properties;
      properties = NULL;
      g_object_set_data(G_OBJECT(greyramp), "GreyrampProperties", NULL);
    }
}

static void greyrampCreateProperties(GtkWidget *greyramp,
                                     GtkWidget *greyrampWindow,
				     DotterWindowContext *dwc,
                                     GdkRectangle *gradientRect,
                                     GtkWidget *whiteSpinButton,
                                     GtkWidget *blackSpinButton,
                                     const int blackPoint,
                                     const int whitePoint,
                                     const int lastBlackPoint,
                                     const int lastWhitePoint,
                                     const gboolean swapValues,
                                     const int greyrampMinDist)
{
  if (greyramp)
    {
      GreyrampProperties *properties = new GreyrampProperties;

      properties->widget = greyramp;
      properties->greyrampWindow = greyrampWindow;
      properties->dwc = dwc;
      properties->gradientRect.x = gradientRect->x;
      properties->gradientRect.y = gradientRect->y;
      properties->gradientRect.width = gradientRect->width;
      properties->gradientRect.height = gradientRect->height;
      properties->whiteSpinButton = whiteSpinButton;
      properties->blackSpinButton = blackSpinButton;
      properties->blackPoint = blackPoint;
      properties->whitePoint = whitePoint;
      properties->lastBlackPoint = lastBlackPoint;
      properties->lastWhitePoint = lastWhitePoint;
      properties->swapValues = swapValues;
      properties->draggingWhite = FALSE;
      properties->draggingBlack = FALSE;
      properties->draggingThreshold = FALSE;
      properties->dragXPos = 0;
      properties->callbackItems = NULL;
      properties->greyrampMinDist = greyrampMinDist;

      g_object_set_data(G_OBJECT(greyramp), "GreyrampProperties", properties);
      g_signal_connect(G_OBJECT(greyramp), "destroy", G_CALLBACK(onDestroyGreyramp), NULL);
    }
}


/* Functions to set the white/black points. Does bounds checking and updates. */
static void greyrampSetWhitePoint(GtkWidget *greyramp, GreyrampProperties *properties, const int whitePoint)
{
  properties->whitePoint = whitePoint;

  const int greyrampMin = GREYRAMP_MIN;
  const int greyrampMax = GREYRAMP_MAX;

  if (properties->whitePoint < greyrampMin)
    {
      properties->whitePoint = greyrampMin;
    }

  if (properties->whitePoint > greyrampMax)
    {
      properties->whitePoint = greyrampMax;
    }

  if (properties->blackPoint - properties->whitePoint < properties->greyrampMinDist)
    {
      properties->whitePoint = properties->blackPoint - properties->greyrampMinDist;
    }

  if (properties->whiteSpinButton)
    {
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->whiteSpinButton), properties->whitePoint);
    }
  else
    {
      gtk_widget_queue_draw(greyramp);
      updateGreyMap(greyramp);
    }
}

static void greyrampSetBlackPoint(GtkWidget *greyramp, GreyrampProperties *properties, const int blackPoint)
{
  properties->blackPoint = blackPoint;

  const int greyrampMin = GREYRAMP_MIN;
  const int greyrampMax = GREYRAMP_MAX;

  if (properties->blackPoint < greyrampMin)
    {
      properties->blackPoint = greyrampMin;
    }

  if (properties->blackPoint > greyrampMax)
    {
      properties->blackPoint = greyrampMax;
    }

  if (properties->blackPoint - properties->whitePoint < properties->greyrampMinDist)
    {
      properties->blackPoint = properties->whitePoint + properties->greyrampMinDist;
    }

  if (properties->blackSpinButton)
    {
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->blackSpinButton), properties->blackPoint);
    }
  else
    {
      gtk_widget_queue_draw(greyramp);
      updateGreyMap(greyramp);
    }
}


/* This registers callbacks so that the given function is called on the given widget whenever the
 * greyramp changes */
void registerGreyrampCallback(GtkWidget *greyramp, GtkWidget *widget, GtkCallback func)
{
  DEBUG_ENTER("registerGreyrampCallback");

  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  CallbackItem *callbackItem = new CallbackItem;
  callbackItem->widget = widget;
  callbackItem->func = func;

  properties->callbackItems = g_slist_append(properties->callbackItems, callbackItem);

  DEBUG_EXIT("registerGreyrampCallback returning ");
}


/***********************************************************
 *                         Utilities                       *
 ***********************************************************/

static gboolean onDeleteGreyrampTool(GtkWidget *widget, GdkEvent *event, gpointer data)
{
  GtkWidget *greyrampTool = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyrampTool);

  if (properties && properties->dwc)
    dotterSetToggleMenuStatus(properties->dwc, "ToggleGreyramp", FALSE);

  return TRUE;
}


/* This should be called whenever the black- or white- point has changed. It causes
 * the greymap for all graphs that have one to be updated. */
void updateGreyMap(GtkWidget *greyramp)
{
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  unsigned char *ramp = (unsigned char*)g_malloc(256 * sizeof(unsigned char));

  int whitePoint = properties->whitePoint;
  int blackPoint = properties->blackPoint;

  /* If black and white point are the same, make the white point less than the black point
   * (But make sure it's still in bounds) */
  if (blackPoint == whitePoint && blackPoint < GREYRAMP_MAX)
    {
      ++blackPoint;
    }
  else if (blackPoint == whitePoint)
    {
      --whitePoint;
    }

  if (blackPoint < whitePoint)
    {
      int i = 0;

      for (i = 0 ; i < blackPoint ; ++i)
	ramp[i] = 0x0 ;

      gdouble fac = 0xff / (gdouble)(whitePoint - blackPoint) ;

      for ( ; i < whitePoint && i < 256 ; ++i)
	ramp[i] =  fac * (i - blackPoint) ;

      for ( ; i < 256 ; ++i)
	ramp[i] = 0xff ;
    }
  else
    {
      int i = 0;

      for (i = 0 ; i < whitePoint ; ++i)
	ramp[i] = 0xff ;

      gdouble fac = 0xff / (gdouble)(blackPoint - whitePoint) ;

      for ( ; i < blackPoint && i < 256 ; ++i)
	ramp[i] = fac * (blackPoint - i) ;

      for ( ; i < 256 ; ++i)
	ramp[i] = 0x0 ;
    }


  /* Loop through our list of widgets that require their callbacks to be called */
  GSList *item = properties->callbackItems;

  for ( ; item; item = item->next)
    {
      CallbackItem *callbackItem = (CallbackItem*)(item->data);
      callbackItem->func(callbackItem->widget, ramp);
    }

  g_free(ramp);
}


/* Determine the extent of the white marker and return the coords in the given result rectangle */
static void getWhiteMarkerRect(GreyrampProperties *properties, GdkRectangle *result)
{
  result->x = properties->gradientRect.x + properties->whitePoint - (GRADIENT_RECT_MARKER_HEIGHT / 2);
  result->y = properties->gradientRect.y - GRADIENT_RECT_MARKER_HEIGHT - (GRADIENT_RECT_Y_PADDING / 2);
  result->width = GRADIENT_RECT_MARKER_HEIGHT;
  result->height = GRADIENT_RECT_MARKER_HEIGHT;
}

/* Determine the extent of the black marker and return the coords in the given result rectangle */
static void getBlackMarkerRect(GreyrampProperties *properties, GdkRectangle *result)
{
  result->x = properties->gradientRect.x + properties->blackPoint - (GRADIENT_RECT_MARKER_HEIGHT / 2);
  result->y = properties->gradientRect.y + properties->gradientRect.height + (GRADIENT_RECT_Y_PADDING / 2);
  result->width = GRADIENT_RECT_MARKER_HEIGHT;
  result->height = GRADIENT_RECT_MARKER_HEIGHT;
}

/* Determine the extent of the threshold marker and return the coords in the given result rectangle */
static void getThresholdMarkerRect(GreyrampProperties *properties, GdkRectangle *result)
{
  result->x = properties->gradientRect.x + properties->whitePoint + (properties->blackPoint - properties->whitePoint) / 2 - (GRADIENT_RECT_MARKER_HEIGHT / 2);
  result->y = properties->gradientRect.y + (properties->gradientRect.height / 2) - (GRADIENT_RECT_MARKER_HEIGHT / 2);
  result->width = GRADIENT_RECT_MARKER_HEIGHT;
  result->height = GRADIENT_RECT_MARKER_HEIGHT;
}


/* Determine whether the given point is inside the white marker */
static gboolean pointInWhiteMarker(GreyrampProperties *properties, const int x, const int y)
{
  GdkRectangle markerRect;
  getWhiteMarkerRect(properties, &markerRect);
  return pointInRect(x, y, &markerRect);
}

/* Determine whether the given point is inside the black marker */
static gboolean pointInBlackMarker(GreyrampProperties *properties, const int x, const int y)
{
  GdkRectangle markerRect;
  getBlackMarkerRect(properties, &markerRect);
  return pointInRect(x, y, &markerRect);
}

/* Determine whether the given point is inside the threshold marker */
static gboolean pointInThresholdMarker(GreyrampProperties *properties, const int x, const int y)
{
  GdkRectangle markerRect;
  getThresholdMarkerRect(properties, &markerRect);
  return pointInRect(x, y, &markerRect);
}


/***********************************************************
 *                         Drawing                         *
 ***********************************************************/

/* Draw the gradient rectangle for the greyramp tool */
static void drawGradient(GdkDrawable *drawable, GtkWidget *greyramp)
{
  GreyrampProperties *properties = greyrampGetProperties(greyramp);
  GdkRectangle *rect = &properties->gradientRect;

  /* Create a gradient */
  cairo_pattern_t *pattern = cairo_pattern_create_linear(rect->x, rect->y, rect->x + rect->width, rect->y);

  /* If the positions are the same, make the black pos slightly smaller (to emulate historic behaviour) */
  int whitePoint = properties->whitePoint;
  int blackPoint = properties->blackPoint;

  if (whitePoint == blackPoint)
    {
      if (blackPoint <= 1)
        {
          whitePoint += 1;
        }
      else
        {
          blackPoint -= 1;
        }
    }

  /* Get the white and black points as fractions of 1.0 */
  gdouble blackPos = blackPoint / ((gdouble)GREYRAMP_MAX);
  gdouble whitePos = whitePoint / ((gdouble)GREYRAMP_MAX);

  /* Add a stop at the start and end. The start is white if the whitePos is less than blackPos */
  const int whiteEnd = whitePos < blackPos ? 0 : 1;
  const int blackEnd = whitePos < blackPos ? 1 : 0;
  cairo_pattern_add_color_stop_rgb(pattern, whiteEnd, 1.0, 1.0, 1.0);
  cairo_pattern_add_color_stop_rgb(pattern, blackEnd, 0.0, 0.0, 0.0);

  /* Add stop points at the black/white positions */
  cairo_pattern_add_color_stop_rgb(pattern, whitePos, 1.0, 1.0, 1.0);
  cairo_pattern_add_color_stop_rgb(pattern, blackPos, 0.0, 0.0, 0.0);

  if (cairo_pattern_status(pattern) == CAIRO_STATUS_SUCCESS)
    {
      cairo_t *cr = gdk_cairo_create(drawable);

      gdk_cairo_rectangle(cr, rect);
      cairo_set_source(cr, pattern);
      cairo_fill(cr);

      cairo_destroy(cr);
    }

  cairo_pattern_destroy(pattern);
}


/* Draw the marker that indicates where the white point is */
static void drawWhiteMarker(GdkDrawable *drawable, GtkWidget *greyramp)
{
  GreyrampProperties *properties = greyrampGetProperties(greyramp);
  DotterContext *dc = properties->dwc->dotterCtx;

  GdkRectangle markerRect;
  getWhiteMarkerRect(properties, &markerRect);

  /* Draw a triangle. Line color is the 'selected' color if dragging */
  GdkColor *fillColor = getGdkColor(DOTCOLOR_MARKER_FILL, dc->defaultColors, FALSE, properties->dwc->usePrintColors);
  GdkColor *lineColor = getGdkColor(DOTCOLOR_MARKER_LINE, dc->defaultColors, properties->draggingWhite, properties->dwc->usePrintColors);

  int numPoints = 3;
  GdkPoint points[numPoints];
  points[0].x = markerRect.x;
  points[0].y = markerRect.y;
  points[1].x = markerRect.x + markerRect.width;
  points[1].y = markerRect.y;
  points[2].x = markerRect.x + (markerRect.width / 2);
  points[2].y = markerRect.y + markerRect.height;

  /* Draw and outline of marker */
  GdkGC *gc = gdk_gc_new(drawable);

  gdk_gc_set_foreground(gc, fillColor);
  gdk_draw_polygon(drawable, gc, TRUE, points, numPoints);
  gdk_gc_set_foreground(gc, lineColor);
  gdk_draw_polygon(drawable, gc, FALSE, points, numPoints);

  g_object_unref(gc);
}


/* Draw the marker that indicates where the black point is */
static void drawBlackMarker(GdkDrawable *drawable, GtkWidget *greyramp)
{
  GreyrampProperties *properties = greyrampGetProperties(greyramp);
  DotterContext *dc = properties->dwc->dotterCtx;

  GdkRectangle markerRect;
  getBlackMarkerRect(properties, &markerRect);

  /* Draw a triangle. Line color is the 'selected' color if dragging */
  GdkColor *fillColor = getGdkColor(DOTCOLOR_MARKER_FILL, dc->defaultColors, FALSE, properties->dwc->usePrintColors);
  GdkColor *lineColor = getGdkColor(DOTCOLOR_MARKER_LINE, dc->defaultColors, properties->draggingBlack, properties->dwc->usePrintColors);

  int numPoints = 3;
  GdkPoint points[numPoints];
  points[0].x = markerRect.x + (markerRect.width / 2);
  points[0].y = markerRect.y;
  points[1].x = markerRect.x;
  points[1].y = markerRect.y + markerRect.height;
  points[2].x = markerRect.x + markerRect.width;
  points[2].y = markerRect.y + markerRect.height;

  /* Draw fill in white and line in black (or green if dragging) */
  GdkGC *gc = gdk_gc_new(drawable);

  gdk_gc_set_foreground(gc, fillColor);
  gdk_draw_polygon(drawable, gc, TRUE, points, numPoints);
  gdk_gc_set_foreground(gc, lineColor);
  gdk_draw_polygon(drawable, gc, FALSE, points, numPoints);

  g_object_unref(gc);
}


/* Draw the marker that indicates where the threshold is */
static void drawThresholdMarker(GdkDrawable *drawable, GtkWidget *greyramp)
{
  GreyrampProperties *properties = greyrampGetProperties(greyramp);
  DotterContext *dc = properties->dwc->dotterCtx;

  GdkRectangle markerRect;
  getThresholdMarkerRect(properties, &markerRect);

  /* Draw the threshold marker outline (there's no fill because we want the background greyramp to show through) */
  GdkGC *gc = gdk_gc_new(drawable);

  GdkColor *lineColor = getGdkColor(DOTCOLOR_THRESHOLD_MARKER, dc->defaultColors, properties->draggingThreshold, properties->dwc->usePrintColors);
  gdk_gc_set_foreground(gc, lineColor);
  gdk_draw_rectangle(drawable, gc, FALSE, markerRect.x, markerRect.y, markerRect.width, markerRect.height);

  g_object_unref(gc);
}


/***********************************************************
 *                       Events                        *
 ***********************************************************/

/* Handler for when the 'Close' button is pressed. Hides the greyramp window. */
static void onCloseGreyramp(GtkWidget *greyramp, gpointer data)
{
  GtkWidget *greyrampTool = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyrampTool);

  if (properties && properties->dwc)
    dotterSetToggleMenuStatus(properties->dwc, "ToggleGreyramp", FALSE);
}


/* Expose function for the gradient area */
static gboolean onExposeGradient(GtkWidget *gradient, GdkEventExpose *event, gpointer data)
{
  GdkWindow *drawable = GTK_LAYOUT(gradient)->bin_window;
  GtkWidget *greyramp = GTK_WIDGET(data);

  if (drawable && greyramp)
    {
      drawGradient(drawable, greyramp);
      drawWhiteMarker(drawable, greyramp);
      drawBlackMarker(drawable, greyramp);
      drawThresholdMarker(drawable, greyramp);
    }

  return TRUE;
}


/* Expose function for the minimised version of the gradient area */
static gboolean onExposeGradientMinimised(GtkWidget *gradient, GdkEventExpose *event, gpointer data)
{
  GdkWindow *drawable = GTK_LAYOUT(gradient)->bin_window;
  GtkWidget *greyramp = GTK_WIDGET(data);

  if (drawable && greyramp)
    {
      drawGradient(drawable, greyramp);
      drawThresholdMarker(drawable, greyramp);
    }

  return TRUE;
}


/* Mouse-button press handler */
static gboolean onButtonPressGradient(GtkWidget *gradient, GdkEventButton *event, gpointer data)
{
  gboolean handled = FALSE;

  if (event->button == 1)
    {
      GtkWidget *greyramp = GTK_WIDGET(data);
      GreyrampProperties *properties = greyrampGetProperties(greyramp);

      /* See if we clicked inside one of the markers */
      if (pointInWhiteMarker(properties, event->x, event->y))
        {
          properties->draggingWhite = TRUE;
          properties->dragXPos = event->x;
        }
      else if (pointInBlackMarker(properties, event->x, event->y))
        {
          properties->draggingBlack = TRUE;
          properties->dragXPos = event->x;
        }
      else if (pointInThresholdMarker(properties, event->x, event->y))
        {
          properties->draggingThreshold = TRUE;
          properties->dragXPos = event->x;
        }

      if (properties->draggingWhite || properties->draggingBlack || properties->draggingThreshold)
        {
          /* Remember the current values so we can undo if this is the start of a drag */
          properties->lastBlackPoint =  properties->blackPoint;
          properties->lastWhitePoint =  properties->whitePoint;

          gtk_widget_queue_draw(greyramp);
        }

      handled = TRUE;
    }

  return handled;
}

/* Mouse-button release handler */
static gboolean onButtonReleaseGradient(GtkWidget *gradient, GdkEventButton *event, gpointer data)
{
  GtkWidget *greyramp = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  properties->draggingWhite = FALSE;
  properties->draggingBlack = FALSE;
  properties->draggingThreshold = FALSE;
  properties->dragXPos = 0;

  gtk_widget_queue_draw(greyramp);

  return TRUE;
}


/* Mouse-motion event handler. */
static gboolean onMouseMoveGradient(GtkWidget *gradient, GdkEventMotion *event, gpointer data)
{
  gboolean handled = FALSE;

  if (event->state & GDK_BUTTON1_MASK) /* left button pressed */
    {
      GtkWidget *greyramp = GTK_WIDGET(data);
      GreyrampProperties *properties = greyrampGetProperties(greyramp);

      if (properties->draggingWhite)
        {
          /* Move the white point by the amount offset since the last move/button-press */
          const int offset = event->x - properties->dragXPos;
          properties->dragXPos = event->x;
          greyrampSetWhitePoint(greyramp, properties, properties->whitePoint + offset);
        }
      else if (properties->draggingBlack)
        {
          /* Move the black point */
          const int offset = event->x - properties->dragXPos;
          properties->dragXPos = event->x;
          greyrampSetBlackPoint(greyramp, properties, properties->blackPoint + offset);
        }
      else if (properties->draggingThreshold)
        {
          /* Move both points together */
          const int offset = event->x - properties->dragXPos;
          properties->dragXPos = event->x;
          greyrampSetWhitePoint(greyramp, properties, properties->whitePoint + offset);
          greyrampSetBlackPoint(greyramp, properties, properties->blackPoint + offset);
        }

      handled = TRUE;
    }

  return handled;
}


/* Called when the black-point spin button's value has changed */
static void onBlackSpinButtonChanged(GtkSpinButton *blackSpinButton, gpointer data)
{
  GtkWidget *greyramp = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  properties->blackPoint = gtk_spin_button_get_value_as_int(blackSpinButton);

  gtk_widget_queue_draw(greyramp);
  updateGreyMap(greyramp);
}


/* Called when the white-point spin button's value has changed */
static void onWhiteSpinButtonChanged(GtkSpinButton *whiteSpinButton, gpointer data)
{
  GtkWidget *greyramp = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  properties->whitePoint = gtk_spin_button_get_value_as_int(whiteSpinButton);

  gtk_widget_queue_draw(greyramp);
  updateGreyMap(greyramp);
}


/* Undo: this restores the previous spin button values that were saved before the last drag. */
static gint onPressUndoButton(GtkWidget *button, gpointer data)
{
  GtkWidget *greyramp = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  int temp = properties->blackPoint;
  properties->blackPoint = properties->lastBlackPoint;
  properties->lastBlackPoint = temp;

  temp = properties->whitePoint;
  properties->whitePoint = properties->lastWhitePoint;
  properties->lastWhitePoint = temp;

  if (properties->blackSpinButton)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->blackSpinButton), (float)properties->blackPoint) ;

  if (properties->whiteSpinButton)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->whiteSpinButton), (float)properties->whitePoint);

  return TRUE;
}


/* Swap the min and max values - has the effect of inverting the colors */
static gint onPressSwapButton(GtkWidget *button, gpointer data)
{
  GtkWidget *greyramp = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyramp);

  int temp = properties->blackPoint;
  properties->blackPoint = properties->whitePoint;
  properties->whitePoint = temp;

  properties->swapValues = !properties->swapValues;

  if (properties->blackSpinButton)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->blackSpinButton), (float)properties->blackPoint) ;

  if (properties->whiteSpinButton)
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(properties->whiteSpinButton), (float)properties->whitePoint);

  return TRUE;
}


static void onMaximiseGreyrampClicked(GtkWidget *button, gpointer data)
{
  GtkWidget *greyrampTool = GTK_WIDGET(data);
  GreyrampProperties *properties = greyrampGetProperties(greyrampTool);

  dotterSetToggleMenuStatus(properties->dwc, "ToggleGreyramp", TRUE);
}


/***********************************************************
 *                    Initialisation                       *
 ***********************************************************/

/* Utility to create a spin button with the range 0 to 255 and various other default values. */
static GtkWidget* createSpinButton(const gdouble initValue, GCallback callbackFunc, GtkWidget *greyramp)
{
  GtkObject *adjustment = gtk_adjustment_new(initValue,
                                             0.0,       /* lower */
                                             255.0,     /* upper */
                                             1.0,       /* step increment */
                                             1.0,       /* page increment */
                                             0.0);      /* page size */

  GtkWidget *spinButton = gtk_spin_button_new(GTK_ADJUSTMENT(adjustment),
                                              0.5, /* climb rate, i.e. amount to increase per button click */
                                              0);  /* num decimal places to display */

  g_signal_connect(spinButton,"value_changed", callbackFunc, greyramp);

  return spinButton;
}


/* Create the toolbar for the greyramp window. */
static GtkWidget* createGreyrampToolbar(GtkTable *table,
                                        const int col,
                                        GtkWidget *greyramp,
                                        GtkWidget *whiteSpinButton,
                                        GtkWidget *blackSpinButton,
                                        const int blackPoint,
                                        const int whitePoint)
{
  GtkWidget *quitButton = gtk_button_new_with_label("Close");
  GtkWidget *swapButton = gtk_button_new_with_label("Swap");
  GtkWidget *undoButton = gtk_button_new_with_label("Undo");

  gtk_widget_set_tooltip_text(quitButton, "Close the full greyramp tool and show the minimised version instead (Ctrl-G)");
  gtk_widget_set_tooltip_text(swapButton, "Swap the black and white point to invert the colours");
  gtk_widget_set_tooltip_text(undoButton, "Undo/redo the last change to the sliders");

  g_signal_connect(G_OBJECT(quitButton), "pressed", G_CALLBACK(onCloseGreyramp), greyramp);
  g_signal_connect(G_OBJECT(undoButton), "pressed", G_CALLBACK(onPressUndoButton), greyramp);
  g_signal_connect(G_OBJECT(swapButton), "pressed", G_CALLBACK(onPressSwapButton), greyramp);

  /* Pack them vertically into the 'toolbar' (which is just a vbox for now) */
  GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_border_width(GTK_CONTAINER(vbox), 5);

  const int padding = 5;
  int row = 0;

  gtk_table_attach(table, whiteSpinButton, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, padding, padding);
  ++row;
  gtk_table_attach(table, quitButton, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, padding, padding);
  ++row;
  gtk_table_attach(table, swapButton, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, padding, padding);
  ++row;
  gtk_table_attach(table, undoButton, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, padding, padding);
  ++row;
  gtk_table_attach(table, blackSpinButton, col, col + 1, row, row + 1, GTK_SHRINK, GTK_SHRINK, padding, padding);
  ++row;

  /* Make sure neither spin button is focused at the start because if it is
   * then its text will be selected and we will inadvertently overwrite the
   * contents of the primary clipboard. */
  gtk_container_set_focus_child(GTK_CONTAINER(vbox), quitButton);

  return vbox;
}


/* Create the rectangle area that displays the gradient */
static GtkWidget* createGradientRect(GtkWidget *greyramp, GdkRectangle *rect)
{
  /* Put it in a frame with a border */
  GtkWidget *frame = gtk_frame_new(NULL);
  gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_ETCHED_IN);
  gtk_container_set_border_width(GTK_CONTAINER(frame), GRADIENT_RECT_FRAME_PADDING);

  /* Get the total size of the gradient area, including markers and padding */
  const int totalWidth = GRADIENT_RECT_WIDTH + (2 * GRADIENT_RECT_X_PADDING) ;
  const int totalHeight = GRADIENT_RECT_HEIGHT + (2 * (GRADIENT_RECT_MARKER_HEIGHT + GRADIENT_RECT_Y_PADDING));

  /* We'll draw the gradient and the child markers onto a layout */
  GtkWidget *layout = gtk_layout_new(NULL, NULL);
  gtk_layout_set_size(GTK_LAYOUT(layout), totalWidth, totalHeight);
  gtk_widget_set_size_request(layout, totalWidth, totalHeight);
  gtk_container_add(GTK_CONTAINER(frame), layout);

  gtk_widget_add_events(layout, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(layout, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(layout, GDK_POINTER_MOTION_MASK);

  g_signal_connect(G_OBJECT(layout), "expose-event", G_CALLBACK(onExposeGradient), greyramp);
  g_signal_connect(G_OBJECT(layout), "button-press-event", G_CALLBACK(onButtonPressGradient), greyramp);
  g_signal_connect(G_OBJECT(layout), "button-release-event", G_CALLBACK(onButtonReleaseGradient), greyramp);
  g_signal_connect(G_OBJECT(layout), "motion-notify-event", G_CALLBACK(onMouseMoveGradient), greyramp);

  /* Set the size of the gradient rectangle to be drawn */
  rect->x = GRADIENT_RECT_X_PADDING;
  rect->y = GRADIENT_RECT_MARKER_HEIGHT + GRADIENT_RECT_Y_PADDING;
  rect->width = GRADIENT_RECT_WIDTH;
  rect->height = GRADIENT_RECT_HEIGHT;

  return frame;
}


/* Create the minimised version of the rectangle area that displays the gradient */
static GtkWidget* createGradientRectMinimised(GdkRectangle *rect, GtkWidget *greyrampTool)
{
  /* Get the total size of the gradient area, including markers and padding */
  const int totalWidth = GRADIENT_RECT_WIDTH + (2 * GRADIENT_RECT_X_PADDING) ;
  const int totalHeight = GRADIENT_RECT_HEIGHT_MIN + (2 * GRADIENT_RECT_Y_PADDING_MIN);

  /* We'll draw the gradient and the child markers onto a layout */
  GtkWidget *greyramp = gtk_layout_new(NULL, NULL);
  gtk_layout_set_size(GTK_LAYOUT(greyramp), totalWidth, totalHeight);
  gtk_widget_set_size_request(greyramp, totalWidth, totalHeight);

  gtk_widget_add_events(greyramp, GDK_BUTTON_PRESS_MASK);
  gtk_widget_add_events(greyramp, GDK_BUTTON_RELEASE_MASK);
  gtk_widget_add_events(greyramp, GDK_POINTER_MOTION_MASK);

  g_signal_connect(G_OBJECT(greyramp), "expose-event", G_CALLBACK(onExposeGradientMinimised), greyrampTool);
  g_signal_connect(G_OBJECT(greyramp), "button-press-event", G_CALLBACK(onButtonPressGradient), greyrampTool);
  g_signal_connect(G_OBJECT(greyramp), "button-release-event", G_CALLBACK(onButtonReleaseGradient), greyrampTool);
  g_signal_connect(G_OBJECT(greyramp), "motion-notify-event", G_CALLBACK(onMouseMoveGradient), greyrampTool);

  /* Set the size of the gradient rectangle to be drawn */
  rect->x = GRADIENT_RECT_X_PADDING_MIN;
  rect->y = GRADIENT_RECT_Y_PADDING_MIN;
  rect->width = GRADIENT_RECT_WIDTH;
  rect->height = GRADIENT_RECT_HEIGHT_MIN;

  return greyramp;
}


/* Create a window to hold the greyramp tool when it is un-docked */
static GtkWidget *createGreyrampToolWindow(DotterWindowContext *dwc, GtkWidget *greyrampTool)
{
  GtkWidget *greyrampWindow = gtk_window_new(GTK_WINDOW_TOPLEVEL);

  char *title = g_strdup_printf("%sGreyramp Tool", dotterGetTitlePrefix(dwc->dotterCtx));
  gtk_window_set_title(GTK_WINDOW(greyrampWindow), title);
  g_free(title);

  /* Set event handlers */
  gtk_widget_add_events(greyrampTool, GDK_BUTTON_PRESS_MASK);
  g_signal_connect(G_OBJECT(greyrampWindow), "delete-event", G_CALLBACK(onDeleteGreyrampTool), greyrampTool);

  return greyrampWindow;
}


/* Create the minimised version of the greyramp widget. */
GtkWidget* createGreyrampToolMinimised(DotterWindowContext *dwc,
                                       const int whitePoint,
                                       const int blackPoint)
{
  DEBUG_ENTER("createGreyrampToolMinimised");

  GtkWidget *greyrampTool = gtk_hbox_new(FALSE, 0) ;

  /* Create a layout for drawing the greyramp gradient onto. This will be in the first column,
   * spanning all rows */
  GdkRectangle gradientRect;
  GtkWidget *gradient = createGradientRectMinimised(&gradientRect, greyrampTool);
  gtk_box_pack_start(GTK_BOX(greyrampTool), gradient, TRUE, TRUE, 0);

  /* Create a 'maximise' button to toggle to the full sized greyramp tool */
  GtkWidget *button = gtk_button_new_with_label("Maximise");
  gtk_box_pack_start(GTK_BOX(greyrampTool), button, FALSE, FALSE, 0);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(onMaximiseGreyrampClicked), greyrampTool);
  gtk_widget_set_tooltip_text(button, "Maximise the greyramp tool for more options (Ctrl-G)");

  greyrampCreateProperties(greyrampTool,
                           NULL,
			   dwc,
                           &gradientRect,
                           NULL,
                           NULL,
                           blackPoint,
                           whitePoint,
                           0,
                           255,
                           FALSE,
                           blackPoint - whitePoint);

  gtk_widget_show_all(greyrampTool);

  updateGreyMap(greyrampTool);

  DEBUG_EXIT("createGreyrampToolMinimised returning ");
  return greyrampTool;
}


/* Create the greyramp widget. Pass the initial value for the top and bottom spin buttons. If
 * 'swapValues' is true these will be swapped. */
GtkWidget* createGreyrampTool(DotterWindowContext *dwc,
                              const int whitePointIn,
                              const int blackPointIn,
                              const gboolean swapValues,
                              GtkWidget **greyrampWindow_out)
{
  DEBUG_ENTER("createGreyrampTool");

  /* Create the window */
  GtkWidget *greyrampTool = gtk_frame_new(NULL);
  gtk_frame_set_shadow_type(GTK_FRAME(greyrampTool), GTK_SHADOW_ETCHED_IN);

  /* Outer container is a table */
  const int numRows = 5;
  const int numCols = 2;
  const int padding = 0;
  GtkTable *table = GTK_TABLE(gtk_table_new(numRows, numCols, FALSE));
  gtk_container_add(GTK_CONTAINER(greyrampTool), GTK_WIDGET(table));

  /* Create a layout for drawing the greyramp gradient onto. This will be in the first column,
   * spanning all rows */
  GdkRectangle gradientRect;
  GtkWidget *gradientWidget = createGradientRect(greyrampTool, &gradientRect);
  gtk_table_attach(table, gradientWidget, 0, 1, 0, numRows, GTK_SHRINK, GTK_SHRINK, padding, padding);

  /* Create the toolbar and pack it into the parent */
  const int whitePoint = swapValues ? blackPointIn : whitePointIn;
  const int blackPoint = swapValues ? whitePointIn : blackPointIn;

  GtkWidget *whiteSpinButton = createSpinButton(whitePoint + 1, G_CALLBACK(onWhiteSpinButtonChanged), greyrampTool);
  GtkWidget *blackSpinButton = createSpinButton(blackPoint + 1, G_CALLBACK(onBlackSpinButtonChanged), greyrampTool);

  createGreyrampToolbar(table, 1, greyrampTool, whiteSpinButton, blackSpinButton, blackPoint, whitePoint);

  GtkWidget *greyrampWindow = createGreyrampToolWindow(dwc, greyrampTool);

  greyrampCreateProperties(greyrampTool,
                           greyrampWindow,
			   dwc,
                           &gradientRect,
                           whiteSpinButton,
                           blackSpinButton,
                           blackPoint,
                           whitePoint,
                           0,
                           255,
                           swapValues,
                           0);

  gtk_spin_button_set_value(GTK_SPIN_BUTTON(whiteSpinButton), (float)whitePoint);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(blackSpinButton), (float)blackPoint);

  gtk_widget_show_all(greyrampTool);

  updateGreyMap(greyrampTool);

  if (greyrampWindow_out)
    *greyrampWindow_out = greyrampWindow;

  DEBUG_EXIT("createGreyrampTool returning ");
  return greyrampTool;
}
