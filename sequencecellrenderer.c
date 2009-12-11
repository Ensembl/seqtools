/*
 *  SequenceCellRenderer.c
 *  GtkTest
 *
 *  Created by Gemma Barson on 15/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <SeqTools/sequencecellrenderer.h>
#include <SeqTools/blixem_.h>
#include <SeqTools/detailviewtree.h>

#include <wh/smap.h>

#include <gtk/gtkcellrenderertext.h>


#define MATCH_SEQ_GAP_CHAR	"."	    /* Character to display when there is a gap in the match sequence */
#define YELLOW			"ffff00"
#define DARK_YELLOW		"808000"
#define BLUE			"00ffff"
#define DARK_BLUE		"008080"
#define GREY			"bebebe"
#define DARK_GREY		"404040"

GdkColor yellow = {16776960};
GdkColor dark_yellow = {8421376};
GdkColor cyan = {65535};
GdkColor dark_cyan = {32896};
GdkColor grey = {12500670};
GdkColor dark_grey = {4210752};
GdkColor black = {0};

/* Properties */
enum 
{
  PROP_0,
  
  PROP_DATA,
  PROP_TEXT
};


/* Some boring function declarations: GObject type system stuff */
static void     sequence_cell_renderer_init       (SequenceCellRenderer      *cellrenderersequence);
static void     sequence_cell_renderer_class_init (SequenceCellRendererClass *klass);
static void     sequence_cell_renderer_finalize (GObject *gobject);


/* These functions are the heart of our custom cell renderer: */
static void     sequence_cell_renderer_get_size   (GtkCellRenderer            *cell,
                                                          GtkWidget                  *widget,
                                                          GdkRectangle               *cell_area,
                                                          gint                       *x_offset,
                                                          gint                       *y_offset,
                                                          gint                       *width,
                                                          gint                       *height);
static void     sequence_cell_renderer_render     (GtkCellRenderer            *cell,
                                                          GdkWindow                  *window,
                                                          GtkWidget                  *widget,
                                                          GdkRectangle               *background_area,
                                                          GdkRectangle               *cell_area,
                                                          GdkRectangle               *expose_area,
                                                          guint                       flags);
static void sequence_cell_renderer_get_property (GObject      *object,
					 guint         param_id,
					 GValue	      *value,
					 GParamSpec   *pspec);
static void sequence_cell_renderer_set_property (GObject      *object,
					 guint         param_id,
					 const GValue *value,
					 GParamSpec   *pspec);

static   gpointer parent_class;


/***************************************************************************
 *
 *  sequence_cell_renderer_get_type: here we register our type with
 *                                          the GObject type system if we
 *                                          haven't done so yet. Everything
 *                                          else is done in the callbacks.
 *
 ***************************************************************************/

GType
sequence_cell_renderer_get_type (void)
{
  static GType cell_sequence_type = 0;
  
  if (cell_sequence_type == 0)
    {
      static const GTypeInfo cell_sequence_info =
      {
	sizeof (SequenceCellRendererClass),
	NULL,                                                     /* base_init */
	NULL,                                                     /* base_finalize */
	(GClassInitFunc) sequence_cell_renderer_class_init,
	NULL,                                                     /* class_finalize */
	NULL,                                                     /* class_data */
	sizeof (SequenceCellRenderer),
	0,                                                        /* n_preallocs */
	(GInstanceInitFunc) sequence_cell_renderer_init,
      };
      
      /* Derive from GtkCellRenderer */
      cell_sequence_type = g_type_register_static (GTK_TYPE_CELL_RENDERER,
						   "SequenceCellRenderer",
						   &cell_sequence_info,
						   0);
    }
  
  return cell_sequence_type;
}


/***************************************************************************
 *
 *  sequence_cell_renderer_init: set some default properties of the
 *                                      parent (GtkCellRendererText).
 *
 ***************************************************************************/

static void
sequence_cell_renderer_init (SequenceCellRenderer *cellrenderersequence)
{
  GTK_CELL_RENDERER(cellrenderersequence)->mode = GTK_CELL_RENDERER_MODE_INERT;
  GTK_CELL_RENDERER(cellrenderersequence)->xpad = 2;
  GTK_CELL_RENDERER(cellrenderersequence)->ypad = 0;

  cellrenderersequence->msp = NULL;
  cellrenderersequence->text = NULL;
  cellrenderersequence->charWidth = 0;
  cellrenderersequence->charHeight = 0;
}


/***************************************************************************
 *
 *  sequence_cell_renderer_class_init:
 *
 *  If you want cells that can be activated on their own (ie. not
 *  just the whole row selected) or cells that are editable, you
 *  will need to override 'activate' and 'start_editing'.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_class_init (SequenceCellRendererClass *klass)
{
  GtkCellRendererClass *cell_class   = GTK_CELL_RENDERER_CLASS(klass);
  GObjectClass         *object_class = G_OBJECT_CLASS(klass);
  
  parent_class           = g_type_class_peek_parent (klass);
  object_class->finalize = sequence_cell_renderer_finalize;

  /* Override the two crucial functions that are the heart
   *   of a cell renderer in the parent class */
  cell_class->get_size = sequence_cell_renderer_get_size;
  cell_class->render   = sequence_cell_renderer_render;
  
  /* Override functions to get/set properties */
  object_class->get_property = sequence_cell_renderer_get_property;
  object_class->set_property = sequence_cell_renderer_set_property;
  
  g_object_class_install_property (object_class,
                                   PROP_DATA,
                                   g_param_spec_pointer ("msp",
                                                         "MSP",
                                                         "Pointer to an msp whose sequence to display",
                                                         G_PARAM_WRITABLE));
  
  g_object_class_install_property (object_class,
                                   PROP_TEXT,
                                   g_param_spec_string ("text",
                                                        "Text",
                                                        "Text to display",
							NULL,
                                                        G_PARAM_WRITABLE));
  
}


/***************************************************************************
 *
 *  sequence_cell_renderer_finalize: free any resources here
 *
 ***************************************************************************/

static void
sequence_cell_renderer_finalize (GObject *object)
{
  /*
   SequenceCellRenderer *cellrenderersequence = SEQUENCE_CELL_RENDERER(object);
   */
  
  /* Free any dynamically allocated resources here */
  
  (* G_OBJECT_CLASS (parent_class)->finalize) (object);
}


/***************************************************************************
 *
 *  sequence_cell_renderer_new: return a new cell renderer instance
 *
 ***************************************************************************/

GtkCellRenderer *
sequence_cell_renderer_new (void)
{
  return g_object_new(SEQUENCE_CELL_RENDERER_TYPE, NULL);
}


/***************************************************************************
 *
 *  sequence_cell_renderer_set_property: 
 *
 ***************************************************************************/

static void
sequence_cell_renderer_set_property (GObject      *object,
				     guint         param_id,
				     const GValue *value,
				     GParamSpec   *pspec)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(object);
  
  switch (param_id)
  {
    case PROP_DATA:
      renderer->msp = (MSP*)g_value_get_pointer(value);
      renderer->useMsp = TRUE;
      break;
      
    case PROP_TEXT:
      renderer->text = g_strdup(g_value_get_string(value));
      renderer->useMsp = FALSE;
      break;
      
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
      break;
  };
}


/***************************************************************************
 *
 *  sequence_cell_renderer_get_property: 
 *
 ***************************************************************************/

static void
sequence_cell_renderer_get_property (GObject      *object,
					 guint         param_id,
					 GValue	      *value,
					 GParamSpec   *pspec)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(object);
  switch (param_id)
  {
    case PROP_DATA:
      g_value_set_pointer(value, renderer->msp);
      break;
      
    case PROP_TEXT:
      g_value_set_string(value, renderer->text);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
      break;
  };
}


/* Returns the current strand that this renderer is displaying */
static int getStrand(SequenceCellRenderer *renderer)
{
  return treeGetStrand(renderer->tree);
}


/* Returns the reference sequence */
static char* getRefSeq(SequenceCellRenderer *renderer)
{
  return treeGetRefSeq(renderer->tree);
}


/* Get the detail view display range */
static IntRange* getDisplayRange(SequenceCellRenderer *renderer)
{
  return treeGetDisplayRange(renderer->tree);
}


/* Get the main window strands-toggled status */
static gboolean getStrandsToggled(SequenceCellRenderer *renderer)
{
  return treeGetStrandsToggled(renderer->tree);
}


/* Get the number of reading frames in the detail view */
static int getNumReadingFrames(SequenceCellRenderer *renderer)
{
  return treeGetNumReadingFrames(renderer->tree);
}


/* Get the currently-selected base in the reference sequence */
static int getSelectedBaseIdx(SequenceCellRenderer *renderer)
{
  return treeGetSelectedBaseIdx(renderer->tree);
}


/* input is co-ord on query sequence, find corresonding base in subject sequence */
int gapCoord(const MSP *msp, int qIdx, int numFrames, Strand strand)
{
  int result = UNSET_INT ;
  
  int qSeqMin, qSeqMax, sSeqMin, sSeqMax;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
  
  BOOL qForward = ((strchr(msp->qframe, '+'))) ? TRUE : FALSE ;
  BOOL sForward = ((strchr(msp->sframe, '+'))) ? TRUE : FALSE ;
  BOOL sameDirection = (qForward == sForward);
  
  //BOOL rightToLeft = (qForward && plusmin < 0) || (!qForward && plusmin > 0); 
  
  Array gaps = msp->gaps ;
  if (!gaps || arrayMax(gaps) < 1)
    {
      /* If strands are in the same direction, find the offset from qSeqMin and add it to sSeqMin.
       * If strands are in opposite directions, find the offset from qSeqMin and subtract it from sSeqMax. */
      int offset = (qIdx - qSeqMin)/numFrames ;
      result = (sameDirection) ? sSeqMin + offset : sSeqMax - offset ;
    }
  else
    {
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
		}
  //		  else if (i == 0)
  //		    {
  //		      /* qIdx lies before the first range. Use the start of the first range. */
  //		      result = curRange->s1;
  //		      
  //		      if (!rightToLeft)
  //			{
  //			  /* This is a special case where for normal left-to-right display we want to display the
  //			   * gap before the base rather than after it.  Use 1 beyond the edge of the range to do this. */
  //			  result = sForward ? result - 1 : result + 1 ;
  //			}
  //		    }
	      
	      break;
	    }
  //	      else
  //		{
  //		  /* Remember the end of the current range (which is the result we want if qIdx lies after the last range). */
  //		  result = curRange->s2;
  //		  
  //		  if (rightToLeft)
  //		    {
  //		      /* For right-to-left display, we want to display the gap before the base rather 
  //		       * than after it.  Use 1 index beyond the edge of the range to do this. */
  //		      result = sForward ? result + 1 : result - 1 ;
  //		    }
  //		}
	}
    }
  
  return result ;
}


static void
get_size (GtkCellRenderer *cell,
	  GtkWidget       *widget,
	  GdkRectangle    *cell_area,
	  PangoLayout     *layout,
	  gint            *x_offset,
	  gint            *y_offset,
	  gint            *width,
	  gint            *height,
	  gint		  *vertical_separator)
{
  if (height)
    *height = cell->height;
  
  if (width)
    *width = cell->width;
  
  if (cell_area)
    cell_area->height = cell->height;
  
  if (vertical_separator)
    {
      gtk_widget_style_get (widget, "vertical-separator", vertical_separator, NULL);
    }
}


static void drawText(SequenceCellRenderer *renderer, 
		     GtkWidget *widget,
		     GdkWindow *window, 
		     GtkStateType state, 
		     GdkRectangle *cell_area)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, renderer->text);
  PangoFontDescription *font_desc = pango_font_description_copy(widget->style->font_desc);
  pango_layout_set_font_description(layout, font_desc);

  gint x_offset = 0, y_offset = 0, vertical_separator = 0;
  get_size (GTK_CELL_RENDERER(renderer), widget, cell_area, NULL, &x_offset, &y_offset, NULL, NULL, &vertical_separator);

  gtk_paint_layout (widget->style,
		    window,
		    state,
		    TRUE,
		    NULL,
		    widget,
		    NULL, //"cellrenderertext",
		    cell_area->x + x_offset + GTK_CELL_RENDERER(renderer)->xpad,
		    cell_area->y + y_offset + GTK_CELL_RENDERER(renderer)->ypad - vertical_separator,
		    layout);
}


static char getRefSeqBase(char *refSeq, int qIdx, Strand qStrand)
{
  char result;
  char base = tolower(refSeq[qIdx - 1]);
  
  if (qStrand == FORWARD_STRAND)
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
  
  return result;
}


/* Place the given character into the given string at the given index, then increment the index.
 * The given string may be null, in which case this function does nothing. */
static void insertChar(char *text1, int *i, char charToAdd, MSP *msp)
{
  if (text1)
    {
      text1[*i] = charToAdd;
      *i = *i + 1;
    }
}


/* Colour in the background of a particular base in the given match sequence. Returns
 * the calculated index into the match sequence (for effiency, so that we don't have to
 * recalculate it later on). */
static void drawBaseBackgroundColour(MSP *msp, 
				     int qIdx, 
				     char *refSeq, 
				     int selectedBaseIdx, 
				     int numReadingFrames, 
				     Strand qStrand,
				     int qSeqMin,
				     int qSeqMax,
				     int x, 
				     int y,
				     int width,
				     int height,
				     GdkWindow *window,
				     GdkGC *gc,
				     char *displayText,
				     int *strIdx)
{
  gboolean selected = (qIdx == selectedBaseIdx);
  
  char charToDisplay = 'x';
  GdkColor *baseBgColour;
  
  if (qIdx >= qSeqMin && qIdx <= qSeqMax)
    {
      if (msp->score < 0)
	{
	  /* Exon. Draw a blank space with yellow background */
	  charToDisplay = ' ';
	  baseBgColour = selected ? &dark_yellow : &yellow;
	}
      else if (msp->score == 0 && msp->id == -1) /* (not a great way of identifying the ref seq...) */
	{
	  /* Reference sequence. */
	  charToDisplay = getRefSeqBase(refSeq, qIdx, qStrand);
	  baseBgColour = selected ? &dark_yellow : &yellow;
	}
      else
	{
	  int sIdx = gapCoord(msp, qIdx, numReadingFrames, qStrand);
	  
	  if (sIdx != UNSET_INT)
	    {
	      /* Find the background colour depending on whether this is a match or not. */
	      charToDisplay = tolower(msp->sseq[sIdx - 1]);
	      char qBase = getRefSeqBase(refSeq, qIdx, qStrand);
	      
	      if (charToDisplay == qBase)
		{
		  /* Match. Blue background */
		  baseBgColour = selected ? &dark_cyan : &cyan;
		}
	      else
		{
		  /* Mismatch. Grey background */
		  baseBgColour = selected ? &dark_grey : &grey;
		}
	    }
	  else
	    {
	      /* There is no equivalent base in the match sequence so draw a gap */
	      charToDisplay = '.';
	      baseBgColour = selected ? &dark_grey : &grey;
	    }
	}
    }
    else if (selected)
    {
      /* Base does not exist in this match sequence but it is selected, so colour the background */
      charToDisplay = ' ';
      baseBgColour = &grey;
    }

  if (charToDisplay != 'x')
    {
      insertChar(displayText, strIdx, charToDisplay, msp);
      
      if (GDK_IS_GC(gc))
	{
	  gdk_gc_set_foreground(gc, baseBgColour);
	  gdk_draw_rectangle(window, gc, TRUE, x, y, width, height);
	}
      else
	{
	  messerror("Invalid graphics context");
	}
    }  
}


static PangoLayout* getLayoutFromText(char *displayText, GtkWidget *widget, PangoFontDescription *font_desc)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, displayText);
  pango_layout_set_font_description(layout, font_desc);
  return layout;
}


/* Return the x/y coords for the top-left corner where we want to draw the base
 * with the given index in the reference sequence */
static void getCoordsForBaseIdx(int qIdx, 
				SequenceCellRenderer *renderer, 
				IntRange *displayRange, 
				gboolean rightToLeft,
				GdkRectangle *cell_area, 
				int x_offset, 
				int y_offset, 
				int vertical_separator, 
				int *x, 
				int* y)
{
  /* Calculate the index of the character within the cell. */
  int charIdx = UNSET_INT;
  if (rightToLeft)
    {
      /* Strands have been toggled, so display sequences reversed (i.e. so they read right-to-left) */
      charIdx = displayRange->max - qIdx - 1;
    }
  else
    {
      /* Normal left-to-right display */
      charIdx = qIdx - displayRange->min;
    }

  /* Calculate the coords */
  *x = cell_area->x + x_offset + GTK_CELL_RENDERER(renderer)->xpad + (charIdx * renderer->charWidth);
  *y = cell_area->y + y_offset + GTK_CELL_RENDERER(renderer)->ypad - vertical_separator;
}


static void drawBases(SequenceCellRenderer *renderer,
		      GtkWidget *widget,
		      GdkWindow *window, 
		      GtkStateType state,
		      GdkRectangle *cell_area)
{
  gint x_offset = 0, y_offset = 0, vertical_separator = 0;
  get_size (GTK_CELL_RENDERER(renderer), widget, cell_area, NULL, &x_offset, &y_offset, NULL, NULL, &vertical_separator);
  
  PangoFontDescription *font_desc = pango_font_description_copy(widget->style->font_desc);
  GdkGC *gc = gdk_gc_new(window);
  
  IntRange *displayRange = getDisplayRange(renderer);
  gboolean rightToLeft = getStrandsToggled(renderer);
  int selectedBaseIdx = getSelectedBaseIdx(renderer);
  
  int qSeqMin, qSeqMax;
  getMspRangeExtents(renderer->msp, &qSeqMin, &qSeqMax, NULL, NULL);
  
  /* We'll through each index in the match sequence that is within the display range. 
   * First, find the min and max values to display. */
  int qMin = qSeqMin >= displayRange->min ? qSeqMin : displayRange->min;
  int qMax = qSeqMax <= displayRange->max ? qSeqMax : displayRange->max;

  /* We'll start drawing from the left. In normal left-to-right display, low values
   * are on the left and high values on the right. If the display is toggled, this
   * is reversed. */
  int qStart = rightToLeft ? qMax : qMin;
  
  /* We'll populate a string with the characters to display as we loop through the indices. */
  char displayText[qMax - qMin + 2];
  int strIdx = 0;
  
  gboolean doneSelectedBase = FALSE;

  /* Ok, let's loop through the bases we're interested in displaying. */
  int qIdx;
  for (qIdx = qStart; qIdx >= qMin && qIdx <= qMax; (rightToLeft ? --qIdx : ++qIdx))
    {
      int x, y;
      getCoordsForBaseIdx(qIdx, renderer, displayRange, rightToLeft, cell_area, x_offset, y_offset, vertical_separator, &x, &y);
      
      drawBaseBackgroundColour(renderer->msp, 
				qIdx, 
				getRefSeq(renderer), 
				selectedBaseIdx,
				getNumReadingFrames(renderer), 
				getStrand(renderer), 
				qSeqMin, 
				qSeqMax, 
				x, 
				y, 
				renderer->charWidth, 
				renderer->charHeight,
				window, 
				gc,
				displayText,
				&strIdx);
      
      doneSelectedBase = doneSelectedBase || (qIdx == selectedBaseIdx);
    }

  /* Null-terminate the string */
  insertChar(displayText, &strIdx, '\0', renderer->msp);


  /* Make sure we always colour the selected base index's background for every row */
  if (!doneSelectedBase && selectedBaseIdx >= displayRange->min && selectedBaseIdx <= displayRange->max)
    {
      int x, y;
      getCoordsForBaseIdx(selectedBaseIdx, renderer, displayRange, rightToLeft, cell_area, x_offset, y_offset, vertical_separator, &x, &y);

      drawBaseBackgroundColour(renderer->msp, 
			       selectedBaseIdx, 
			       getRefSeq(renderer), 
			       selectedBaseIdx, 
			       getNumReadingFrames(renderer), 
			       getStrand(renderer), 
			       qSeqMin, 
			       qSeqMax, 
			       x, 
			       y, 
			       renderer->charWidth, 
			       renderer->charHeight,
			       window, 
			       gc,
			       NULL,
			       NULL);
    }
  
  /* Draw the sequence text over the top of the coloured backgrounds. */
  /* Get the coords for the first base we're displaying. The display text was 
   * constructed such that everything else will line up from here. */
  if (g_utf8_validate(displayText, -1, NULL))
    {
      int x, y;
      getCoordsForBaseIdx(qStart, renderer, displayRange, rightToLeft, cell_area, x_offset, y_offset, vertical_separator, &x, &y);

      PangoLayout *layout = getLayoutFromText(displayText, widget, font_desc);
      gtk_paint_layout (widget->style, window, state, TRUE, NULL, widget, NULL, x, y, layout);
    }
  else
    {
      messerror("Invalid string constructed when trying to display sequence:\n%s\n", displayText);
    }
}    


static GtkStateType getState(GtkWidget *widget, guint flags)
{
  GtkStateType state;
  if ((flags & GTK_CELL_RENDERER_SELECTED) == GTK_CELL_RENDERER_SELECTED)
    {
      state = GTK_WIDGET_HAS_FOCUS(widget) ? GTK_STATE_SELECTED : GTK_STATE_ACTIVE;
    }
  else if ((flags & GTK_CELL_RENDERER_PRELIT) == GTK_CELL_RENDERER_PRELIT && GTK_WIDGET_STATE (widget) == GTK_STATE_PRELIGHT)
    {
      state = GTK_STATE_PRELIGHT;
    }
  else if (GTK_WIDGET_STATE (widget) == GTK_STATE_INSENSITIVE)
    {
      state = GTK_STATE_INSENSITIVE;
    }
  else
    {
      state = GTK_STATE_NORMAL;
    }
  
  return state;
}

/***************************************************************************
 *
 *  sequence_cell_renderer_get_size: crucial - calculate the size
 *                                          of our cell, taking into account
 *                                          padding and alignment properties
 *                                          of parent.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_get_size (GtkCellRenderer *cell,
                                        GtkWidget       *widget,
                                        GdkRectangle    *cell_area,
                                        gint            *x_offset,
                                        gint            *y_offset,
                                        gint            *width,
                                        gint            *height)
{
  get_size(cell, widget, cell_area, NULL, x_offset, y_offset, width, height, NULL);  
}


/***************************************************************************
 *
 *  sequence_cell_renderer_render: crucial - do the rendering.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_render (GtkCellRenderer *cell,
			       GdkWindow       *window,
			       GtkWidget       *widget,
			       GdkRectangle    *background_area,
			       GdkRectangle    *cell_area,
			       GdkRectangle    *expose_area,
			       guint            flags)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  if (renderer->useMsp)
    {
      drawBases(renderer, 
		widget, 
		window, 
		getState(widget, flags),
		cell_area);
    }
  else
    {
      drawText(renderer,
	       widget, 
	       window, 
	       getState(widget, flags),
	       cell_area);
    }
}



