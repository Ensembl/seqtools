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


static void appendChar(char *text1, int *i, char text2)
{
  text1[*i] = text2;
  *i = *i + 1;
}


static void appendText(char *text1, int *i, char *text2)
{
  int j = 0;
  int len = strlen(text2);
  for ( ; j < len; ++j)
    {
      text1[*i] = text2[j];
      *i = *i + 1;
    }
}


//static char* getMarkupText(SequenceCellRenderer *renderer)
//{
//  char *markupText = NULL;
//  
//  MSP *msp = renderer->msp;
//  char *refSeq = getRefSeq(renderer);
//  IntRange *displayRange = getDisplayRange(renderer);
//  int selectedBaseIdx = getSelectedBaseIdx(renderer);
//  
//  if (msp && msp->sseq)
//    {
//      /* Loop through the displayed chars in the reference sequence and compare to the s seq */
//      markupText = g_malloc(sizeof(char) * (displayRange->max - displayRange->min + 50000)); //extra space to allow for markup
//      markupText[0] = '\0';
//      
//      int qIdx = displayRange->min;
//      int i = 0;
//      
//      appendText(markupText, &i, "<span letter_spacing='0'>");
//      
//      for ( ; qIdx < displayRange->max; ++qIdx)
//	{
//	  int sIdx = gapCoord(msp, qIdx, getNumReadingFrames(renderer), getStrand(renderer));
//	  gboolean mismatch = FALSE;
//	  
//	  if (sIdx > 0)
//	    {
//	      if (msp->score < 0)
//		{
//		  /* Ref seq. Draw as it is, just with a yellow background */
//		  if (qIdx == selectedBaseIdx)
//		    appendText(markupText, &i, "<span background='#808000'>");
//		  else
//		    appendText(markupText, &i, "<span background='#ffff00'>");
//		  
//		  appendChar(markupText, &i, msp->sseq[sIdx - 1]);
//		  appendText(markupText, &i, "</span>");
//		}
//	      else if (tolower(msp->sseq[sIdx - 1]) == tolower(refSeq[qIdx - 1]))
//		{
//		  /* Match. Blue background */
//		  if (qIdx == selectedBaseIdx)
//		    appendText(markupText, &i, "<span background='#008080'>");
//		  else
//		    appendText(markupText, &i, "<span background='#00ffff'>");
//		  
//		  appendChar(markupText, &i, msp->sseq[sIdx - 1]);
//		  appendText(markupText, &i, "</span>");
//		}
//	      else
//		{
//		  /* Mismatch. Grey background */
//		  mismatch = TRUE;
//		  if (qIdx == selectedBaseIdx)
//		    appendText(markupText, &i, "<span background='#404040'>");
//		  else
//		    appendText(markupText, &i, "<span background='#bebebe'>");
//		  
//		  appendChar(markupText, &i, msp->sseq[sIdx - 1]);
//		  appendText(markupText, &i, "</span>");
//		}
//	    }
//	  else
//	    {
//	      if (qIdx == selectedBaseIdx)
//		{
//		  appendText(markupText, &i, "<span background='#808080'>");
//		  appendChar(markupText, &i, ' ');
//		  appendText(markupText, &i, "</span>");
//		}
//	      else
//		appendChar(markupText, &i, ' ');
//	    }
//	}
//      
//      appendText(markupText, &i, "</span>");
//      markupText[i] = '\0';
//    }
//  
//  return markupText;
//}


//static PangoLayout* get_layout(GtkWidget *widget, SequenceCellRenderer *renderer)
//{
//  PangoLayout *layout = NULL;
//
//  if (renderer->useMsp)
//    {
//      /* Get the marked-up text of the sequence */
//      char *markupText = getMarkupText(renderer);
//
//      if (markupText)
//	{
//	  /* Parse the markup to create the displayed text and a list of attributes */
//	  char *displayText;
//	  GError *error = NULL;
//	  PangoAttrList *attr_list = pango_attr_list_new ();
//	  
//	  if (markupText && !pango_parse_markup(markupText, -1, 0, &attr_list, &displayText, NULL, &error))
//	    {
//	      g_warning ("Failed to set cell text from markup due to error parsing markup: %s\n\nmarkupText:\n%s\n",
//			 error->message, markupText);
//	      g_error_free (error);
//	    }
//	  else
//	    {
//	      /* Create the layout and set the attributes*/
//	      layout = gtk_widget_create_pango_layout (widget, displayText);
//	      pango_layout_set_attributes (layout, attr_list);
//	      PangoFontDescription *font_desc = pango_font_description_copy(widget->style->font_desc);
//	      pango_layout_set_font_description(layout, font_desc);
//	    }      
//	  
//	  /* clean up */
//	  g_free(displayText);
//	}
//    }
//  else
//    {
//      layout = gtk_widget_create_pango_layout(widget, renderer->text);
//      PangoFontDescription *font_desc = pango_font_description_copy(widget->style->font_desc);
//      pango_layout_set_font_description(layout, font_desc);
//    }
//  
//  
//  return layout;
//}


static void
get_size (GtkCellRenderer *cell,
	  GtkWidget       *widget,
	  GdkRectangle    *cell_area,
	  PangoLayout     *layout,
	  gint            *x_offset,
	  gint            *y_offset,
	  gint            *width,
	  gint            *height)
{
  if (height)
    *height = cell->height;
  
  if (width)
    *width = cell->width;
  
  if (cell_area)
    cell_area->height = cell->height;
}


static void drawText(SequenceCellRenderer *renderer, 
		     GtkWidget *widget,
		     GdkWindow *window, 
		     GtkStateType state, 
		     GdkRectangle *cell_area, 
		     int x_offset, 
		     int y_offset, 
		     int vertical_separator)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(widget, renderer->text);
  PangoFontDescription *font_desc = pango_font_description_copy(widget->style->font_desc);
  pango_layout_set_font_description(layout, font_desc);
  
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


static void appendBaseWithBgColour(char *markupText, MSP *msp, char charToDisplay, int *i, char *bgColour)
{
  char spanStartText[100];
  sprintf(spanStartText, "<span background='#%s'>", bgColour);
  appendText(markupText, i, spanStartText);
  
  appendChar(markupText, i, tolower(charToDisplay)); /* adjust for zero-indexing */
  
  appendText(markupText, i, "</span>");
}


static char* getMarkupTextForBase(MSP *msp, 
				  int qIdx, 
				  char *refSeq, 
				  int selectedBaseIdx, 
				  int numReadingFrames, 
				  Strand strand,
				  int qSeqMin,
				  int qSeqMax)
{
  char *markupText = g_malloc(sizeof(char) * 1000);
  markupText[0] = '\0';
  
  int i = 0;
  appendText(markupText, &i, "<span letter_spacing='0'>");
  
  if (qIdx >= qSeqMin && qIdx <= qSeqMax)
    {
      if (msp->score < 0)
	{
	  /* Exon. Draw a blank space with yellow background */
	  char *bgColour = (qIdx == selectedBaseIdx ? DARK_YELLOW : YELLOW);
	  appendBaseWithBgColour(markupText, msp, ' ', &i, bgColour);
	}
      else if (msp->score == 0 && msp->id == -1) /* (not a great way of identifying the ref seq...) */
	{
	  /* Reference sequence. */
	  char *bgColour = (qIdx == selectedBaseIdx ? DARK_YELLOW : YELLOW);
	  appendBaseWithBgColour(markupText, msp, msp->sseq[qIdx - 1], &i, bgColour);
	}
      else
	{
	  int sIdx = gapCoord(msp, qIdx, numReadingFrames, strand);
	  
	  if (sIdx != UNSET_INT)
	    {
	      /* Find the background colour depending on whether this is a match or not. */
	      char *bgColour;
	      if (tolower(msp->sseq[sIdx - 1]) == tolower(refSeq[qIdx - 1]))
		{
		  /* Match. Blue background */
		  bgColour = (qIdx == selectedBaseIdx ? DARK_BLUE : BLUE);
		}
	      else
		{
		  /* Mismatch. Grey background */
		  bgColour = (qIdx == selectedBaseIdx ? DARK_GREY : GREY);
		}
	      
	      appendBaseWithBgColour(markupText, msp, msp->sseq[sIdx - 1], &i, bgColour);
	    }
	  else
	    {
	      /* There is no equivalent base in the match sequence, i.e. draw a gap. */
	      char *bgColour =  (qIdx == selectedBaseIdx ? DARK_GREY : GREY);
	      appendBaseWithBgColour(markupText, msp, MATCH_SEQ_GAP_CHAR[0], &i, bgColour);
	    }
	}
    }
    else if (qIdx == selectedBaseIdx)
    {
      /* Base does not exist in this match sequence but it is selected, so colour the background */
      appendBaseWithBgColour(markupText, msp, ' ', &i, GREY);
    }
  
  appendText(markupText, &i, "</span>");
  markupText[i] = '\0';
  
  return markupText;
}


static PangoLayout* getLayoutFromMarkup(char *markupText, GtkWidget *widget)
{
  PangoLayout *layout = NULL;
  
  if (markupText)
    {
      /* Parse the markup to create the displayed text and a list of attributes */
      char *displayText;
      GError *error = NULL;
      PangoAttrList *attr_list = pango_attr_list_new ();
      
      if (markupText && !pango_parse_markup(markupText, -1, 0, &attr_list, &displayText, NULL, &error))
	{
	  g_warning ("Failed to set cell text from markup due to error parsing markup: %s\n\nmarkupText:\n%s\n",
		     error->message, markupText);
	  g_error_free (error);
	}
      else
	{
	  /* Create the layout and set the attributes */
	  layout = gtk_widget_create_pango_layout (widget, displayText);
	  pango_layout_set_attributes (layout, attr_list);
	  PangoFontDescription *font_desc = pango_font_description_copy(widget->style->font_desc);
	  pango_layout_set_font_description(layout, font_desc);
	}
      
      g_free(displayText);
      g_free(markupText);
    }
  
  return layout;
}


static void drawBases(SequenceCellRenderer *renderer,
		      GtkWidget *widget,
		      GdkWindow *window, 
		      GtkStateType state, 
		      GdkRectangle *cell_area, 
		      int x_offset, 
		      int y_offset, 
		      int vertical_separator)
{
  /* Extract some info we'll use again and again */
  MSP *msp = renderer->msp;
  char *refSeq = getRefSeq(renderer);
  IntRange *displayRange = getDisplayRange(renderer);
  int selectedBaseIdx = getSelectedBaseIdx(renderer);
  int numReadingFrames = getNumReadingFrames(renderer);
  Strand strand = getStrand(renderer);
//  gboolean qForward = strchr(msp->qframe, '+') ? TRUE : FALSE ;
  
  int qSeqMin, qSeqMax;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, NULL, NULL);
  
  /* Loop through each index in the match sequence that is within the display range. */
  gboolean doneSelectedBase = FALSE;
  int qIdx = qSeqMin >= displayRange->min ? qSeqMin : displayRange->min;
  int lastIdx = qSeqMax <= displayRange->max ? qSeqMax : displayRange->max;

  for ( ; qIdx <= lastIdx; ++qIdx)
    {
      doneSelectedBase = (qIdx == selectedBaseIdx);
      
      char *markupText = getMarkupTextForBase(msp, qIdx, refSeq, selectedBaseIdx, numReadingFrames, strand, qSeqMin, qSeqMax);
      PangoLayout *layout = getLayoutFromMarkup(markupText, widget);
      
      if (layout)
	{
	  /* Calculate the position */
	  int charIdx = qIdx - displayRange->min;
	  int x = cell_area->x + x_offset + GTK_CELL_RENDERER(renderer)->xpad + (charIdx * renderer->charWidth);
	  int y = cell_area->y + y_offset + GTK_CELL_RENDERER(renderer)->ypad - vertical_separator;
	  
	  /* Draw it */
	  gtk_paint_layout (widget->style, window, state, TRUE, NULL, widget, NULL, x, y, layout);
	}
    }
  
  /* Do the selected base index, if not already done, because we always need to colour its background */
  if (selectedBaseIdx != UNSET_INT && !doneSelectedBase)
    {
      char *markupText = getMarkupTextForBase(msp, selectedBaseIdx, refSeq, selectedBaseIdx, numReadingFrames, strand, qSeqMin, qSeqMax);
      PangoLayout *layout = getLayoutFromMarkup(markupText, widget);
      
      if (layout)
	{
	  /* Calculate the position */
	  int charIdx = selectedBaseIdx - displayRange->min;
	  int x = cell_area->x + x_offset + GTK_CELL_RENDERER(renderer)->xpad + (charIdx * renderer->charWidth);
	  int y = cell_area->y + y_offset + GTK_CELL_RENDERER(renderer)->ypad - vertical_separator;
	  
	  /* Draw it */
	  gtk_paint_layout (widget->style, window, state, TRUE, NULL, widget, NULL, x, y, layout);
	}
    }
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
  get_size(cell, widget, cell_area, NULL, x_offset, y_offset, width, height);  
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
  gint x_offset = 0, y_offset = 0, width = 0, height = 0;
  get_size (cell, widget, cell_area, NULL, &x_offset, &y_offset, &width, &height);

  gint vertical_separator;
  gtk_widget_style_get (widget, "vertical-separator", &vertical_separator, NULL);
  height += vertical_separator * 2;

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
  
  
  if (renderer->useMsp)
    drawBases(renderer, widget, window, state, cell_area, x_offset, y_offset, vertical_separator);
  else
    drawText(renderer, widget, window, state, cell_area, x_offset, y_offset, vertical_separator);
  
  
//  PangoLayout *layout = get_layout(widget, renderer);
//  if (layout)
//    {
//      printf("expose area ht = %d; background area ht = %d\n", expose_area->height, background_area->height);
      
      /* 
      // if (renderer->background_set && (flags & GTK_CELL_RENDERER_SELECTED) == 0)
	{
	  GdkColor color;
	  color.red = 0; //renderer->background.red;
	  color.green = 0; //renderer->background.green;
	  color.blue = 0; //renderer->background.blue;
	  
	  GdkGC *gc = gdk_gc_new (window);
	  gdk_gc_set_rgb_fg_color (gc, &color);
	  
	  if (expose_area)               
	    gdk_gc_set_clip_rectangle (gc, expose_area);
	  //background_area->y -= vertical_separator;
	  //background_area->height += vertical_separator;
	  gdk_draw_rectangle (window,
			      gc,
			      TRUE,
			      background_area->x,
			      background_area->y,
			      background_area->width,
			      background_area->height);
	  if (expose_area)               
	    gdk_gc_set_clip_rectangle (gc, NULL);
	  g_object_unref (gc);
	} 
       */
      
//      gtk_paint_layout (widget->style,
//			window,
//			state,
//			TRUE,
//			NULL,
//			widget,
//			NULL, //"cellrenderertext",
//			cell_area->x + x_offset + cell->xpad,
//			cell_area->y + y_offset + cell->ypad - vertical_separator,
//			layout);
//      
//      g_object_unref (layout);
      
      
//	  GdkColor color;
//	  color.red = 0; //renderer->background.red;
//	  color.green = 255; //renderer->background.green;
//	  color.blue = 255; //renderer->background.blue;
//	  
//	  GdkGC *gc = gdk_gc_new (window);
//	  gdk_gc_set_rgb_fg_color (gc, &color);
//	  
//	  if (expose_area)               
//	    gdk_gc_set_clip_rectangle (gc, expose_area);
//	  //background_area->y -= vertical_separator;
//	  //background_area->height += vertical_separator;
//	  gdk_draw_rectangle (window,
//			      gc,
//			      TRUE,
//			      background_area->x,
//			      background_area->y,
//			      background_area->width,
//			      background_area->height);
//	  if (expose_area)               
//	    gdk_gc_set_clip_rectangle (gc, NULL);
//	  g_object_unref (gc);
//    }
}



