/*
 *  SequenceCellRenderer.c
 *  GtkTest
 *
 *  Created by Gemma Barson on 15/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <SeqTools/sequencecellrenderer.h>
#include <SeqTools/detailview.h>
#include <SeqTools/detailviewtree.h>
#include <SeqTools/blxviewMainWindow.h>
#include <SeqTools/utilities.h>
#include <wh/smap.h>
#include <gtk/gtkcellrenderertext.h>


#define SEQUENCE_CELL_RENDERER_NAME	"SequenceCellRenderer"
#define GAP_WIDTH_AS_FRACTION		0.375	/* multiplier used to get the width of the "gap" marker based on a fraction of char width */
#define MIN_GAP_WIDTH			2

typedef struct _RenderData
  {
    GdkRectangle *cell_area;
    GdkWindow *window;
    GtkStateType state;
    GdkGC *gc;
    const gboolean rightToLeft;
    const Strand qStrand;
    const int qFrame;
    const int selectedBaseIdx;
    const int cellXPadding;
    const int cellYPadding;
    const int numFrames;
    const int charWidth;
    const int charHeight;
    const BlxSeqType seqType;
    const BlxBlastMode blastMode;
    const IntRange const *displayRange;
    const IntRange const *refSeqRange;
    const gboolean highlightDiffs;
    GdkDrawable *drawable;
    GtkWidget *mainWindow;
    GdkColor *exonColour;
    GdkColor *exonColourSelected;
    GdkColor *insertionColour;
    GdkColor *insertionColourSelected;
    GdkColor *matchColour;
    GdkColor *matchColourSelected;
    GdkColor *consColour;
    GdkColor *consColourSelected;
    GdkColor *mismatchColour;
    GdkColor *mismatchColourSelected;
    GdkColor *exonBoundaryColourStart;
    GdkColor *exonBoundaryColourEnd;
    int exonBoundaryWidth;
    GdkLineStyle exonBoundaryStyleStart;
    GdkLineStyle exonBoundaryStyleEnd;
  } RenderData;



/* Properties */
enum 
{
  PROP_0,
  
  PROP_NAME,
  PROP_SCORE,
  PROP_ID,
  PROP_START,
  PROP_MSP,
  PROP_END,
  PROP_DATA
};


/* Some boring function declarations: GObject type system stuff */
static void     sequence_cell_renderer_init       (SequenceCellRenderer      *cellrenderersequence);
static void     sequence_cell_renderer_class_init (SequenceCellRendererClass *klass);
static void     sequence_cell_renderer_finalize (GObject *gobject);

static IntRange getVisibleMspRange(MSP *msp, RenderData *data);
void		drawVisibleExonBoundaries(GtkWidget *tree, RenderData *data);

static void drawSequenceText(GtkWidget *tree,
			     gchar *displayText, 
			     const IntRange const *segmentRange,
			     RenderData *data);


static void getCoordsForBaseIdx(const int segmentIdx, 
				const IntRange const *segmentRange,
				RenderData *data,
				int *x, 
				int* y);


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
                                                          GtkWidget                  *tree,
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



/***************************************************************************/

/* Custom version of gtk_paint_layout that draws the same thing to two GdkDrawables */
static void paintLayout2(GtkStyle *style,
			 GdkDrawable *drawable1,
			 GdkDrawable *drawable2,
			 GtkStateType state_type,
			 gboolean use_text,
			 GdkRectangle *area,
			 GtkWidget *widget,
			 gchar *detail,
			 gint x,
			 gint y,
			 PangoLayout *layout)
{
  if (drawable1)
    gtk_paint_layout (style, drawable1, state_type, use_text, area, widget, detail, x, y, layout);
  
  if (drawable2)
    gtk_paint_layout (style, drawable2, state_type, use_text, area, widget, detail, x, y, layout);
}

/* Custom version of gdk_draw_rectangle that draws the same thing to two GdkDrawables */
void drawRectangle2(GdkDrawable *drawable1,
		    GdkDrawable *drawable2,
		    GdkGC *gc,
		    gboolean filled,
		    gint x,
		    gint y,
		    gint width,
		    gint height)
{
  if (drawable1)
    gdk_draw_rectangle(drawable1, gc, filled, x, y, width, height);
  
  if (drawable2)
    gdk_draw_rectangle(drawable2, gc, filled, x, y, width, height);
}

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
  GTK_CELL_RENDERER(cellrenderersequence)->xalign = 0;
  GTK_CELL_RENDERER(cellrenderersequence)->yalign = 0;
  GTK_CELL_RENDERER(cellrenderersequence)->xpad = 0;
  GTK_CELL_RENDERER(cellrenderersequence)->ypad = 0;

  cellrenderersequence->data = NULL;
  cellrenderersequence->mspGList = NULL;
  cellrenderersequence->name = NULL;
  cellrenderersequence->score = NULL;
  cellrenderersequence->id = NULL;
  cellrenderersequence->start = NULL;
  cellrenderersequence->end = NULL;
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
                                   g_param_spec_pointer ("data",
                                                         "Data",
                                                         "Pointer to the msp",
                                                         G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_MSP,
                                   g_param_spec_pointer (SEQ_COLUMN_PROPERTY_NAME,
                                                         SEQ_COLUMN_HEADER_TEXT,
                                                         "Pointer to an msp whose sequence to display",
                                                         G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_NAME,
                                   g_param_spec_string (NAME_COLUMN_PROPERTY_NAME,
                                                        NAME_COLUMN_HEADER_TEXT,
                                                        "Sequence name",
							NULL,
                                                        G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_SCORE,
                                   g_param_spec_string (SCORE_COLUMN_PROPERTY_NAME,
                                                        SCORE_COLUMN_HEADER_TEXT,
                                                        "Score",
							NULL,
                                                        G_PARAM_WRITABLE));
							
  g_object_class_install_property (object_class,
                                   PROP_ID,
                                   g_param_spec_string (ID_COLUMN_PROPERTY_NAME,
                                                        ID_COLUMN_HEADER_TEXT,
                                                        "Identity",
							NULL,
                                                        G_PARAM_WRITABLE));
							
  g_object_class_install_property (object_class,
                                   PROP_START,
                                   g_param_spec_string (START_COLUMN_PROPERTY_NAME,
                                                        START_COLUMN_HEADER_TEXT,
                                                        "Start index of match",
							NULL,
                                                        G_PARAM_WRITABLE));
							
  g_object_class_install_property (object_class,
                                   PROP_END,
                                   g_param_spec_string (END_COLUMN_PROPERTY_NAME,
                                                        END_COLUMN_HEADER_TEXT,
                                                        "End index of match",
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

static void setAllPropertiesNull(SequenceCellRenderer *renderer)
{
  renderer->name = renderer->score = renderer->id = renderer->start = renderer->end = NULL;
  renderer->mspGList = NULL;
}

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
      /* Additional data. DON'T reset other properties. */
      renderer->data = (GList*)g_value_get_pointer(value);
      break;
      
    case PROP_MSP:
      setAllPropertiesNull(renderer);
      renderer->mspGList = (GList*)g_value_get_pointer(value);
      break;

    case PROP_NAME:
      setAllPropertiesNull(renderer);
      renderer->name = g_strdup(g_value_get_string(value));
      break;
      
    case PROP_SCORE:
      setAllPropertiesNull(renderer);
      renderer->score = g_strdup(g_value_get_string(value));
      break;

    case PROP_ID:
      setAllPropertiesNull(renderer);
      renderer->id = g_strdup(g_value_get_string(value));
      break;

    case PROP_START:
      setAllPropertiesNull(renderer);
      renderer->start = g_strdup(g_value_get_string(value));
      break;

    case PROP_END:
      setAllPropertiesNull(renderer);
      renderer->end = g_strdup(g_value_get_string(value));
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
      g_value_set_pointer(value, renderer->data);
      break;

    case PROP_MSP:
      g_value_set_pointer(value, renderer->mspGList);
      break;
      
    case PROP_NAME:
      g_value_set_string(value, renderer->name);
      break;

    case PROP_SCORE:
      g_value_set_string(value, renderer->score);
      break;

    case PROP_ID:
      g_value_set_string(value, renderer->id);
      break;

    case PROP_START:
      g_value_set_string(value, renderer->start);
      break;

    case PROP_END:
      g_value_set_string(value, renderer->end);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, param_id, pspec);
      break;
  };
}


static void
get_size (GtkCellRenderer *cell,
	  GtkWidget       *widget,
	  GdkRectangle    *cell_area,
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


static char* getText(SequenceCellRenderer *renderer)
{
  if (renderer->name)
    return renderer->name;
  else if (renderer->score)
    return renderer->score;
  else if (renderer->id)
    return renderer->id;
  else if (renderer->start)
    return renderer->start;
  else if (renderer->end)
    return renderer->end;
    
  return NULL;
}


/* Render function for a cell that contains simple text */
static void drawText(SequenceCellRenderer *renderer, 
		     GtkWidget *tree,
		     GdkWindow *window, 
		     GtkStateType state, 
		     GdkRectangle *cell_area)
{
  gchar *displayText = getText(renderer);
  PangoLayout *layout = gtk_widget_create_pango_layout(tree, displayText);

  PangoFontDescription *font_desc = pango_font_description_copy(tree->style->font_desc);
  pango_layout_set_font_description(layout, font_desc);
  pango_font_description_free(font_desc);

  paintLayout2(tree->style,
	       window,
	       widgetGetDrawable(tree),
	       state,
	       TRUE,
	       NULL,
	       tree,
	       NULL, 
	       cell_area->x - treeGetCellXPadding(tree),
	       cell_area->y - treeGetCellYPadding(tree),
	       layout);
  
  g_object_unref(layout);
}


/* The given renderer is an MSP. This function checks if there is a base index
 * selected and, if so, colours the background for that base with the given colour. */
static void highlightSelectedBase(const int selectedBaseIdx,
				  GdkColor *highlightColour,
				  RenderData *data)
{
  if (selectedBaseIdx != UNSET_INT && valueWithinRange(selectedBaseIdx, data->displayRange))
    {
      /* Convert the display-range index to a 0-based index for the section of sequence displayed */
      const int segmentIdx = selectedBaseIdx - data->displayRange->min;
      
      int x, y;
      getCoordsForBaseIdx(segmentIdx, data->displayRange, data, &x, &y);

      gdk_gc_set_foreground(data->gc, highlightColour);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, data->charWidth, data->charHeight);
    }
}


/* The given renderer is an MSP that is an exon. This function draws the exon
 * or part of the exon that is in view, if it is within the current display range. */
static void drawExon(SequenceCellRenderer *renderer,
		     MSP *msp,
		     GtkWidget *tree,
		     RenderData *data)
{
  IntRange segmentRange = getVisibleMspRange(msp, data);
  const int segmentLen = segmentRange.max - segmentRange.min + 1;

  int x, y;
  getCoordsForBaseIdx(0, &segmentRange, data, &x, &y);
  const int width = segmentLen * data->charWidth;

  /* Just draw one big rectangle the same colour for the whole thing */
  gdk_gc_set_foreground(data->gc, data->exonColour);
  drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, width, data->charHeight);
  
  /* If a base is selected, highlight it. Its colour depends on whether it the base is within the exon range or not. */
  if (data->selectedBaseIdx != UNSET_INT && valueWithinRange(data->selectedBaseIdx, &segmentRange))
    {
      highlightSelectedBase(data->selectedBaseIdx, data->exonColourSelected, data);
    }
  
  drawVisibleExonBoundaries(tree, data);
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


/* Get the base from the given match sequence at the given index. Converts to
 * upper or lower case as appropriate for the sequence type. The given index is
 * assumed to be 1-based. */
static char getMatchSeqBase(char *matchSeq, const int sIdx, const BlxSeqType seqType)
{
  char result = '-';
  
  if (matchSeq)
  {
    result = matchSeq[sIdx - 1];
    result = convertBaseToCorrectCase(result, seqType);
  }

  return result;
}


/* Given a 0-based index into a segment of the reference sequence, find the
 * index into the full reference sequence. */
static int getRefSeqIndexFromSegment(const int segmentIdx,
				     const IntRange const *segmentRange)
{
  return (segmentRange->min + segmentIdx);
}


/* Colour in the background of a particular base in the given match sequence. Returns
 * the calculated index into the match sequence (for effiency, so that we don't have to
 * recalculate it later on). Returns the equivalent index in the subject sequence, or 
 * UNSET_INT if there is none. */
static void drawBase(MSP *msp,
		     const int segmentIdx, 
		     const IntRange const *segmentRange,
		     char *refSeqSegment, 
		     RenderData *data,
		     const int x,
		     const int y,
		     gchar *displayText,
		     int *sIdx,
		     int *qIdx)
{
  char sBase = '\0';
  GdkColor *baseBgColour;
  
  *qIdx = getRefSeqIndexFromSegment(segmentIdx, segmentRange);
  *sIdx = getMatchIdxFromDisplayIdx(msp, *qIdx, data->qFrame, data->qStrand, data->rightToLeft, data->seqType, data->numFrames, data->refSeqRange);
  
  /* Highlight the base if its base index is selected */
  gboolean selected = (*qIdx == data->selectedBaseIdx);
  
  if (*sIdx == UNSET_INT)
    {
      /* There is no equivalent base in the match sequence so draw a gap */
      sBase = '.';
      baseBgColour = selected ? data->mismatchColourSelected : data->mismatchColour;
    }
  else
    {
      /* There is a base in the match sequence. See if it matches the ref sequence */
      sBase = getMatchSeqBase(msp->sseq, *sIdx, data->seqType);
      char qBase = refSeqSegment[segmentIdx];

      if (sBase == qBase)
	{
	  /* Match */
	  baseBgColour = selected ? data->matchColourSelected : data->matchColour;
	  
	  /* If we're highlighting differences, don't show this base (show a dash instead) */
	  if (data->highlightDiffs)
	    {
	      sBase = '-';
	    }
	}
      else if (data->blastMode != BLXMODE_BLASTN && PAM120[aa_atob[(unsigned int)qBase]-1 ][aa_atob[(unsigned int)sBase]-1 ] > 0)
	{
	  /* 'Conserved' match (i.e. similar amino acid) */
	  baseBgColour = selected ? data->consColourSelected : data->consColour;
	}
      else
	{
	  /* Mismatch */
	  baseBgColour = selected ? data->mismatchColourSelected : data->mismatchColour;
	}
    }

  /* Draw the background colour */
  gdk_gc_set_foreground(data->gc, baseBgColour);
  drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, data->charWidth, data->charHeight);

  /* Add this character into the display text */
  displayText[segmentIdx] = sBase;
}


static PangoLayout* getLayoutFromText(gchar *displayText, GtkWidget *tree, PangoFontDescription *font_desc)
{
  PangoLayout *layout = gtk_widget_create_pango_layout(tree, displayText);
  pango_layout_set_font_description(layout, font_desc);
  return layout;
}


/* Return the x/y coords for the top-left corner where we want to draw the base
 * with the given index in the segment, where the segment starts at the given 
 * index in the display. */
static void getCoordsForBaseIdx(const int segmentIdx, 
				const IntRange const *segmentRange,
				RenderData *data,
				int *x, 
				int* y)
{
  /* Find the start of the segment with respect to the display range */
  const int startPos = segmentRange->min - data->displayRange->min;
  
  /* Find the position of the character within the segment */
  int charIdx = startPos + segmentIdx;

  /* Calculate the coords */
  *x = data->cell_area->x - data->cellXPadding + (charIdx * data->charWidth);
  *y = data->cell_area->y - data->cellYPadding;
}


/* Draw the start/end boundaries for the given MSP if it is an exon whose start/end
 * coords are within the current display range */
static gboolean drawExonBoundary(const MSP *msp, RenderData *rd)
{
  if (msp && mspIsExon(msp))
    {
      const int minIdx = min(msp->displayStart, msp->displayEnd);
      const int maxIdx = max(msp->displayStart, msp->displayEnd);
      
      if (minIdx >= rd->displayRange->min && minIdx <= rd->displayRange->max)
	{
	  /* Draw the lower index. The colour and line style depend on whether it's the start or end index. */
	  gdk_gc_set_foreground(rd->gc, rd->exonBoundaryColourStart);
	  gdk_gc_set_line_attributes(rd->gc, rd->exonBoundaryWidth, rd->exonBoundaryStyleStart, GDK_CAP_BUTT, GDK_JOIN_MITER);

	  const int idx = rd->rightToLeft ? rd->displayRange->max - minIdx + 1 : minIdx - rd->displayRange->min;

	  int x, y;
	  getCoordsForBaseIdx(idx, rd->displayRange, rd, &x, &y);
	  
	  gdk_draw_line(rd->window, rd->gc, x, y, x, y + rd->charHeight);
	  gdk_draw_line(rd->drawable, rd->gc, x, y, x, y + rd->charHeight);
	}
      
      if (maxIdx >= rd->displayRange->min && maxIdx <= rd->displayRange->max)
	{
	  /* Draw the upper index. The colour and line style depend on whether it's the start or end index. */
	  gdk_gc_set_foreground(rd->gc, rd->exonBoundaryColourEnd);
	  gdk_gc_set_line_attributes(rd->gc, rd->exonBoundaryWidth, rd->exonBoundaryStyleEnd, GDK_CAP_BUTT, GDK_JOIN_MITER);
	  
	  const int idx = rd->rightToLeft ? rd->displayRange->max - maxIdx : maxIdx + 1 - rd->displayRange->min;

	  int x, y;
	  getCoordsForBaseIdx(idx, rd->displayRange, rd, &x, &y);

	  gdk_draw_line(rd->window, rd->gc, x, y, x, y + rd->charHeight);
	  gdk_draw_line(rd->drawable, rd->gc, x, y, x, y + rd->charHeight);
	}
    }
  
  return FALSE;
}


/* Draw the boundaries of all exons in the given tree that are within the current
 * display range */
void drawVisibleExonBoundaries(GtkWidget *tree, RenderData *data)
{
  /* Loop through all MSPs. */
  const MSP *msp = mainWindowGetMspList(data->mainWindow);

  for ( ; msp; msp = msp->next)
    {
      if (mspGetRefFrame(msp, data->seqType) == data->qFrame && mspGetRefStrand(msp) == data->qStrand)
	{
	  drawExonBoundary(msp, data);
	}
    }
}


/* Draw a vertical yellow line to indicate an insertion between two bases in the match sequence */
static void drawInsertionMarker(int sIdx, 
				int lastFoundSIdx, 
				int qIdx, 
				int lastFoundQIdx, 
				int x, 
				int y, 
				RenderData *data)
{
  if (sIdx != UNSET_INT && lastFoundSIdx != UNSET_INT && abs(sIdx - lastFoundSIdx) > 1)
    {
      /* There is a gap between this index and the previous one (in a left-to-right sense),
       * so draw an insertion marker between these two bases (at the position given by x and y). */
      
      /* This is not very sophisticated - just uses a fudge factor to find a suitable width and
       * draws it half over the current base and half over the previous one. */
      int gapWidth = roundNearest((gdouble)(data->charWidth) * GAP_WIDTH_AS_FRACTION);

      if (gapWidth < MIN_GAP_WIDTH)
	{
	  gapWidth = MIN_GAP_WIDTH;
        }

      /* See if either of the bases is selected */
      if ((qIdx == data->selectedBaseIdx && qIdx != UNSET_INT) ||
	  (lastFoundQIdx == data->selectedBaseIdx && lastFoundQIdx != UNSET_INT))
	{
	  /* Draw the two halves separately, so we can colour them in the selected/unselected 
	   * colour as necessary to match the base they are lying over. */
	  GdkColor *colour1 = (lastFoundQIdx == data->selectedBaseIdx) ? data->insertionColourSelected : data->insertionColour;
	  gdk_gc_set_foreground(data->gc, colour1);
	  const int offset1 = gapWidth / 2;
	  drawRectangle2(data->window, data->drawable, data->gc, TRUE, x - offset1, y, offset1, data->charHeight);

	  GdkColor *colour2 = (qIdx == data->selectedBaseIdx) ? data->insertionColourSelected : data->insertionColour;
	  gdk_gc_set_foreground(data->gc, colour2);
	  const int offset2 = gapWidth - offset1; /* may not be the same as offset1 if gapWidth is an odd number */
	  drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, offset2, data->charHeight);
	}
      else
	{
	  /* No selections to worry about. Draw the whole thing in the normal colour */
	  gdk_gc_set_foreground(data->gc, data->insertionColour);
	  drawRectangle2(data->window, data->drawable, data->gc, TRUE, x - gapWidth/2, y, gapWidth, data->charHeight);
	}
    }
}


static void drawSequenceText(GtkWidget *tree,
			     gchar *displayText, 
			     const IntRange const *segmentRange,
			     RenderData *data)
{
  if (g_utf8_validate(displayText, -1, NULL))
    {
      /* Get the coords for the first base. The display text should have been
       * was constructed such that everything else will line up from here. */
      int x, y;
      getCoordsForBaseIdx(0, segmentRange, data, &x, &y);
      
      PangoFontDescription *font_desc = pango_font_description_copy(tree->style->font_desc);
      PangoLayout *layout = getLayoutFromText(displayText, tree, font_desc);
      pango_font_description_free(font_desc);

      if (layout)
	{
	  paintLayout2(tree->style, data->window, data->drawable, data->state, TRUE, NULL, tree, NULL, x, y, layout);
	  g_object_unref(layout);
	}
      else
	{
	  messerror("Error creating layout while trying to display sequence:\n%s\n", displayText);
	}
    }
  else
    {
      messerror("Invalid string constructed when trying to display sequence:\n%s\n", displayText);
    }
  
}


/* Get the range of the msp that is inside the currently displayed range. The returned
 * range has UNSET_INTs if no part of the msp is visible. */
static IntRange getVisibleMspRange(MSP *msp, RenderData *data)
{
  IntRange result = {UNSET_INT, UNSET_INT};
  
  /* Find the start/end of the MSP in terms of the display coords */
  const int coord1 = convertDnaToPeptide(msp->qstart, data->qFrame, data->numFrames, data->rightToLeft, data->refSeqRange, NULL);
  const int coord2 = convertDnaToPeptide(msp->qend, data->qFrame, data->numFrames, data->rightToLeft, data->refSeqRange, NULL);
  
  int minIdx = min(coord1, coord2);
  int maxIdx = max(coord1, coord2);

  /* Check it's in the visible range. */
  if (maxIdx >= data->displayRange->min && minIdx <= data->displayRange->max)
    {
      /* Limit the range of indices in the match sequence to those within the display range. */

      if (minIdx < data->displayRange->min)
	{
	  minIdx = data->displayRange->min;
	}
      
      if (maxIdx > data->displayRange->max)
	{
	  maxIdx = data->displayRange->max;
	}
      
      result.min = minIdx;
      result.max = maxIdx;
    }

  return result;
}


static void drawDnaSequence(SequenceCellRenderer *renderer,
			    MSP *msp,
			    GtkWidget *tree,
			    RenderData *data)
{
  /* Extract the section of the reference sequence that we're interested in. */
  IntRange segmentRange = getVisibleMspRange(msp, data);
  
  /* Nothing to do if this msp is not in the visible range */
  if (segmentRange.min == UNSET_INT)
    {
      return;
    }
  
  gchar *refSeqSegment = getSequenceSegment(data->mainWindow,
					    mainWindowGetRefSeq(data->mainWindow),
					    mainWindowGetRefSeqRange(data->mainWindow),
					    segmentRange.min, 
					    segmentRange.max, 
					    data->qStrand, 
					    data->seqType,
					    data->qFrame, 
					    data->numFrames,
					    data->rightToLeft,
					    data->rightToLeft,
					    TRUE);

  if (refSeqSegment)
    {
      /* We'll populate a string with the characters we want to display as we loop through the indices. */
      const int segmentLen = segmentRange.max - segmentRange.min + 1;
      gchar displayText[segmentLen + 1];
      
      int lastFoundSIdx = UNSET_INT;  /* remember the last index where we found a valid base */
      int lastFoundQIdx = UNSET_INT;  /* remember the last index where we found a valid base */

      int segmentIdx = 0;
      for ( ; segmentIdx < segmentLen; ++segmentIdx)
	{
	  int x, y;
	  getCoordsForBaseIdx(segmentIdx, &segmentRange, data, &x, &y);
	  
	  /* Find the base in the match sequence and draw the background colour according to how well it matches */
	  int sIdx = UNSET_INT, qIdx = UNSET_INT;
	  drawBase(msp, segmentIdx, &segmentRange, refSeqSegment, data, x, y, displayText, &sIdx, &qIdx);
	  
	  /* If there is an insertion (i.e. extra bases on the match sequence) between this 
	   * and the previous coord, draw a marker */
	  drawInsertionMarker(sIdx, lastFoundSIdx, qIdx, lastFoundQIdx, x, y, data);

	  if (sIdx != UNSET_INT)
	    {
	      lastFoundSIdx = sIdx;
	      lastFoundQIdx = qIdx;
	    }
	}

      /* Null-terminate the string */
      insertChar(displayText, &segmentIdx, '\0', msp);

      /* Draw the sequence text */
      drawSequenceText(tree, displayText, &segmentRange, data);
      
      g_free(refSeqSegment);
    }
}    


/* There can be multiple MSPs in the same cell. This function loops through them
 * and draws each one. */
static void drawMsps(SequenceCellRenderer *renderer,
		     GtkWidget *tree,
		     GdkWindow *window, 
		     GtkStateType state,
		     GdkRectangle *cell_area)
{
  /* Extract all the info from the tree that we'll need repeatedly. */
  TreeProperties *treeProperties = treeGetProperties(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(treeProperties->detailView);
  MainWindowProperties *mainWindowProperties = mainWindowGetProperties(detailViewProperties->mainWindow);
  
  const gboolean highlightDiffs = detailViewProperties->highlightDiffs; /* swap match/mismatch colours if this is true */
  
  RenderData data = {
    cell_area,
    window,
    state,
    gdk_gc_new(window),
    mainWindowProperties->strandsToggled,
    treeGetStrand(tree),
    treeProperties->readingFrame,
    detailViewProperties->selectedBaseIdx,
    detailViewProperties->cellXPadding,
    detailViewProperties->cellYPadding,
    mainWindowProperties->numReadingFrames,
    detailViewProperties->charWidth,
    detailViewProperties->charHeight,
    mainWindowProperties->seqType,
    mainWindowProperties->blastMode,
    &detailViewProperties->displayRange,
    &mainWindowProperties->refSeqRange,
    highlightDiffs,
    widgetGetDrawable(tree),
    detailViewProperties->mainWindow,
    &detailViewProperties->exonColour,
    &detailViewProperties->exonColourSelected,
    &detailViewProperties->insertionColour,
    &detailViewProperties->insertionColourSelected,
    highlightDiffs ? &detailViewProperties->mismatchColour : &detailViewProperties->matchColour,
    highlightDiffs ? &detailViewProperties->mismatchColourSelected : &detailViewProperties->matchColourSelected,
    &detailViewProperties->consColour,
    &detailViewProperties->consColourSelected,
    highlightDiffs ? &detailViewProperties->matchColour : &detailViewProperties->mismatchColour,
    highlightDiffs ? &detailViewProperties->matchColourSelected : &detailViewProperties->mismatchColourSelected,
    &detailViewProperties->exonBoundaryColourStart,
    &detailViewProperties->exonBoundaryColourEnd,
    detailViewProperties->exonBoundaryLineWidth,
    detailViewProperties->exonBoundaryLineStyleStart,
    detailViewProperties->exonBoundaryLineStyleEnd
  };  
  
  /* If a base is selected highlight it now in case we don't come to process it (in 
   * which case this will get drawn over). */
  highlightSelectedBase(data.selectedBaseIdx, &detailViewProperties->mismatchColourSelected, &data);

  /* Draw all MSPs in this row */
  GList *mspListItem = renderer->mspGList;
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);
      
      if (mspIsExon(msp))
	{
	  drawExon(renderer, msp, tree, &data);
	}
      else if (mspIsBlastMatch(msp))
	{
	  drawDnaSequence(renderer, msp, tree, &data);
	}
    }
  
  drawVisibleExonBoundaries(tree, &data);
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


/* This function checks whether the cell is in a group that should be higlighted
 * and, if so, sets the cell background colour accordingly. */
static void setBackgroundColour(GtkCellRenderer *cell, GtkWidget *tree, GdkWindow *window, GdkRectangle *background_area)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  
  if (renderer->data)
    {
      /* Find out whether the MSP(s) that this cell is displaying are in 
       * a grouped sequence. */
      MSP *msp = (MSP*)(renderer->data->data);
      GtkWidget *mainWindow = treeGetMainWindow(tree);
      SequenceGroup *group = mainWindowGetSequenceGroup(mainWindow, msp->sname);
      
      if (group && group->highlighted)
	{
	  GdkGC *gc = gdk_gc_new(window);
	  gdk_gc_set_foreground(gc, &group->highlightColour);
	  drawRectangle2(window, widgetGetDrawable(tree), gc, TRUE, background_area->x, background_area->y, background_area->width, background_area->height);
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
  get_size(cell, widget, cell_area, x_offset, y_offset, width, height);  
}


/***************************************************************************
 *
 *  sequence_cell_renderer_render: crucial - do the rendering.
 *
 ***************************************************************************/

static void
sequence_cell_renderer_render (GtkCellRenderer *cell,
			       GdkWindow       *window,
			       GtkWidget       *tree,
			       GdkRectangle    *background_area,
			       GdkRectangle    *cell_area,
			       GdkRectangle    *expose_area,
			       guint            flags)
{
  setBackgroundColour(cell, tree, window, background_area);
  
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  
  MSP *msp = renderer->mspGList ? (MSP*)(renderer->mspGList->data) : NULL;
  
  if (msp)
    {
      drawMsps(renderer, tree, window, getState(tree, flags), cell_area);
    }
  else
    {
      drawText(renderer, tree, window, getState(tree, flags), cell_area);
    }
}



