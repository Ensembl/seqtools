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
#include <SeqTools/blxwindow.h>
#include <SeqTools/utilities.h>
#include <gtk/gtkcellrenderertext.h>
#include <stdlib.h>


#define SEQUENCE_CELL_RENDERER_NAME	"SequenceCellRenderer"
#define GAP_WIDTH_AS_FRACTION		0.375	/* multiplier used to get the width of the "gap" marker based on a fraction of char width */
#define MIN_GAP_WIDTH			2

typedef struct _RenderData
  {
    BlxViewContext *bc;
    GdkRectangle *cell_area;
    GdkWindow *window;
    GtkStateType state;
    GdkGC *gc;
    const BlxStrand qStrand;
    const int qFrame;
    const int selectedBaseIdx;
    const gboolean seqSelected;
    const int cellXPadding;
    const int cellYPadding;
    const int charWidth;
    const int charHeight;
    const IntRange const *displayRange;
    const gboolean highlightDiffs;
    GdkDrawable *drawable;
    GtkWidget *blxWindow;
    GdkColor *exonColorCds;
    GdkColor *exonColorCdsSelected;
    GdkColor *exonColorUtr;
    GdkColor *exonColorUtrSelected;
    GdkColor *insertionColor;
    GdkColor *insertionColorSelected;
    GdkColor *matchColor;
    GdkColor *matchColorSelected;
    GdkColor *consColor;
    GdkColor *consColorSelected;
    GdkColor *mismatchColor;
    GdkColor *mismatchColorSelected;
    GdkColor *unalignedSeqColor;
    GdkColor *unalignedSeqColorSelected;
    GdkColor *exonBoundaryColorStart;
    GdkColor *exonBoundaryColorEnd;
    GdkColor *polyAColor;
    GdkColor *polyAColorSelected;
    int exonBoundaryWidth;
    GdkLineStyle exonBoundaryStyleStart;
    GdkLineStyle exonBoundaryStyleEnd;
    gboolean showUnalignedSeq;
    gboolean limitUnalignedBases;
    int numUnalignedBases;
  } RenderData;



/* Properties */
enum 
{
  PROP_0,
  
  PROP_TEXT,
  PROP_MSP,
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

/* Custom version of gdk_draw_line that draws the same thing to two GdkDrawables */
void drawLine2(GdkDrawable *drawable1,
               GdkDrawable *drawable2,
               GdkGC *gc,
               gint x1,
               gint y1,
               gint x2,
               gint y2)
{
  if (drawable1)
    gdk_draw_line(drawable1, gc, x1, y1, x2, y2);
  
  if (drawable2)
    gdk_draw_line(drawable2, gc, x1, y1, x2, y2);
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
  cellrenderersequence->text = NULL;
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
                                   g_param_spec_pointer (RENDERER_DATA_PROPERTY,
                                                         "Data",
                                                         "Pointer to the msp",
                                                         G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_MSP,
                                   g_param_spec_pointer (RENDERER_SEQUENCE_PROPERTY,
                                                         "Sequence",
                                                         "Pointer to an msp whose sequence to display",
                                                         G_PARAM_WRITABLE));

  g_object_class_install_property (object_class,
                                   PROP_TEXT,
                                   g_param_spec_string (RENDERER_TEXT_PROPERTY,
                                                        "Text",
                                                        "Text field",
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
  renderer->text = NULL;
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

    case PROP_TEXT:
      setAllPropertiesNull(renderer);
      renderer->text = g_strdup(g_value_get_string(value));
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
      
    case PROP_TEXT:
      g_value_set_string(value, renderer->text);
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


/* Render function for a cell that contains simple text */
static void drawText(SequenceCellRenderer *renderer, 
		     GtkWidget *tree,
		     GdkWindow *window, 
		     GtkStateType state, 
		     GdkRectangle *cell_area)
{
  gchar *displayText = renderer->text;
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
 * selected and, if so, colors the background for that base with the given color. */
static void highlightSelectedBase(const int selectedBaseIdx,
				  GdkColor *highlightColor,
				  RenderData *data)
{
  if (selectedBaseIdx != UNSET_INT && valueWithinRange(selectedBaseIdx, data->displayRange))
    {
      /* Convert the display-range index to a 0-based index for the section of sequence displayed */
      const int segmentIdx = selectedBaseIdx - data->displayRange->min;
      
      int x, y;
      getCoordsForBaseIdx(segmentIdx, data->displayRange, data, &x, &y);

      gdk_gc_set_foreground(data->gc, highlightColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, data->charWidth, data->charHeight);
    }
}


/* Utility to get the exon color based on whether it is CDS/UTR and whether it is selected */
static GdkColor* getExonFillColor(const MSP const *msp, const gboolean isSelected, RenderData *data)
{
  if (msp->type == BLXMSP_EXON_CDS)
    return (isSelected ? data->exonColorCdsSelected : data->exonColorCds);
  else /* UTR and undefined exon types */
    return (isSelected ? data->exonColorUtrSelected : data->exonColorUtr);
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

  /* Just draw one big rectangle the same color for the whole thing. Color depends if row selected. */
  GdkColor *color = getExonFillColor(msp, data->seqSelected, data);
  gdk_gc_set_foreground(data->gc, color);
  drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, width, data->charHeight);
  
  /* If a base is selected, highlight it. Its color depends on whether it the base is within the exon range or not. */
  if (data->selectedBaseIdx != UNSET_INT && valueWithinRange(data->selectedBaseIdx, &segmentRange))
    {
      /* Negate the color if double-selected (i.e. if the row is selected as well) */
      GdkColor *color = getExonFillColor(msp, !data->seqSelected, data);
      highlightSelectedBase(data->selectedBaseIdx, color, data);
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


/* Get the base from the given sequence at the given index. Converts to
 * upper or lower case as appropriate for the sequence type. The given index is
 * assumed to be 1-based. If the given BlxSequence has null sequence data, then
 * this function returns a padding character instead */
static char getMatchSeqBase(BlxSequence *blxSeq, const int sIdx, const BlxSeqType seqType)
{
  char result = SEQUENCE_CHAR_PAD;
  
  if (blxSeq && blxSeq->sequence && blxSeq->sequence->str && sIdx <= blxSeq->sequence->len)
    {
      result = blxSeq->sequence->str[sIdx - 1];
      result = convertBaseToCorrectCase(result, seqType);
    }

  return result;
}


/* Color in the background of a particular base in the given match sequence. Returns
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
  GdkColor *baseBgColor = NULL;
  
  /* From the segment index, find the display index and the ref seq index */
  const int displayIdx = segmentRange->min + segmentIdx;
  *qIdx = convertDisplayIdxToDnaIdx(displayIdx, data->bc->seqType, data->qFrame, 1, data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange);
  
  *sIdx = gapCoord(msp, *qIdx, data->bc->numFrames, data->qStrand, data->bc->displayRev, 
		   data->showUnalignedSeq, data->limitUnalignedBases, data->numUnalignedBases);
  
  /* Highlight the base if its base index is selected, or if its sequence is selected.
   * (If it is selected in both, show it in the normal color) */
  gboolean selected = (displayIdx == data->selectedBaseIdx) != data->seqSelected;

  if (*qIdx < msp->qRange.min || *qIdx > msp->qRange.max)
    {
      /* We're outside the alignment range. There might still be a base to display if
       * we're displaying unaligned parts of the match sequence; otherwise, we show nothing. */
      if (*sIdx != UNSET_INT)
	{
	  sBase = getMatchSeqBase(msp->sSequence, *sIdx, data->bc->seqType);
	  baseBgColor = selected ? data->unalignedSeqColorSelected : data->unalignedSeqColor;
	}
    }
  else if (*sIdx == UNSET_INT)
    {
      /* We're inside the alignment range but there is no base to display: we must be in a gap. */
      sBase = SEQUENCE_CHAR_GAP;
      baseBgColor = selected ? data->mismatchColorSelected : data->mismatchColor;
    }
  else if (mspIsPolyATail(msp))
    {
      /* Special treatment for bases inside a polyA tail: only show polyA tails if the sequence
       * is selected. */
      if (data->seqSelected)
	{
	  baseBgColor = selected ? data->polyAColorSelected : data->polyAColor;
	  sBase = getMatchSeqBase(msp->sSequence, *sIdx, data->bc->seqType);
	}
    }
  else
    {
      /* There is a base in the match sequence. See if it matches the ref sequence */
      sBase = getMatchSeqBase(msp->sSequence, *sIdx, data->bc->seqType);
      char qBase = refSeqSegment[segmentIdx];

      if (sBase == qBase)
	{
	  /* Match */
	  baseBgColor = selected ? data->matchColorSelected : data->matchColor;
	  
	  /* If we're highlighting differences, don't show this base (show a dash instead) */
	  if (data->highlightDiffs)
	    {
	      sBase = SEQUENCE_CHAR_BLANK;
	    }
	}
      else if (data->bc->blastMode != BLXMODE_BLASTN && PAM120[aa_atob[(unsigned int)qBase]-1 ][aa_atob[(unsigned int)sBase]-1 ] > 0)
	{
	  /* 'Conserved' match (i.e. similar amino acid) */
	  baseBgColor = selected ? data->consColorSelected : data->consColor;
	}
      else
	{
	  /* Mismatch */
	  baseBgColor = selected ? data->mismatchColorSelected : data->mismatchColor;
	}
    }

  /* Draw the background color */
  if (baseBgColor)
    {
      gdk_gc_set_foreground(data->gc, baseBgColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x, y, data->charWidth, data->charHeight);
    }
  
  if (sBase != '\0')
    {
      /* Add this character into the display text */
      displayText[segmentIdx] = sBase;
    }
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
      /* Get the msp's start/end in terms of the display coords */
      const int coord1 = convertDnaIdxToDisplayIdx(msp->qRange.min, rd->bc->seqType, rd->qFrame, rd->bc->numFrames, rd->bc->displayRev, &rd->bc->refSeqRange, NULL);
      const int coord2 = convertDnaIdxToDisplayIdx(msp->qRange.max, rd->bc->seqType, rd->qFrame, rd->bc->numFrames, rd->bc->displayRev, &rd->bc->refSeqRange, NULL);
      const int minIdx = min(coord1, coord2);
      const int maxIdx = max(coord1, coord2);
      
      if (minIdx >= rd->displayRange->min && minIdx <= rd->displayRange->max)
	{
	  /* Draw the lower index. The color and line style depend on whether it's the start or end index. */
	  GdkColor *color = rd->bc->displayRev ? rd->exonBoundaryColorEnd : rd->exonBoundaryColorStart;
	  gdk_gc_set_foreground(rd->gc, color);
	  
	  GdkLineStyle lineStyle = rd->bc->displayRev ? rd->exonBoundaryStyleEnd : rd->exonBoundaryStyleStart;
	  gdk_gc_set_line_attributes(rd->gc, rd->exonBoundaryWidth, lineStyle, GDK_CAP_BUTT, GDK_JOIN_MITER);

	  const int idx = minIdx - rd->displayRange->min;

	  int x, y;
	  getCoordsForBaseIdx(idx, rd->displayRange, rd, &x, &y);
	  
          drawLine2(rd->window, rd->drawable, rd->gc, x, y, x, y + rd->charHeight);
	}
      
      if (maxIdx >= rd->displayRange->min && maxIdx <= rd->displayRange->max)
	{
	  /* Draw the upper index. The color and line style depend on whether it's the start or end index. */
	  GdkColor *color = rd->bc->displayRev ? rd->exonBoundaryColorStart : rd->exonBoundaryColorEnd;
	  gdk_gc_set_foreground(rd->gc, color);
	  
	  GdkLineStyle lineStyle = rd->bc->displayRev ? rd->exonBoundaryStyleStart : rd->exonBoundaryStyleEnd;
	  gdk_gc_set_line_attributes(rd->gc, rd->exonBoundaryWidth, lineStyle, GDK_CAP_BUTT, GDK_JOIN_MITER);
	  
	  const int idx = maxIdx + 1 - rd->displayRange->min;

	  int x, y;
	  getCoordsForBaseIdx(idx, rd->displayRange, rd, &x, &y);

          drawLine2(rd->window, rd->drawable, rd->gc, x, y, x, y + rd->charHeight);
	}
    }
  
  return FALSE;
}


/* Draw the boundaries of all exons in the given tree that are within the current
 * display range */
void drawVisibleExonBoundaries(GtkWidget *tree, RenderData *data)
{
  /* Loop through all MSPs. */
  const MSP *msp = blxWindowGetMspList(data->blxWindow);

  for ( ; msp; msp = msp->next)
    {
      if (mspIsExon(msp) && mspGetRefFrame(msp, data->bc->seqType) == data->qFrame && mspGetRefStrand(msp) == data->qStrand)
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

      /* No selections to worry about. Draw the whole thing in the normal color */
      gdk_gc_set_foreground(data->gc, data->insertionColor);
      drawRectangle2(data->window, data->drawable, data->gc, TRUE, x - gapWidth/2, y, gapWidth, data->charHeight);
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
	  g_warning("Error creating layout while trying to display sequence:\n%s\n", displayText);
	}
    }
  else
    {
      g_critical("Invalid string constructed when trying to display sequence.\n");
    }
  
}


/* Get the range of the msp that is inside the currently displayed range. The returned
 * range has UNSET_INTs if no part of the msp is visible. */
static IntRange getVisibleMspRange(MSP *msp, RenderData *data)
{
  IntRange result = {UNSET_INT, UNSET_INT};
  const gboolean sameDirection = mspGetRefStrand(msp) == mspGetMatchStrand(msp);
  
  /* Find the start/end of the MSP in terms of the display coords */
  const int idx1 = convertDnaIdxToDisplayIdx(msp->qRange.min, data->bc->seqType, data->qFrame, data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange, NULL);
  const int idx2 = convertDnaIdxToDisplayIdx(msp->qRange.max, data->bc->seqType, data->qFrame, data->bc->numFrames, data->bc->displayRev, &data->bc->refSeqRange, NULL);
  
  /* Get the min/max coord */
  int minIdx = min(idx1, idx2);
  int maxIdx = max(idx1, idx2);

  /* If we're displaying additional bases from the unaligned part of the match sequence, add those 
   * bases here (not applicable to exons) */
  if (data->showUnalignedSeq && mspIsBlastMatch(msp))
    {
      const int seqLen = mspGetMatchSeqLen(msp);
      
      /* Find the offset to the start/end of the sequence */
      int startOffset = msp->sRange.min - 1;
      int endOffset = seqLen - msp->sRange.max;
      
      if (data->limitUnalignedBases)
        {
          /* Limit the offset to the given number of bases */
          startOffset = min(startOffset, data->numUnalignedBases);
          endOffset = min(endOffset, data->numUnalignedBases);
        }

      if (sameDirection != data->bc->displayRev)
	{
	  minIdx -= startOffset;
	  maxIdx += endOffset;
	}
      else
	{
	  maxIdx += startOffset;
	  minIdx -= endOffset;
	}
    }
  
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
  
  GError *error = NULL;
  gchar *refSeqSegment = getSequenceSegment(data->bc,
					    data->bc->refSeq,
					    segmentRange.min, 
					    segmentRange.max, 
					    data->qStrand, 
					    data->bc->seqType,
					    data->qFrame, 
					    data->bc->displayRev,
					    data->bc->displayRev,
					    TRUE,
					    TRUE,
					    &error);

  if (!refSeqSegment)
    {
      g_assert(error);
      prefixError(error, "Could not draw alignment for sequence '%s'. ", msp->sname);
      reportAndClearIfError(&error, G_LOG_LEVEL_WARNING);
      return;
    }
    
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
      
      /* Find the base in the match sequence and draw the background color according to how well it matches */
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


/* There can be multiple MSPs in the same cell. This function loops through them
 * and draws each one (IF it is in the correct strand/frame for this tree). */
static void drawMsps(SequenceCellRenderer *renderer,
		     GtkWidget *tree,
		     GdkWindow *window, 
		     GtkStateType state,
		     GdkRectangle *cell_area)
{
  /* Extract all the info from the tree that we'll need repeatedly. */
  TreeProperties *treeProperties = treeGetProperties(tree);
  DetailViewProperties *detailViewProperties = detailViewGetProperties(treeProperties->detailView);
  BlxViewContext *bc = blxWindowGetContext(detailViewProperties->blxWindow);
  
  const gboolean highlightDiffs = blxContextGetFlag(bc, BLXFLAG_HIGHLIGHT_DIFFS); /* swap match/mismatch colors if this is true */
  const MSP *firstMsp = (const MSP*)(renderer->mspGList->data);
  const BlxSequence *seq = firstMsp ? firstMsp->sSequence : NULL;

  GdkColor *matchColor = getGdkColor(BLXCOLOR_MATCH, bc->defaultColors, FALSE, bc->usePrintColors);
  GdkColor *matchColorSelected = getGdkColor(BLXCOLOR_MATCH, bc->defaultColors, TRUE, bc->usePrintColors);
  GdkColor *mismatchColor = getGdkColor(BLXCOLOR_MISMATCH, bc->defaultColors, FALSE, bc->usePrintColors);
  GdkColor *mismatchColorSelected = getGdkColor(BLXCOLOR_MISMATCH, bc->defaultColors, TRUE, bc->usePrintColors);
  
  RenderData data = {
    bc,
    cell_area,
    window,
    state,
    gdk_gc_new(window),
    treeGetStrand(tree),
    treeProperties->readingFrame,
    detailViewProperties->selectedBaseIdx,
    blxWindowIsSeqSelected(detailViewProperties->blxWindow, seq),
    detailViewProperties->cellXPadding,
    detailViewProperties->cellYPadding,
    detailViewProperties->charWidth,
    detailViewProperties->charHeight,
    &detailViewProperties->displayRange,
    highlightDiffs,
    widgetGetDrawable(tree),
    detailViewProperties->blxWindow,
    getGdkColor(BLXCOLOR_EXON_CDS, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_CDS, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_UTR, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_UTR, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_INSERTION, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_INSERTION, bc->defaultColors, TRUE, bc->usePrintColors),
    highlightDiffs ? mismatchColor : matchColor,
    highlightDiffs ? mismatchColorSelected : matchColorSelected,
    getGdkColor(BLXCOLOR_CONS, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_CONS, bc->defaultColors, TRUE, bc->usePrintColors),
    highlightDiffs ? matchColor : mismatchColor,
    highlightDiffs ? matchColorSelected : mismatchColorSelected,
    getGdkColor(BLXCOLOR_UNALIGNED_SEQ, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_UNALIGNED_SEQ, bc->defaultColors, TRUE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_START, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_EXON_END, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_POLYA_TAIL, bc->defaultColors, FALSE, bc->usePrintColors),
    getGdkColor(BLXCOLOR_POLYA_TAIL, bc->defaultColors, TRUE, bc->usePrintColors),
    detailViewProperties->exonBoundaryLineWidth,
    detailViewProperties->exonBoundaryLineStyleStart,
    detailViewProperties->exonBoundaryLineStyleEnd,
    blxContextGetFlag(bc, BLXFLAG_SHOW_UNALIGNED_SEQ),
    blxContextGetFlag(bc, BLXFLAG_LIMIT_UNALIGNED_BASES),
    detailViewProperties->numUnalignedBases
  };  
  
  /* If a base is selected highlight it now in case we don't come to process it (in 
   * which case this will get drawn over). */
  highlightSelectedBase(data.selectedBaseIdx, mismatchColorSelected, &data);

  /* Draw all MSPs in this row */
  GList *mspListItem = renderer->mspGList;
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);
      
      if (mspGetRefStrand(msp) == data.qStrand && mspGetRefFrame(msp, data.bc->seqType) == data.qFrame)
	{
	  if (mspIsExon(msp))
	    {
	      drawExon(renderer, msp, tree, &data);
	    }
	  else if (mspIsBlastMatch(msp) || mspIsPolyATail(msp))
	    {
	      drawDnaSequence(renderer, msp, tree, &data);
	    }
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


/* This function sets the background color for the row. */
static void setBackgroundColor(GtkCellRenderer *cell, GtkWidget *tree, GdkWindow *window, GdkRectangle *background_area)
{
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  
  if (renderer->data)
    {
      /* Find out whether the MSP(s) that this cell is displaying are in 
       * a grouped sequence or are selected. */
      MSP *msp = (MSP*)(renderer->data->data);
      GtkWidget *blxWindow = treeGetBlxWindow(tree);
      SequenceGroup *group = blxWindowGetSequenceGroup(blxWindow, msp->sSequence);
      const gboolean isSelected = blxWindowIsSeqSelected(blxWindow, msp->sSequence);
      
      GdkGC *gc = gdk_gc_new(window);

      if (isSelected || (group && group->highlighted))
	{
	  if (group && group->highlighted && isSelected)
	    {
	      /* Use the group's highlight color but darken it because it is also selected */
	      GdkColor color;
	      getSelectionColor(&group->highlightColor, &color);
	      gdk_gc_set_foreground(gc, &color);
	    }
	  else if (group && group->highlighted)
	    {
	      /* Use the group's highlight color */
	      gdk_gc_set_foreground(gc, &group->highlightColor);
	    }
	  else
	    {
	      /* Not in a group but is selected. Use the background color but darken it. */
	      GdkColor color;
	      getSelectionColor(&tree->style->bg[GTK_STATE_NORMAL], &color);
	      gdk_gc_set_foreground(gc, &color);
	    }
	}
      else
        {
  	  gdk_gc_set_foreground(gc, &tree->style->bg[GTK_STATE_NORMAL]);
        }

      drawRectangle2(window, widgetGetDrawable(tree), gc, TRUE, background_area->x, background_area->y, background_area->width, background_area->height);
    }
}


/* Draw the tree grid lines (i.e. column separators) */
static void drawGridLines(GtkWidget *tree, GdkWindow *window, GdkRectangle *cell_area)
{
  GdkDrawable *drawable = widgetGetDrawable(tree);
  BlxViewContext *bc = treeGetContext(tree);

  /* Set the line color */
  GdkGC *gc = gdk_gc_new(window);
  GdkColor *color = getGdkColor(BLXCOLOR_TREE_GRID_LINES, bc->defaultColors, FALSE, bc->usePrintColors);
  gdk_gc_set_foreground(gc, color);
  
  /* We want dashed lines */
  const int lineWidth = 1;
  gdk_gc_set_line_attributes(gc, lineWidth, GDK_LINE_ON_OFF_DASH, GDK_CAP_BUTT, GDK_JOIN_MITER);

  /* Make the dashes very short and closely packed (i.e. dash length of 1 pixel and gaps of 1 pixel) */
  int listLen = 1;
  gint8 dashList[listLen];
  dashList[0] = 1;
  gdk_gc_set_dashes(gc, 0, dashList, listLen);
  
  /* Draw vertical lines. Draw the right edge of each cell (because we don't really want a
   * line at the leftmost edge of the first cell.) */
  const int x = cell_area->x + cell_area->width - lineWidth;
  drawLine2(window, drawable, gc, x, cell_area->y, x, cell_area->y + cell_area->height);
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
  setBackgroundColor(cell, tree, window, background_area);
  
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
  
  drawGridLines(tree, window, background_area);
}



