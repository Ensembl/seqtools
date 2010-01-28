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
#define GAP_WIDTH_AS_FRACTION		0.25	/* multiplier used to get the width of the "gap" marker based on a fraction of char width */
#define MIN_GAP_WIDTH			2

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
  
//  PROP_DATA,
//  PROP_TEXT
};


/* Some boring function declarations: GObject type system stuff */
static void     sequence_cell_renderer_init       (SequenceCellRenderer      *cellrenderersequence);
static void     sequence_cell_renderer_class_init (SequenceCellRendererClass *klass);
static void     sequence_cell_renderer_finalize (GObject *gobject);

static IntRange getVisibleMspRange(MSP *msp, GtkWidget *tree);

static void drawSequenceText(GtkWidget *tree,
			     gchar *displayText, 
			     const IntRange const *segmentRange,
			     GdkRectangle *cell_area, 
			     GdkWindow *window, 
			     GtkStateType state,
			     GdkGC *gc);


static void getCoordsForBaseIdx(const int segmentIdx, 
				const IntRange const *segmentRange,
				const IntRange const *displayRange,
				const gboolean rightToLeft,
				const int charWidth,
				const int cellXPadding,
				const int cellYPadding,
				GdkRectangle *cell_area, 
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
  cellrenderersequence->msp = NULL;
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
//  g_object_class_install_property (object_class,
//                                   PROP_TEXT,
//                                   g_param_spec_string ("text",
//                                                        "Text",
//                                                        "Text to display",
//							NULL,
//                                                        G_PARAM_WRITABLE));
  
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
  renderer->msp = NULL;
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
      renderer->data = (MSP*)g_value_get_pointer(value);
      break;
      
    case PROP_MSP:
      setAllPropertiesNull(renderer);
      renderer->msp = (MSP*)g_value_get_pointer(value);
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
      g_value_set_pointer(value, renderer->msp);
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


/* Given a base index on the query sequence, find the corresonding base 
 * in subject sequence. The return value is always UNSET_INT if there is not
 * a corresponding base at this position. However, in this case the start/end
 * of the sequence (or nearest gap) is returned in the nearestIdx argument. If
 * the base exists then the return idx and nearestIdx will be the same. 
 * nearestIdx can be null if not required. */
int gapCoord(const MSP *msp, 
	     const int qIdx, 
	     const int numFrames, 
	     const Strand strand, 
	     const gboolean rightToLeft, 
	     int *nearestIdx)
{
  int result = UNSET_INT;
  
  int qSeqMin = UNSET_INT, qSeqMax = UNSET_INT, sSeqMin = UNSET_INT, sSeqMax = UNSET_INT;
  getMspRangeExtents(msp, &qSeqMin, &qSeqMax, &sSeqMin, &sSeqMax);
  
  BOOL qForward = ((strchr(msp->qframe, '+'))) ? TRUE : FALSE ;
  BOOL sForward = ((strchr(msp->sframe, '+'))) ? TRUE : FALSE ;
  BOOL sameDirection = (qForward == sForward);
  
  Array gaps = msp->gaps ;
  if (!gaps || arrayMax(gaps) < 1)
    {
      /* If strands are in the same direction, find the offset from qSeqMin and add it to sSeqMin.
       * If strands are in opposite directions, find the offset from qSeqMin and subtract it from sSeqMax. */
      int offset = (qIdx - qSeqMin)/numFrames ;
      result = (sameDirection) ? sSeqMin + offset : sSeqMax - offset ;
      
      if (result < sSeqMin)
	{
	  result = UNSET_INT;
	  
	  if (nearestIdx)
	    {
	      *nearestIdx = sSeqMin;
	    }
	}
      else if (result > sSeqMax)
	{
	  result = UNSET_INT;
	  
	  if (nearestIdx)
	    {
	      *nearestIdx = sSeqMax;
	    }
	}
      else if (nearestIdx)
	{
	  *nearestIdx = result;
	}
    }
  else
    {
      gboolean gapsReversedWrtDisplay = (rightToLeft == qForward);
      
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
		  
		  if (nearestIdx)
		    {
		      *nearestIdx = result;
		    }
		}
		else if (nearestIdx && (i == 0 || gapsReversedWrtDisplay))
		  {
		    /* Remember the start of the current range. This will be the value we
		     * return if the index is before of the first range, or if we're in a gap
		     * and gaps are ordered in the opposite direction to the display (i.e. gap
		     * coords go from low to high and display shows high to low or v.v.) */
		    *nearestIdx = curRange->s1;
		  }
	      
	      break;
	    }
	  else if (nearestIdx && !gapsReversedWrtDisplay)
	    {
	      /* Remember the end of the current range (which is the result we will return
	       * if qIdx lies after this range but before the next - unless gaps are reversed
	       * with respect to the display). */
	      *nearestIdx = curRange->s2;
	    }
	}
    }
  
  return result ;
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

  gtk_paint_layout (tree->style,
		    window,
		    state,
		    TRUE,
		    NULL,
		    tree,
		    NULL, //"cellrenderertext",
		    cell_area->x - treeGetCellXPadding(tree),
		    cell_area->y - treeGetCellYPadding(tree),
		    layout);
  
  g_object_unref(layout);
}


/* The given renderer is an MSP. This function checks if there is a base index
 * selected and, if so, colours the background for that base with the given colour. */
static void highlightSelectedBase(const int selectedBaseIdx,
				  GdkColor *highlightColour,
				  const IntRange const *displayRange,
				  const gboolean rightToLeft,
				  const int charWidth,
				  const int charHeight,
				  const int cellXPadding,
				  const int cellYPadding,
				  GdkGC *gc,
				  GdkWindow *window,
				  GdkRectangle *cell_area)
{
  if (selectedBaseIdx != UNSET_INT && indexWithinRange(selectedBaseIdx, displayRange))
    {
      /* Convert the display-range index to a 0-based index for the section of sequence displayed */
      const int segmentIdx = rightToLeft ? displayRange->max - selectedBaseIdx : selectedBaseIdx - displayRange->min;
      
      int x, y;
      getCoordsForBaseIdx(segmentIdx, displayRange, displayRange, rightToLeft, charWidth, cellXPadding, cellYPadding, cell_area, &x, &y);

      gdk_gc_set_foreground(gc, highlightColour);
      gdk_draw_rectangle(window, gc, TRUE, x, y, charWidth, charHeight);
    }
}


/* The given renderer is an MSP that is an exon. This function draws the exon
 * or part of the exon that is in view, if it is within the current display range. */
static void drawExon(SequenceCellRenderer *renderer,
		     GtkWidget *tree,
		     GdkWindow *window, 
		     GtkStateType state,
		     GdkRectangle *cell_area)
{
  GdkGC *gc = gdk_gc_new(window);
  
  IntRange segmentRange = getVisibleMspRange(renderer->msp, tree);
  const int segmentLen = segmentRange.max - segmentRange.min + 1;

  const IntRange const *displayRange = treeGetDisplayRange(tree);
  const int charWidth = treeGetCharWidth(tree);
  const int charHeight = treeGetCharHeight(tree);
  const int cellXPadding = treeGetCellXPadding(tree);
  const int cellYPadding = treeGetCellYPadding(tree);
  const gboolean rightToLeft = treeGetStrandsToggled(tree);
  
  int x, y;
  getCoordsForBaseIdx(0, &segmentRange, displayRange, rightToLeft, charWidth, cellXPadding, cellYPadding, cell_area, &x, &y);
  const int width = segmentLen * charWidth;

  /* Just draw one big rectangle the same colour for the whole thing */
  GdkColor *baseBgColour = treeGetExonColour(tree, FALSE);
  gdk_gc_set_foreground(gc, baseBgColour);
  gdk_draw_rectangle(window, gc, TRUE, x, y, width, charHeight);
  
  /* If a base is selected, highlight it. The colour depends on whether it is within our exon range
   * or not. */
  const int selectedBaseIdx = treeGetSelectedBaseIdx(tree);
  if (selectedBaseIdx != UNSET_INT)
    {
      GdkColor *bgColour = indexWithinRange(selectedBaseIdx, &segmentRange) 
	? treeGetExonColour(tree, TRUE) 
	: treeGetGapColour(tree, TRUE);
      
      highlightSelectedBase(selectedBaseIdx, bgColour, displayRange, rightToLeft, charWidth, charHeight, cellXPadding, cellYPadding, gc, window, cell_area);
    }
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
  char result = matchSeq[sIdx - 1];
  result = convertBaseToCorrectCase(result, seqType);
  return result;
}


/* Given a 0-based index into a segment of the reference sequence, find the
 * index into the full reference sequence. */
static int getRefSeqIndexFromSegment(const int segmentIdx,
				     const IntRange const *segmentRange,
				     const gboolean rightToLeft)
{
  return (rightToLeft ? segmentRange->max - segmentIdx : segmentRange->min + segmentIdx);
}


/* Given an index into the full displayed reference sequence, find the equivalent base
 * in the match sequence (converting the ref sequence coord from peptide to dna if
 * necessary). */
static int getMatchIdxFromRefIdx(MSP *msp,
				 const int refSeqIdx,
				 const int qFrame,
				 const Strand qStrand,
				 const gboolean rightToLeft,
				 const BlxSeqType seqType,
				 const int numFrames)
{
  /* If we have a peptide sequence we must convert the ref seq coord to the DNA sequence coord */
  int qIdx = refSeqIdx;

  if (seqType == BLXSEQ_PEPTIDE)
    {
      qIdx = convertPeptideToDna(qIdx, qFrame, numFrames);
    }
  
  /* Find the s index */
  return gapCoord(msp, qIdx, numFrames, qStrand, rightToLeft, NULL);
}
	   

/* Colour in the background of a particular base in the given match sequence. Returns
 * the calculated index into the match sequence (for effiency, so that we don't have to
 * recalculate it later on). Returns the equivalent index in the subject sequence, or 
 * UNSET_INT if there is none. */
static int drawBase(MSP *msp,
		    const int segmentIdx, 
		    const IntRange const *segmentRange,
		    char *refSeqSegment, 

		    const int selectedBaseIdx,
		    const int qFrame,
		    const Strand qStrand,
		    const BlxSeqType seqType,
		    const gboolean rightToLeft,
		    const int numFrames,
		    const int charWidth,
		    const int charHeight,
		    
		    GtkWidget *tree,
		    GdkWindow *window,
		    GdkGC *gc,
		    
		    const int x,
		    const int y,
		    gchar *displayText)
{
  char charToDisplay = '\0';
  GdkColor *baseBgColour;
  
  const int refSeqIdx = getRefSeqIndexFromSegment(segmentIdx, segmentRange, rightToLeft);
  const int sIdx = getMatchIdxFromRefIdx(msp, refSeqIdx, qFrame, qStrand, rightToLeft, seqType, numFrames);
  
  gboolean selected = (refSeqIdx == selectedBaseIdx);
  
  if (sIdx == UNSET_INT)
    {
      /* There is no equivalent base in the match sequence so draw a gap */
      charToDisplay = '.';
      baseBgColour = treeGetGapColour(tree, selected);
    }
  else
    {
      /* There is a base in the match sequence. See if it matches the ref sequence */
      charToDisplay = getMatchSeqBase(msp->sseq, sIdx, seqType);
      char qBase = refSeqSegment[segmentIdx];

      if (charToDisplay == qBase)
	{
	  baseBgColour = treeGetMatchColour(tree, selected); 
	}
      else
	{
	  baseBgColour = treeGetMismatchColour(tree, selected);
	}
    }

  /* Draw the background colour */
  gdk_gc_set_foreground(gc, baseBgColour);
  gdk_draw_rectangle(window, gc, TRUE, x, y, charWidth, charHeight);

  /* Add this character into the display text */
  displayText[segmentIdx] = charToDisplay;

  return sIdx;
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
				const IntRange const *displayRange,
				const gboolean rightToLeft,
				const int charWidth,
				const int cellXPadding,
				const int cellYPadding,
				GdkRectangle *cell_area, 
				int *x, 
				int* y)
{
  /* Find the start of the segment with respect to the display range */
  const int startPos = rightToLeft 
    ? displayRange->max - segmentRange->max
    : segmentRange->min - displayRange->min;
  
  /* Find the position of the character within the segment */
  int charIdx = startPos + segmentIdx;

  /* Calculate the coords */
  *x = cell_area->x - cellXPadding + (charIdx * charWidth);
  *y = cell_area->y - cellYPadding;
}


static void drawGap(int sIdx, 
		    int sIdxLastFound, 
		    int x, 
		    int y, 
		    const int charWidth,
		    const int charHeight,
		    GdkWindow *window, 
		    GdkGC *gc)
{
  if (sIdx != UNSET_INT && sIdxLastFound != UNSET_INT && abs(sIdx - sIdxLastFound) > 1)
    {
      /* There is a gap between this index and the previous one (in a left-to-right sense) */
      GdkColor col = {GDK_YELLOW};
      gdk_gc_set_foreground(gc, &col);
      
      /* This is not very sophisticated - just uses a fudge factor to find a suitable width and
       * draws it half over the current base and half over the previous one. */
      int gapWidth = round((gdouble)charWidth * GAP_WIDTH_AS_FRACTION);

      if (gapWidth < MIN_GAP_WIDTH)
	{
	  gapWidth = MIN_GAP_WIDTH;
        }
      
      gdk_draw_rectangle(window, gc, TRUE, x - gapWidth/2, y, gapWidth, charHeight);
    }
}


static void drawSequenceText(GtkWidget *tree,
			     gchar *displayText, 
			     const IntRange const *segmentRange,
			     GdkRectangle *cell_area, 
			     GdkWindow *window, 
			     GtkStateType state,
			     GdkGC *gc)
{
  if (g_utf8_validate(displayText, -1, NULL))
    {
      /* Get the coords for the first base. The display text should have been
       * was constructed such that everything else will line up from here. */
      int x, y;
      getCoordsForBaseIdx(0, segmentRange, treeGetDisplayRange(tree), treeGetStrandsToggled(tree), treeGetCharWidth(tree), treeGetCellXPadding(tree), treeGetCellYPadding(tree), cell_area, &x, &y);
      
      PangoFontDescription *font_desc = pango_font_description_copy(tree->style->font_desc);
      PangoLayout *layout = getLayoutFromText(displayText, tree, font_desc);
      pango_font_description_free(font_desc);

      if (layout)
	{
	  gtk_paint_layout (tree->style, window, state, TRUE, NULL, tree, NULL, x, y, layout);
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
static IntRange getVisibleMspRange(MSP *msp, GtkWidget *tree)
{
  IntRange result = {UNSET_INT, UNSET_INT};
  
  /* Find the start/end of the MSP in terms of the display coords */
  int minIdx = min(msp->displayStart, msp->displayEnd);
  int maxIdx = max(msp->displayStart, msp->displayEnd);

  IntRange *displayRange = treeGetDisplayRange(tree);

  /* Check it's in the visible range. */
  if (maxIdx > displayRange->min && minIdx < displayRange->max)
    {
      /* Limit the range of indices in the match sequence to those within the display range. */

      if (minIdx < displayRange->min)
	{
	  minIdx = displayRange->min;
	}
      
      if (maxIdx > displayRange->max)
	{
	  maxIdx = displayRange->max;
	}
      
      result.min = minIdx;
      result.max = maxIdx;
    }

  return result;
}


static void drawDnaSequence(SequenceCellRenderer *renderer,
			    GtkWidget *tree,
			    GdkWindow *window, 
			    GtkStateType state,
			    GdkRectangle *cell_area)
{
  GdkGC *gc = gdk_gc_new(window);

  /* Extract the section of the reference sequence that we're interested in. */
  IntRange segmentRange = getVisibleMspRange(renderer->msp, tree);
  
  /* Nothing to do if this msp is not in the visible range */
  if (segmentRange.min == UNSET_INT)
    {
      return;
    }
  
  const gboolean rightToLeft = treeGetStrandsToggled(tree);
  const Strand qStrand = treeGetStrand(tree);
  const int qFrame = treeGetFrame(tree);
  const int selectedBaseIdx = treeGetSelectedBaseIdx(tree);
  const int cellXPadding = treeGetCellXPadding(tree);
  const int cellYPadding = treeGetCellYPadding(tree);
  const int numFrames = treeGetNumReadingFrames(tree);
  const int charWidth = treeGetCharWidth(tree);
  const int charHeight = treeGetCharHeight(tree);
  const BlxSeqType seqType = treeGetSeqType(tree);
  const IntRange const *displayRange = treeGetDisplayRange(tree);
  
  GtkWidget *mainWindow = treeGetMainWindow(tree);
  
  gchar *refSeqSegment = getSequenceSegment(mainWindow, 
					    mainWindowGetRefSeq(mainWindow),
					    mainWindowGetRefSeqRange(mainWindow),
					    segmentRange.min, 
					    segmentRange.max, 
					    qStrand, 
					    seqType,
					    qFrame, 
					    numFrames,
					    rightToLeft);

  if (refSeqSegment)
    {
      /* We'll populate a string with the characters we want to display as we loop through the indices. */
      const int segmentLen = segmentRange.max - segmentRange.min + 1;
      gchar displayText[segmentLen + 1];
      
      int lastFoundSIdx = UNSET_INT;  /* remember the last index in the match sequence where we found a valid base */

      int segmentIdx = 0;
      for ( ; segmentIdx < segmentLen; ++segmentIdx)
	{
	  int x, y;
	  getCoordsForBaseIdx(segmentIdx, &segmentRange, displayRange, rightToLeft, charWidth, cellXPadding, cellYPadding, cell_area, &x, &y);
	  
	  /* Find the base in the match sequence and draw the background colour according to how well it matches */
	  int sIdx = drawBase(renderer->msp, segmentIdx, &segmentRange, refSeqSegment, selectedBaseIdx, qFrame, qStrand, seqType, rightToLeft, numFrames, charWidth, charHeight, tree, window, gc, x, y, displayText);
	  
	  /* If there is an insertion/deletion between this and the previous match, draw it now */
	  drawGap(sIdx, lastFoundSIdx, x, y, charWidth, charHeight, window, gc);

	  if (sIdx != UNSET_INT)
	    {
	      lastFoundSIdx = sIdx;
	    }
	}

      /* Null-terminate the string */
      insertChar(displayText, &segmentIdx, '\0', renderer->msp);

      /* Draw the sequence text */
      drawSequenceText(tree, displayText, &segmentRange, cell_area, window, state, gc);
      
      g_free(refSeqSegment);
    }
    
  /* If a base is selected and we've not already processed it, highlight it now */
  if (selectedBaseIdx != UNSET_INT && !indexWithinRange(selectedBaseIdx, &segmentRange))
    {
      highlightSelectedBase(selectedBaseIdx, treeGetGapColour(tree, TRUE), displayRange, rightToLeft, charWidth, charHeight, cellXPadding, cellYPadding, gc, window, cell_area);
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
  SequenceCellRenderer *renderer = SEQUENCE_CELL_RENDERER(cell);
  
  if (renderer->msp && mspIsExon(renderer->msp))
    {
      drawExon(renderer, tree, window, getState(tree, flags), cell_area);
    }
  else if (renderer->msp && treeGetSeqType(tree) == BLXSEQ_PEPTIDE)
    {
      drawDnaSequence(renderer, tree, window, getState(tree, flags), cell_area);
    }
  else if (renderer->msp)
    {
      drawDnaSequence(renderer, tree, window, getState(tree, flags), cell_area);
    }
  else
    {
      drawText(renderer, tree, window, getState(tree, flags), cell_area);
    }
}



