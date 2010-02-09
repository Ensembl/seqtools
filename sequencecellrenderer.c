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
    GdkDrawable *drawable;
    GtkWidget *mainWindow;
    GdkColor *exonColour;
    GdkColor *exonColourSelected;
    GdkColor *gapColour;
    GdkColor *gapColourSelected;
    GdkColor *matchColour;
    GdkColor *matchColourSelected;
    GdkColor *consColour;
    GdkColor *consColourSelected;
    GdkColor *mismatchColour;
    GdkColor *mismatchColourSelected;
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
  
//  PROP_DATA,
//  PROP_TEXT
};


typedef struct
  {
    GtkWidget *tree;
    GdkWindow *window;
    GdkRectangle *cell_area;
  } DrawExonBoundariesData;


/* Some boring function declarations: GObject type system stuff */
static void     sequence_cell_renderer_init       (SequenceCellRenderer      *cellrenderersequence);
static void     sequence_cell_renderer_class_init (SequenceCellRendererClass *klass);
static void     sequence_cell_renderer_finalize (GObject *gobject);

static IntRange getVisibleMspRange(MSP *msp, GtkWidget *tree);
void		drawVisibleExonBoundaries(GtkWidget *tree, GdkWindow *window, GdkRectangle *cell_area);

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



/***************************************************************************/

/* Custom version of gtk_paint_layout that draws the same thing to two GdkDrawables */
static void paintLayout2(GtkStyle *style,
			 GdkDrawable *drawable1,
			 GdkDrawable *drawable2,
			 GtkStateType state_type,
			 gboolean use_text,
			 const GdkRectangle *area,
			 GtkWidget *widget,
			 const gchar *detail,
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
				  RenderData *data,
				  GdkWindow *window,
				  GdkDrawable *drawable,
				  GdkRectangle *cell_area)
{
  if (selectedBaseIdx != UNSET_INT && valueWithinRange(selectedBaseIdx, displayRange))
    {
      /* Convert the display-range index to a 0-based index for the section of sequence displayed */
      const int segmentIdx = data->rightToLeft ? displayRange->max - selectedBaseIdx : selectedBaseIdx - displayRange->min;
      
      int x, y;
      getCoordsForBaseIdx(segmentIdx, displayRange, displayRange, data->rightToLeft, data->charWidth, data->cellXPadding, data->cellYPadding, cell_area, &x, &y);

      gdk_gc_set_foreground(data->gc, highlightColour);
      drawRectangle2(window, data->drawable, data->gc, TRUE, x, y, data->charWidth, data->charHeight);
    }
}


/* The given renderer is an MSP that is an exon. This function draws the exon
 * or part of the exon that is in view, if it is within the current display range. */
static void drawExon(SequenceCellRenderer *renderer,
		     MSP *msp,
		     GtkWidget *tree,
		     GdkWindow *window, 
		     GtkStateType state,
		     GdkRectangle *cell_area,
		     RenderData *data)
{
  IntRange segmentRange = getVisibleMspRange(msp, tree);
  const int segmentLen = segmentRange.max - segmentRange.min + 1;

  int x, y;
  getCoordsForBaseIdx(0, &segmentRange, data->displayRange, data->rightToLeft, data->charWidth, data->cellXPadding, data->cellYPadding, cell_area, &x, &y);
  const int width = segmentLen * data->charWidth;

  /* Just draw one big rectangle the same colour for the whole thing */
  gdk_gc_set_foreground(data->gc, data->exonColour);
  drawRectangle2(window, widgetGetDrawable(tree), data->gc, TRUE, x, y, width, data->charHeight);
  
  /* If a base is selected, highlight it. Its colour depends on whether it the base is within the exon range or not. */
  if (data->selectedBaseIdx != UNSET_INT && valueWithinRange(data->selectedBaseIdx, &segmentRange))
    {
      highlightSelectedBase(data->selectedBaseIdx, data->exonColourSelected, data->displayRange, data, window, widgetGetDrawable(tree), cell_area);
    }
  
  drawVisibleExonBoundaries(tree, window, cell_area);
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


/* Colour in the background of a particular base in the given match sequence. Returns
 * the calculated index into the match sequence (for effiency, so that we don't have to
 * recalculate it later on). Returns the equivalent index in the subject sequence, or 
 * UNSET_INT if there is none. */
static int drawBase(MSP *msp,
		    const int segmentIdx, 
		    const IntRange const *segmentRange,
		    char *refSeqSegment, 

		    RenderData *data,
		    GtkWidget *tree,
		    GdkWindow *window,
		    GdkDrawable *drawable,
		    
		    const int x,
		    const int y,
		    gchar *displayText)
{
  char sBase = '\0';
  GdkColor *baseBgColour;
  
  const int refSeqIdx = getRefSeqIndexFromSegment(segmentIdx, segmentRange, data->rightToLeft);
  const int sIdx = getMatchIdxFromDisplayIdx(msp, refSeqIdx, data->qFrame, data->qStrand, data->rightToLeft, data->seqType, data->numFrames);
  
  /* Highlight the base if its base index is selected */
  gboolean selected = (refSeqIdx == data->selectedBaseIdx);
  
  if (sIdx == UNSET_INT)
    {
      /* There is no equivalent base in the match sequence so draw a gap */
      sBase = '.';
      baseBgColour = selected ? data->gapColourSelected : data->gapColour;
    }
  else
    {
      /* There is a base in the match sequence. See if it matches the ref sequence */
      sBase = getMatchSeqBase(msp->sseq, sIdx, data->seqType);
      char qBase = refSeqSegment[segmentIdx];

      if (sBase == qBase)
	{
	  baseBgColour = selected ? data->matchColourSelected : data->matchColour;
	}
      else if (data->blastMode != BLXMODE_BLASTN && PAM120[aa_atob[(unsigned int)qBase]-1 ][aa_atob[(unsigned int)sBase]-1 ] > 0)
	{
	  baseBgColour = selected ? data->consColourSelected : data->consColour;
	}
      else
	{
	  baseBgColour = selected ? data->mismatchColourSelected : data->mismatchColour;
	}
    }

  /* Draw the background colour */
  gdk_gc_set_foreground(data->gc, baseBgColour);
  drawRectangle2(window, drawable, data->gc, TRUE, x, y, data->charWidth, data->charHeight);

  /* Add this character into the display text */
  displayText[segmentIdx] = sBase;

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


/* Draw the start/end boundaries for the given exon, if the start/end is within the current display range */
static gboolean drawExonBoundary(GtkTreeModel *model, GtkTreePath *path, GtkTreeIter *iter, gpointer data)
{
  /* A tree row can contain multiple MSPs. Process all of them. */
  GList *mspListItem = treeGetMsps(model, iter);
  
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      const MSP *msp = (const MSP*)(mspListItem->data);
      
      if (msp && mspIsExon(msp))
	{
	  DrawExonBoundariesData *drawData = (DrawExonBoundariesData*)data;
	  
	  GtkWidget *tree = drawData->tree;
	  const IntRange const *displayRange = treeGetDisplayRange(tree);
	  const int charWidth = treeGetCharWidth(tree);
	  const int charHeight = treeGetCharHeight(tree);
	  const int cellXPadding = treeGetCellXPadding(tree);
	  const int cellYPadding = treeGetCellYPadding(tree);
	  const gboolean rightToLeft = treeGetStrandsToggled(tree);
	  GdkGC *gc = gdk_gc_new(drawData->window);

	  const int minIdx = min(msp->displayStart, msp->displayEnd);
	  const int maxIdx = max(msp->displayStart, msp->displayEnd);
	  
	  if (minIdx >= displayRange->min && minIdx <= displayRange->max)
	    {
	      /* Draw the lower index. The colour and line style depend on whether it's the start or end index. */
	      gdk_gc_set_foreground(gc, treeGetExonBoundaryColour(tree, TRUE));
	      gdk_gc_set_line_attributes(gc, treeGetExonBoundaryWidth(tree), treeGetExonBoundaryStyle(tree, TRUE), GDK_CAP_BUTT, GDK_JOIN_MITER);

	      const int idx = rightToLeft ? displayRange->max - minIdx + 1 : minIdx - displayRange->min;

	      int x, y;
	      getCoordsForBaseIdx(idx, displayRange, displayRange, treeGetStrandsToggled(tree), charWidth, cellXPadding, cellYPadding, drawData->cell_area, &x, &y);
	      
	      gdk_draw_line(drawData->window, gc, x, y, x, y + charHeight);
	      gdk_draw_line(widgetGetDrawable(tree), gc, x, y, x, y + charHeight);
	    }
	  
	  if (maxIdx >= displayRange->min && maxIdx <= displayRange->max)
	    {
	      /* Draw the upper index. The colour and line style depend on whether it's the start or end index. */
	      gdk_gc_set_foreground(gc, treeGetExonBoundaryColour(tree, FALSE));
	      gdk_gc_set_line_attributes(gc, treeGetExonBoundaryWidth(tree), treeGetExonBoundaryStyle(tree, FALSE), GDK_CAP_BUTT, GDK_JOIN_MITER);
	      
	      const int idx = rightToLeft ? displayRange->max - maxIdx : maxIdx - displayRange->min + 1;

	      int x, y;
	      getCoordsForBaseIdx(idx, displayRange, displayRange, treeGetStrandsToggled(tree), charWidth, cellXPadding, cellYPadding, drawData->cell_area, &x, &y);

	      gdk_draw_line(drawData->window, gc, x, y, x, y + charHeight);
	      gdk_draw_line(widgetGetDrawable(tree), gc, x, y, x, y + charHeight);
	    }
	}
    }
  
  return FALSE;
}


/* Draw the boundaries of all exons within the current display range */
void drawVisibleExonBoundaries(GtkWidget *tree, GdkWindow *window, GdkRectangle *cell_area)
{
  /* Loop through all visible MSPs that are exons and add the coords to the list if they are within range. */
  GtkTreeModel *model = treeGetVisibleDataModel(GTK_TREE_VIEW(tree));
  DrawExonBoundariesData drawData = {tree, window, cell_area};
  gtk_tree_model_foreach(model, drawExonBoundary, &drawData);
}


static void drawGap(int sIdx, 
		    int sIdxLastFound, 
		    int x, 
		    int y, 
		    const int charWidth,
		    const int charHeight,
		    GdkWindow *window, 
		    GdkDrawable *drawable,
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
      
      drawRectangle2(window, drawable, gc, TRUE, x - gapWidth/2, y, gapWidth, charHeight);
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
	  paintLayout2(tree->style, window, widgetGetDrawable(tree), state, TRUE, NULL, tree, NULL, x, y, layout);
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
			    MSP *msp,
			    GtkWidget *tree,
			    GdkWindow *window, 
			    GtkStateType state,
			    GdkRectangle *cell_area,
			    RenderData *data)
{
  /* Extract the section of the reference sequence that we're interested in. */
  IntRange segmentRange = getVisibleMspRange(msp, tree);
  
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
					    TRUE);

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
	  getCoordsForBaseIdx(segmentIdx, &segmentRange, data->displayRange, data->rightToLeft, data->charWidth, data->cellXPadding, data->cellYPadding, cell_area, &x, &y);
	  
	  /* Find the base in the match sequence and draw the background colour according to how well it matches */
	  int sIdx = drawBase(msp, segmentIdx, &segmentRange, refSeqSegment, data, tree, window, data->drawable, x, y, displayText);
	  
	  /* If there is an insertion/deletion between this and the previous match, draw it now */
	  drawGap(sIdx, lastFoundSIdx, x, y, data->charWidth, data->charHeight, window, data->drawable, data->gc);

	  if (sIdx != UNSET_INT)
	    {
	      lastFoundSIdx = sIdx;
	    }
	}

      /* Null-terminate the string */
      insertChar(displayText, &segmentIdx, '\0', msp);

      /* Draw the sequence text */
      drawSequenceText(tree, displayText, &segmentRange, cell_area, window, state, data->gc);
      
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
  RenderData data = {
    gdk_gc_new(window),
    treeGetStrandsToggled(tree),
    treeGetStrand(tree),
    treeGetFrame(tree),
    treeGetSelectedBaseIdx(tree),
    treeGetCellXPadding(tree),
    treeGetCellYPadding(tree),
    treeGetNumReadingFrames(tree),
    treeGetCharWidth(tree),
    treeGetCharHeight(tree),
    treeGetSeqType(tree),
    treeGetBlastMode(tree),
    treeGetDisplayRange(tree),
    widgetGetDrawable(tree),
    treeGetMainWindow(tree),
    treeGetExonColour(tree, FALSE),
    treeGetExonColour(tree, TRUE),
    treeGetGapColour(tree, FALSE),
    treeGetGapColour(tree, TRUE),
    treeGetMatchColour(tree, FALSE),
    treeGetMatchColour(tree, TRUE),
    treeGetConsColour(tree, FALSE),
    treeGetConsColour(tree, TRUE),
    treeGetMismatchColour(tree, FALSE),
    treeGetMismatchColour(tree, TRUE)
  };  
  
  /* If a base is selected and we've not already processed it, highlight it now */
  highlightSelectedBase(treeGetSelectedBaseIdx(tree), 
			treeGetGapColour(tree, TRUE), 
			treeGetDisplayRange(tree), 
			&data,
			window, 
			widgetGetDrawable(tree), 
			cell_area);
  
  
  /* Draw all MSPs in this row */
  GList *mspListItem = renderer->mspGList;
  for ( ; mspListItem; mspListItem = mspListItem->next)
    {
      MSP *msp = (MSP*)(mspListItem->data);
      
      if (mspIsExon(msp))
	{
	  drawExon(renderer, msp, tree, window, state, cell_area, &data);
	}
      else if (mspIsBlastMatch(msp))
	{
	  drawDnaSequence(renderer, msp, tree, window, state, cell_area, &data);
	}
    }
  
  drawVisibleExonBoundaries(tree, window, cell_area);
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



