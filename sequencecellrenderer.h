/*
 *  SequenceCellRenderer.h
 *  GtkTest
 *
 *  Created by Gemma Barson on 15/10/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _sequence_cell_renderer_included_
#define _sequence_cell_renderer_included_

#include <SeqTools/blxview.h>
#include <gtk/gtk.h>
#include <gtk/gtkcellrenderertext.h>

/* Some boilerplate GObject type check and type cast macros.
 *  'klass' is used here instead of 'class', because 'class'
 *  is a c++ keyword */

#define SEQUENCE_CELL_RENDERER_TYPE             (sequence_cell_renderer_get_type())
#define SEQUENCE_CELL_RENDERER(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj),  SEQUENCE_CELL_RENDERER_TYPE, SequenceCellRenderer))
#define SEQUENCE_CELL_RENDERER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass),  SEQUENCE_CELL_RENDERER_TYPE, SequenceCellRendererClass))
#define IS_SEQUENCE_CELL_RENDERER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), SEQUENCE_CELL_RENDERER_TYPE))
#define IS_SEQUENCE_CELL_RENDERER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass),  SEQUENCE_CELL_RENDERER_TYPE))
#define SEQUENCE_CELL_RENDERER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj),  SEQUENCE_CELL_RENDERER_TYPE, SequenceCellRendererClass))

/* Some custom constants */
#define UNSET_INT  -1

/* SequenceCellRenderer: Our custom cell renderer
 *   structure. Extend according to need */

typedef struct _SequenceCellRenderer
{
  GtkCellRenderer   parent;
  
  /* The cell renderer can be used to render a match sequence or plain text */
  MSP *msp;
  char *text;
  gboolean useMsp;
  
  /* Store a pointer to the detail view tree that this cell renderer renders mach sequences for */
  GtkWidget *tree;
} SequenceCellRenderer;


typedef struct _SequenceCellRendererClass
{
  GtkCellRendererTextClass  parent_class;
} SequenceCellRendererClass;

GType                sequence_cell_renderer_get_type (void);

GtkCellRenderer     *sequence_cell_renderer_new (void);

#endif /* _sequence_cell_renderer_included_ */
