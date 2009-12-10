/*
 *  detailviewtree.h
 *  acedb
 *
 *  Created by Gemma Barson on 23/11/2009.
 *
 */

#ifndef _detail_view_tree_included_
#define _detail_view_tree_included_

#include <gtk/gtk.h>
#include <SeqTools/blxview.h>
#include <SeqTools/sequencecellrenderer.h>

#define DETAIL_VIEW_TREE_NAME		"DetailViewTreeName"


enum
{
  S_NAME_COL,
  SCORE_COL,
  ID_COL,
  START_COL,
  MSP_COL,
  END_COL,
  N_COLUMNS
};

typedef struct _TreeProperties 
  {
    GtkWidget *grid;         /* The grid that this tree corresponds to */
    GtkWidget *detailView;   /* The detail view that this tree belongs to */
    GtkCellRenderer *renderer; /* The custom cell renderer to render this tree's match sequences */
    
    char *refSeq;	     /* The reference sequence for this tree. (It is the complemented version if the strand is reversed.) */
  } TreeProperties;


/* Public function declarations */
TreeProperties*	  treeGetProperties(GtkWidget *widget);
GtkAdjustment*	  treeGetAdjustment(GtkWidget *tree);
GtkCellRenderer*  treeGetRenderer(GtkWidget *tree);
GtkWidget*	  treeGetGrid(GtkWidget *tree);
Strand		  treeGetStrand(GtkWidget *tree);
int		  treeGetNumReadingFrames(GtkWidget *tree);
int		  treeGetSelectedBaseIdx(GtkWidget *tree);
char*		  treeGetRefSeq(GtkWidget *tree);
IntRange*	  treeGetDisplayRange(GtkWidget *tree);
const MSP*	  treeGetMsp(GtkTreeModel *model, GtkTreeIter *iter);
GtkTreeModel*	  treeGetVisibleDataModel(GtkTreeView *tree);
GtkTreeModel*	  treeGetBaseDataModel(GtkTreeView *tree);
gboolean	  treePathIsSelected(GtkTreeView *tree, GtkTreePath *path, GtkTreeModel *model);

void		  callFuncOnAllDetailViewTrees(GtkWidget *widget, gpointer data);

void		  refilterTree(GtkWidget *tree, gpointer data);
void		  refreshTreeAndGrid(GtkWidget *tree, gpointer data);

gint		  setCellRendererFont(GtkWidget *tree, 
				      GtkCellRenderer *renderer, 
				      const char *fontName, 
				      const int fontSize);

void		  addMspToTreeModel(GtkTreeModel *model, 
				    const MSP *msp);

GtkWidget*	  createDetailViewTree(GtkWidget *grid, 
				       GtkWidget *detailView,
				       GList **treeList);

void		   treeCreateBaseDataModel(GtkWidget *tree, gpointer data);
void		   treeCreateFilteredDataModel(GtkWidget *tree, gpointer data);

#endif /* _detail_view_tree_included_ */