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
#include <SeqTools/blixem_.h>
#include <SeqTools/sequencecellrenderer.h>

#define DETAIL_VIEW_TREE_NAME		  "DetailViewTreeName"
#define DETAIL_VIEW_TREE_CONTAINER_NAME	  "DetailViewTreeContainerName"


/* This struct holds info about a tree header widget. */
typedef struct _TreeColumnHeaderInfo
  {
    GtkWidget *headerWidget;	/* the actual header */
    GList *columnIds;		/* a list of columns spanned by this header. columns must be next to each other */
    GtkCallback refreshFunc;   /* function to be called when a refresh is requested */
  } TreeColumnHeaderInfo;



typedef struct _TreeProperties 
  {
    GtkWidget *grid;         /* The grid that this tree corresponds to */
    GtkWidget *detailView;   /* The detail view that this tree belongs to */
    GtkCellRenderer *renderer; /* The custom cell renderer to render this tree's match sequences */
    GtkWidget *sequenceColHeader; /* The custom header for the sequence column, or NULL if N/A */
    
    int readingFrame;	      /* Which reading frame this tree displays */
    GList *treeColumnHeaderList;  /* List of info about the tree column headers */
  } TreeProperties;
GdkColor exonColour;


/* Public function declarations */
TreeProperties*	  treeGetProperties(GtkWidget *widget);
GtkAdjustment*	  treeGetAdjustment(GtkWidget *tree);
GtkCellRenderer*  treeGetRenderer(GtkWidget *tree);
int		  treeGetCharWidth(GtkWidget *tree);
GtkWidget*	  treeGetGrid(GtkWidget *tree);
Strand		  treeGetStrand(GtkWidget *tree);
gboolean	  treeGetStrandsToggled(GtkWidget *tree);
int		  treeGetNumReadingFrames(GtkWidget *tree);
int		  treeGetSelectedBaseIdx(GtkWidget *tree);
char*		  treeGetRefSeq(GtkWidget *tree);
IntRange*	  treeGetDisplayRange(GtkWidget *tree);
MSP*		  treeGetMsp(GtkTreeModel *model, GtkTreeIter *iter);
GtkTreeModel*	  treeGetVisibleDataModel(GtkTreeView *tree);
GtkTreeModel*	  treeGetBaseDataModel(GtkTreeView *tree);
GtkWidget*	  treeGetMainWindow(GtkWidget *tree);
//gboolean	  treePathIsSelected(GtkTreeView *tree, GtkTreePath *path, GtkTreeModel *model);
GtkWidget*	  treeGetDetailView(GtkWidget *tree);

void		  callFuncOnAllDetailViewTrees(GtkWidget *widget, gpointer data);

void		  deselectAllSiblingTrees(GtkWidget *tree, gboolean includeCurrent);
void		  deselectAllRows(GtkWidget *tree, gpointer data);

void		  treeSortByName(GtkWidget *tree, gpointer data);
void		  treeSortById(GtkWidget *tree, gpointer data);
void		  treeSortByScore(GtkWidget *tree, gpointer data);
void		  treeSortByPos(GtkWidget *tree, gpointer data);
void		  refilterTree(GtkWidget *tree, gpointer data);
void		  refreshTreeAndGrid(GtkWidget *tree, gpointer data);
void		  refreshTreeHeaders(GtkWidget *tree, gpointer data);
void		  treeUpdateFontSize(GtkWidget *tree, gpointer data);

void		  selectRowsForSelectedMsps(GtkWidget *tree, gpointer data);

void		  addMspToTreeModel(GtkTreeModel *model, 
				    MSP *msp,
				    GtkWidget *tree);

GtkWidget*	  createDetailViewTree(GtkWidget *grid, 
				       GtkWidget *detailView,
				       GtkCellRenderer *renderer,
				       GList **treeList,
				       GList *columnList,
				       BlxSeqType seqType,
				       const char const *refSeqName,
				       const int frame);

void		   treeCreateBaseDataModel(GtkWidget *tree, gpointer data);
void		   treeCreateFilteredDataModel(GtkWidget *tree, gpointer data);

#endif /* _detail_view_tree_included_ */